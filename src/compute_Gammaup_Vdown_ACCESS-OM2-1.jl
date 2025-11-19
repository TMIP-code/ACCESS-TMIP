# qsub -I -P y99 -q normal -l mem=190GB -l storage=scratch/gh0+scratch/y99+gdata/xp65 -l walltime=02:00:00 -l ncpus=48

using Pkg
Pkg.activate(".")
Pkg.instantiate()
const nprocs = 48

ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
using OceanTransportMatrixBuilder
using NetCDF
using YAXArrays
using DataFrames
using DimensionalData
using SparseArrays
using LinearAlgebra
using Unitful
using Unitful: s, yr, d
using Statistics
using Format
using Dates
using FileIO
using LinearSolve
import Pardiso # import Pardiso instead of using (to avoid name clash?)
using NonlinearSolve
using ProgressMeter

model = "ACCESS-OM2-1"
experiment = "1deg_jra55_iaf_omip2_cycle6"
time_window = "Jan1960-Dec1979"
@show inputdir = "/scratch/y99/TMIP/data/$model/$experiment/$time_window"

# preferred diffusivities
κVdeep = 3.0e-5 # m^2/s
κVML = 1.0      # m^2/s
κH = 300.0      # m^2/s (grid-scaling by sqrt(area))
@show κVdeep
@show κVML
@show κH
κVdeep_str = "kVdeep" * format(κVdeep, conversion = "e")
κVML_str = "kVML" * format(κVML, conversion = "e")
κH_str = "kH" * format(κH, conversion = "d")

upwind = false
@show upwind
upwind_str = upwind ? "" : "_centered"
upwind_str2 = upwind ? "upwind" : "centered"

# Load areacello and volcello for grid geometry
areacello_ds = open_dataset(joinpath(inputdir, "area_t.nc"))
dht_ds = open_dataset(joinpath(inputdir, "dht.nc")) # <- (new) cell thickness?
# TODO: caputre correlations between transport and dht
# z* coordinate varies with time in ACCESS-OM2
# volcello_ds = open_dataset(joinpath(fixedvarsinputdir, "volcello.nc")) # <- not in ACCESS-OM2: must be built from dht * area

# Load fixed variables in memory
areacello_OM2 = replace(readcubedata(areacello_ds.area_t), missing => NaN) # This is required for later multiplication
dht = readcubedata(dht_ds.dht)
lon_OM2 = readcubedata(areacello_ds.geolon_t)
lat_OM2 = readcubedata(areacello_ds.geolat_t)
lev = dht_ds.st_ocean

# Unfortunately ACCESS-OM2 raw data does not have coordinates of cell vertices
# So instead I go back to the source: the supergrids
gadi_supergrid_dir = "/g/data/xp65/public/apps/access_moppy_data/grids"
gridfile = "mom1deg.nc"
# CHECK that supergrid matches
supergrid_ds = open_dataset(joinpath(gadi_supergrid_dir, gridfile))
superarea = readcubedata(supergrid_ds.area)
lon = readcubedata(supergrid_ds.x)[2:2:end, 2:2:end]
ilon = .!ismissing.(lon_OM2.data) # Can't check full globe as area_t has missing values over some land
@assert lon_OM2[ilon] ≈ lon[ilon]
lat = readcubedata(supergrid_ds.y)[2:2:end, 2:2:end]
ilat = .!ismissing.(lat_OM2.data)
@assert lat_OM2[ilat] ≈ lat[ilat]
# I am using the supergrid area because I assume it is the ground truth here.
# I am not sure why the areas in the OM2 outputs and the supergrid files are different, by up to 6%!
# See outputs/plots/area_error.png
# TODO send email to ask around?
areacello = YAXArray(
    dims(dht)[1:2],
    [sum(superarea[i:(i + 1), j:(j + 1)]) for i in 1:2:size(superarea, 1) , j in 1:2:size(superarea, 2)],
    Dict("name" => "areacello", "units" => "m^2"),
)
iarea = .!isnan.(areacello_OM2.data) # isnan because I replaced missings with NaNs earlier
@show areacello_OM2[iarea] ≈ areacello[iarea]

# Build vertices from supergrid
# Dimensions of vertices ar (vertex, x, y)
# Note to self: NCO shows it as (y, x, vertex)
SW(x) = x[1:2:(end - 2), 1:2:(end - 2)]
SE(x) = x[3:2:end, 1:2:(end - 2)]
NE(x) = x[3:2:end, 3:2:end]
NW(x) = x[1:2:(end - 2), 3:2:end]
(nx, ny) = size(lon)
vertices(x) = [
    reshape(SW(x), (1, nx, ny))
    reshape(SE(x), (1, nx, ny))
    reshape(NE(x), (1, nx, ny))
    reshape(NW(x), (1, nx, ny))
]
lon_vertices = vertices(supergrid_ds.x)
lat_vertices = vertices(supergrid_ds.y)

@show size(lon_vertices)

volcello = readcubedata(dht .* areacello)

# Make makegridmetrics
gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
(; lon_vertices, lat_vertices, v3D) = gridmetrics

# Make indices
indices = makeindices(v3D)
(; wet3D) = indices

# Make V diagnoal matrix of volumes
V = sparse(Diagonal(v3D[wet3D]))

issrf = let
    issrf3D = falses(size(wet3D))
    issrf3D[:, :, 1] .= true
    issrf3D[wet3D]
end

TMfile = joinpath(inputdir, "yearly_matrix_$(κVdeep_str)_$(κH_str)_$(κVML_str).jld2")
@info "Loading matrix from $TMfile"
T = load(TMfile, "T")
@info "Matrix size: $(size(T)), nnz = $(nnz(T))"

# Following Holzer et al. (2020) or Pasquier et al. (2024),
#     Γꜛᵢ = Tᵃᵢᵢ⁻¹ 1ᵢ
# is the adjoint water age, the mean time until next surface contact.
# The adjoint transport matrix is just
#     Tᵃ = V⁻¹ * Tᵀ * V

idx_interior = findall(.!issrf)
v = v3D[wet3D]
V = sparse(Diagonal(v))
V⁻¹ = sparse(Diagonal(1 ./ v))
Tᵃ = V⁻¹ * transpose(T) * V
Tᵃᵢᵢ = Tᵃ[idx_interior, idx_interior]
Nᵢ = length(idx_interior)
onesᵢ = ones(Nᵢ)

@info "Solve full problem but with MLKPardisoIterate"
matrix_type = Pardiso.REAL_SYM
@show solver = MKLPardisoIterate(; nprocs, matrix_type)
prob = init(LinearProblem(Tᵃᵢᵢ, onesᵢ), solver, rtol = 1.0e-10)
@time "Pardiso solve" Γꜛᵢ = solve!(prob).u

# Check error magnitude
@show sol_error = norm(Tᵃᵢᵢ * Γꜛᵢ - onesᵢ) / norm(onesᵢ)

# turn the age solution vector back into a 3D yax
@show Γꜛyax = YAXArray(
    dims(volcello),
    ustrip.(yr, OceanTransportMatrixBuilder.as3D([zeros(sum(issrf)); Γꜛᵢ], wet3D) * s),
    Dict(
        "description" => "steady-state adjoint mean age (time until next surface contact)",
        "solver" => "MKLPardisoIterate",
        "model" => model,
        "experiment" => experiment,
        "time window" => time_window,
        "upwind" => upwind_str2,
        "units" => "yr",
    )
)
arrays = Dict(:Gammaup => Γꜛyax, :lat => lat, :lon => lon)
ds = Dataset(; properties = Dict(), arrays...)
# Save to netCDF file
outputfile = joinpath(inputdir, "steady_age_$(κVdeep_str)_$(κH_str)_$(κVML_str)_MKLPardisoIterate.nc")
@info "Saving age as netCDF file:\n  $(outputfile)"
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)

# Following Holzer et al. (2020) or Pasquier et al. (2024) the volume 𝒱↓ is given by
#     𝒱ꜜ = −Aₛ⁻¹ Vₛ Tᵃₛᵢ Tᵃᵢᵢ⁻¹ 1ᵢ
# But this is the same as
#     𝒱ꜜ = −Aₛ⁻¹ Vₛ Tᵃₛᵢ Γꜛᵢ
# So I might as well compute 𝒱ꜜ now since I just computed Γꜛᵢ
# Unit is m⁻² m³ s⁻¹ s = interior volume (m³) / surface area (m²)
# Note: In Pasquier et al. (2024) I plot this as %(interior volume) / 10,000km²

idx_surface = findall(issrf)
Tᵃₛᵢ = Tᵃ[idx_surface, idx_interior]
Vₛ = sparse(Diagonal(v3D[wet3D][idx_surface]))
wet2D = wet3D[:, :, 1]
isurface2D = findall(wet2D)
Aₛ⁻¹ = sparse(Diagonal(1 ./ areacello.data[isurface2D]))

𝒱ꜜ = -Aₛ⁻¹ * Vₛ * Tᵃₛᵢ * Γꜛᵢ

# Save 𝒱↑ as netCDF file
𝒱ꜜ2D = fill(NaN, size(wet2D))
𝒱ꜜ2D[isurface2D] .= 𝒱ꜜ
@show 𝒱ꜜyax = YAXArray(
    dims(areacello),
    𝒱ꜜ2D,
    Dict(
        "description" => "steady-state ocean volume ventilated down by unit area",
        "units" => "m^3/m^2",
        "solver" => "MKLPardisoIterate",
        "upwind" => upwind_str2,
    )
)
arrays = Dict(:Vdown => 𝒱ꜜyax, :lat => lat, :lon => lon)
ds = Dataset(; properties = Dict(), arrays...)
# Save to netCDF file
outputfile = joinpath(inputdir, "steady_Vdown_$(κVdeep_str)_$(κH_str)_$(κVML_str)_MKLPardisoIterate.nc")
@info "Saving Vdown as netCDF file:\n  $(outputfile)"
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)


