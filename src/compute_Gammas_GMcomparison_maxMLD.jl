using Pkg
Pkg.activate(".")
Pkg.instantiate()

using OceanTransportMatrixBuilder
using NetCDF
using YAXArrays
using DataFrames
using DimensionalData
using SparseArrays
using LinearAlgebra
using Unitful
using Unitful: s, yr
using NaNStatistics
using Format
using FileIO

# Load functions for GM terms
include("GentMcWilliams.jl")

model = "ACCESS-ESM1-5"
member = "r1i1p1f1"
CMIP_version = "CMIP5"
experiment = "historical"
time_window = "Jan1990-Dec1999"

# Gadi directory for input files
inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)"

# Load umo, vmo, mlotst, volcello, and areacello
umo_ds = open_dataset(joinpath(inputdir, "umo.nc"))
vmo_ds = open_dataset(joinpath(inputdir, "vmo.nc"))
mlotst_ds = open_dataset(joinpath(inputdir, "mlotst_max.nc"))
volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))

# Load variables in memory
umo = readcubedata(umo_ds.umo)
vmo = readcubedata(vmo_ds.vmo)
mlotst = readcubedata(mlotst_ds.mlotst)
areacello = readcubedata(areacello_ds.areacello)
volcello = readcubedata(volcello_ds.volcello)
lon = readcubedata(volcello_ds.lon)
lat = readcubedata(volcello_ds.lat)
lev = volcello_ds.lev

# Load GM terms
ϕᵢGM, ϕⱼGM = ϕfromACCESSESM15("gm", member, time_window)
ϕᵢsubmeso, ϕⱼsubmeso = ϕfromACCESSESM15("submeso", member, time_window)

# Identify the vertices keys (vary across CMIPs / models)
volcello_keys = propertynames(volcello_ds)
lon_vertices_key = volcello_keys[findfirst(x -> occursin("lon", x) & occursin("vert", x), string.(volcello_keys))]
lat_vertices_key = volcello_keys[findfirst(x -> occursin("lat", x) & occursin("vert", x), string.(volcello_keys))]
lon_vertices = readcubedata(getproperty(volcello_ds, lon_vertices_key))
lat_vertices = readcubedata(getproperty(volcello_ds, lat_vertices_key))

# Some parameter values
ρ = 1035.0    # kg/m^3
κH = 500.0    # m^2/s
κVML = 0.1    # m^2/s
κVdeep = 1.0e-5 # m^2/s

# Make makegridmetrics
gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
(; lon_vertices, lat_vertices) = gridmetrics

# Make indices
indices = makeindices(gridmetrics.v3D)

umos = (umo, umo + ϕᵢGM, umo + ϕᵢGM + ϕᵢsubmeso)
vmos = (vmo, vmo + ϕⱼGM, vmo + ϕⱼGM + ϕⱼsubmeso)
strs = ("resolved", "resolved_GM", "resolved_GM_submeso")

for (umo, vmo, str) in zip(umos, vmos, strs)

    # Make fuxes from all directions
    ϕ = facefluxesfrommasstransport(; umo, vmo, gridmetrics, indices)

    # Make transport matrix
    (; T, Tadv, TκH, TκVML, TκVdeep) = transportmatrix(; ϕ, mlotst, gridmetrics, indices, ρ, κH, κVML, κVdeep)

    # Save matrices
    outputfile = joinpath(inputdir, "transportmatrix_$(str)_maxMLD.jld2")
    @info "Saving matrices + metrics as $outputfile"
    save(
        outputfile,
        Dict(
            "T" => T,
            "Tadv" => Tadv,
            "TκH" => TκH,
            "TκVML" => TκVML,
            "TκVdeep" => TκVdeep,
            "gridmetrics" => gridmetrics,
            "indices" => indices,
            "κH" => κH,
            "κVML" => κVML,
            "κVdeep" => κVdeep,
            "model" => model,
            "experiment" => experiment,
            "member" => member,
            "time_window" => time_window,
            "note" => "built from averaging mass transport, $str",
            "OceanTransportMatrixBuilder" => "v$(pkgversion(OceanTransportMatrixBuilder))",
        )
    )

    # unpack model grid
    (; lon, lat, zt, v3D) = gridmetrics
    lev = zt
    # unpack indices
    (; wet3D, N) = indices

    v = v3D[wet3D]

    # surface mask
    issrf3D = copy(wet3D)
    issrf3D[:, :, 2:end] .= false
    issrf = issrf3D[wet3D]
    # Ideal mean age Γ↓ is governed by
    # 	∂Γ↓/∂t + T Γ↓ = 1 - M Γ↓
    # so
    # 	Γ↓ = A⁻¹ 1
    # where A = T + M and M is matrix mask of surface with short timescale (1s)
    sΓin = ones(size(v))
    M = sparse(Diagonal(issrf))
    @info "Factorizing A = T + M"
    A = factorize(T + M)
    @info "Solving ideal mean age"
    Γin = A \ sΓin
    Γinyr = ustrip.(yr, Γin .* s)
    Γinyr3D = OceanTransportMatrixBuilder.as3D(Γinyr, wet3D)

    # global mean of the mean age should be order 1000yr
    (v' * Γinyr) / sum(v)

    # Mean reemergence time is governed by
    # 	-∂Γ↑/∂t + (T̃ + M) Γ↑ = 1    (note the tilde on T̃ = V⁻¹ Tᵀ V, the adjoint of T)
    # For computational efficiency, reuse factors of A = T + M:
    #   Γ↑ = V⁻¹ A⁻ᵀ V 1
    V = sparse(Diagonal(v))
    V⁻¹ = sparse(Diagonal(1 ./ v))
    Γout = V⁻¹ * (transpose(A) \ (V * sΓin))
    Γoutyr = ustrip.(yr, Γout .* s)
    Γoutyr3D = OceanTransportMatrixBuilder.as3D(Γoutyr, wet3D)

    # global mean of the mean reemergence time should also be order 1000yr
    (v' * Γout) / sum(v)

    # Turn Γin into a YAXArray by rebuilding from volcello
    Γinyr_YAXArray = rebuild(
        volcello_ds["volcello"];
        data = Γinyr3D,
        dims = dims(volcello_ds["volcello"]),
        metadata = Dict(
            "origin" => "ideal_mean_age computed from $model $experiment $member $(time_window)",
            "units" => "yr",
        )
    )
    arrays = Dict(:age => Γinyr_YAXArray, :lat => volcello_ds.lat, :lon => volcello_ds.lon)
    Γin_ds = Dataset(; volcello_ds.properties, arrays...)

    # Turn Γout into a YAXArray by rebuilding from volcello
    Γoutyr_YAXArray = rebuild(
        volcello_ds["volcello"];
        data = Γoutyr3D,
        dims = dims(volcello_ds["volcello"]),
        metadata = Dict(
            "origin" => "mean_reemergence_time computed from $model $experiment $member $(time_window)",
            "units" => "yr",
        )
    )
    arrays = Dict(:age => Γoutyr_YAXArray, :lat => volcello_ds.lat, :lon => volcello_ds.lon)
    Γout_ds = Dataset(; volcello_ds.properties, arrays...)

    # Save Γinyr3D to netCDF file
    outputfile = joinpath(inputdir, "ideal_mean_age_$(str)_maxMLD.nc")
    @info "Saving ideal mean age as netCDF file:\n  $(outputfile)"
    savedataset(Γin_ds, path = outputfile, driver = :netcdf, overwrite = true)

    # Save Γoutyr3D to netCDF file
    outputfile = joinpath(inputdir, "mean_reemergence_time_$(str)_maxMLD.nc")
    @info "Saving mean reemergence time as netCDF file:\n  $(outputfile)"
    savedataset(Γout_ds, path = outputfile, driver = :netcdf, overwrite = true)


end
