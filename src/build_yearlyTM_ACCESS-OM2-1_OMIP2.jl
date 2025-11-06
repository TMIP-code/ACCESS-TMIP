# qsub -I -P y99 -N yearlyTM_OM2-1 -l ncpus=24 -l mem=95GB -l jobfs=4GB -l walltime=1:00:00 -l storage=scratch/gh0+gdata/xp65+scratch/xv83+scratch/y99 -l wd -o output/PBS/ -j oe

using Pkg
Pkg.activate(".")
Pkg.instantiate()

ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
using OceanTransportMatrixBuilder
using NetCDF
using YAXArrays
using DataFrames
using DimensionalData
using SparseArrays
using LinearAlgebra
using Unitful
using Unitful: s, yr
using Format
using Dates
using Statistics
using StatsBase
using FileIO

# Making yearly matrices for ACCESS-OM2-1

# script options
inputdir = "/scratch/y99/TMIP/data/ACCESS-OM2-1/1deg_jra55_iaf_omip2_cycle6/Jan1960-Dec1979"

# Load areacello and volcello for grid geometry
areacello_ds = open_dataset(joinpath(inputdir, "area_t.nc"))
dht_ds = open_dataset(joinpath(inputdir, "dht.nc")) # <- (new) cell thickness?
# TODO: caputre correlations between transport and dht
# z* coordinate varies with time in ACCESS-OM2
# volcello_ds = open_dataset(joinpath(fixedvarsinputdir, "volcello.nc")) # <- not in ACCESS-OM2; must be built from dht * area


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
(; lon_vertices, lat_vertices) = gridmetrics

# Make indices
indices = makeindices(gridmetrics.v3D)

# Load transport and mixed layer depth data
umo_ds = open_dataset(joinpath(inputdir, "tx_trans.nc"))
vmo_ds = open_dataset(joinpath(inputdir, "ty_trans.nc"))
ψᵢGM_ds = open_dataset(joinpath(inputdir, "tx_trans_gm.nc"))
ψⱼGM_ds = open_dataset(joinpath(inputdir, "ty_trans_gm.nc"))
# ψᵢsubmeso_ds = open_dataset(joinpath(inputdir, "tx_trans_submeso.nc"))
# ψⱼsubmeso_ds = open_dataset(joinpath(inputdir, "ty_trans_submeso.nc"))
mlotst_ds = open_dataset(joinpath(inputdir, "mld.nc"))

# FIXME no need to average over time for now as the input is already averaged
# FIXME first build monthly climatologies and then average matrices
# mean_days_in_month = umo_ds.mean_days_in_month |> Array
# w = Weights(mean_days_in_month)

# mlotst = dropdims(mean(readcubedata(mlotst_ds.mlotst), w; dims=:month), dims=:month)
# umo = dropdims(mean(readcubedata(umo_ds.umo), w; dims=:month), dims=:month)
# vmo = dropdims(mean(readcubedata(vmo_ds.vmo), w; dims=:month), dims=:month)
umo = readcubedata(umo_ds.tx_trans)
vmo = readcubedata(vmo_ds.ty_trans)
mlotst = readcubedata(mlotst_ds.mld)

# ψᵢGM = dropdims(mean(readcubedata(ψᵢGM_ds.tx_trans_gm), w; dims=:month), dims=:month)
# ψⱼGM = dropdims(mean(readcubedata(ψⱼGM_ds.ty_trans_gm), w; dims=:month), dims=:month)
ψᵢGM = readcubedata(ψᵢGM_ds.tx_trans_gm)
ψⱼGM = readcubedata(ψⱼGM_ds.ty_trans_gm)
# ψᵢsubmeso = dropdims(mean(readcubedata(ψᵢsubmeso_ds.tx_trans_submeso), w; dims=:month), dims=:month)
# ψⱼsubmeso = dropdims(mean(readcubedata(ψⱼsubmeso_ds.ty_trans_submeso), w; dims=:month), dims=:month)

# Replace missing values and convert to arrays
# I think latest YAXArrays converts _FillValues to missing
umo = replace(umo, missing => 0)
vmo = replace(vmo, missing => 0)
ψᵢGM = replace(ψᵢGM |> Array, missing => 0) .|> Float64
ψⱼGM = replace(ψⱼGM |> Array, missing => 0) .|> Float64
# ψᵢsubmeso = replace(ψᵢsubmeso |> Array, missing => 0) .|> Float64
# ψⱼsubmeso = replace(ψⱼsubmeso |> Array, missing => 0) .|> Float64

# Take the vertical diff of zonal/meridional transport diagnostics to get their mass transport
(nx, ny, _) = size(ψᵢGM)
ϕᵢGM = YAXArray(dims(umo), diff([fill(0.0, nx, ny, 1);;; ψᵢGM], dims = 3), umo.properties)
ϕⱼGM = YAXArray(dims(vmo), diff([fill(0.0, nx, ny, 1);;; ψⱼGM], dims = 3), vmo.properties)
# ϕᵢsubmeso = diff([fill(0.0, nx, ny, 1);;; ψᵢsubmeso |> Array], dims=3)
# ϕⱼsubmeso = diff([fill(0.0, nx, ny, 1);;; ψⱼsubmeso |> Array], dims=3)

# TODO fix incompatible dimensions betwewen umo and ϕᵢGM/ϕᵢsubmeso Dim{:i} and Dim{:xu_ocean}
# ϕ = let umo = umo + ϕᵢGM + ϕᵢsubmeso, vmo = vmo + ϕⱼGM + ϕⱼsubmeso
"""
sum YAXArrays a and b, preserving metadata from a
"""
yaxasum(a, b; axlist = dims(a), metadata = a.properties) = YAXArray(axlist, a.data + b.data, metadata)
ϕ = let umo = yaxasum(umo, ϕᵢGM), vmo = yaxasum(vmo, ϕⱼGM)
    facefluxesfrommasstransport(; umo, vmo, gridmetrics, indices)
end

# Some parameter values
ρ = 1035.0    # kg/m^3
# ACCESS-ESM1.5 preferred values
upwind = false
κVdeep = 3.0e-5 # m^2/s
κVML = 1.0      # m^2/s
κH = 300.0      # m^2/s

# Build the matrix
(; T) = transportmatrix(; ϕ, mlotst, gridmetrics, indices, ρ, κH, κVML, κVdeep)

# Save cyclo matrix only (don't save all the metadata in case IO is a bottleneck)
κVdeep_str = "kVdeep" * format(κVdeep, conversion = "e")
κH_str = "kH" * format(κH, conversion = "d")
κVML_str = "kVML" * format(κVML, conversion = "e")
# outputfile = joinpath(inputdir, "cyclo_matrix_$(κVdeep_str)_$(κH_str)_$(κVML_str)_meanflow.jld2")
outputfile = joinpath(inputdir, "yearly_matrix_$(κVdeep_str)_$(κH_str)_$(κVML_str).jld2")
@info "Saving matrix as $outputfile"
save(
    outputfile,
    Dict(
        "T" => T,
        "note" => """Test 1-degree transport matrix built from averaging
            all the transport variables umo and vmo (+GM but no subeso terms) + mlotst.
            """,
    )
)
