# qsub -P y99 -N yearlyTM_OM2-1 -l ncpus=24 -l mem=95GB -l jobfs=4GB -l walltime=1:00:00 -l storage=scratch/gh0+gdata/xv83+scratch/xv83+scratch/y99+gdata/cj50+gdata/ik11 -l wd -o output/PBS/ -j oe

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
areacello_ds = open_dataset(joinpath(inputdir, "area_t_grid.nc"))
dzt_ds = open_dataset(joinpath(inputdir, "dht.nc")) # <- (new) cell thickness?
# TODO: caputre correlations between transport and dzt
# z* coordinate varies with time in ACCESS-OM2
# volcello_ds = open_dataset(joinpath(fixedvarsinputdir, "volcello.nc")) # <- not in ACCESS-OM2; must be built from dzt * area
umo_ds = open_dataset(joinpath(inputdir, "tx_trans.nc"))
vmo_ds = open_dataset(joinpath(inputdir, "ty_trans.nc"))
ψᵢGM_ds = open_dataset(joinpath(inputdir, "tx_trans_gm.nc"))
ψⱼGM_ds = open_dataset(joinpath(inputdir, "ty_trans_gm.nc"))
# ψᵢsubmeso_ds = open_dataset(joinpath(inputdir, "tx_trans_submeso.nc"))
# ψⱼsubmeso_ds = open_dataset(joinpath(inputdir, "ty_trans_submeso.nc"))
mlotst_ds = open_dataset(joinpath(inputdir, "mld.nc"))

# Load fixed variables in memory
areacello = readcubedata(areacello_ds.area_t)
dzt = readcubedata(dzt_ds.dzt)
volcello = dzt .* areacello
lon = readcubedata(areacello_ds.lon)
lat = readcubedata(areacello_ds.lat)
lev = areacello_ds.lev
# Identify the vertices keys (vary across CMIPs / models)
volcello_keys = propertynames(areacello_ds)
lon_vertices_key = volcello_keys[findfirst(x -> occursin("lon", x) & occursin("vert", x), string.(volcello_keys))]
lat_vertices_key = volcello_keys[findfirst(x -> occursin("lat", x) & occursin("vert", x), string.(volcello_keys))]
lon_vertices = readcubedata(getproperty(areacello_ds, lon_vertices_key))
lat_vertices = readcubedata(getproperty(areacello_ds, lat_vertices_key))

# Make makegridmetrics
gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
(; lon_vertices, lat_vertices) = gridmetrics

# Make indices
indices = makeindices(gridmetrics.v3D)

# Some parameter values
ρ = 1035.0    # kg/m^3
# κVML = 0.1    # m^2/s
# κVMLs = [0.001, 0.01, 0.1, 1] # m^2/s
# κVMLs = [1e-7, 1e-6, 1e-5, 1e-4] # m^2/s
# κVdeep = 1e-5 # m^2/s
# κH = 500.0    # m^2/s
κVML = 1e-7    # m^2/s
κVdeep = 1e-7 # m^2/s
κH = 5    # m^2/s
# κVdeeps = [1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4] # m^2/s
# κHs = [50, 150, 500, 1500, 5000] # m^2/s
# κVdeeps = [1e-7, 3e-7] # m^2/s
# κHs = [50, 150, 500, 1500, 5000] # m^2/s
# κVdeeps = [1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4] # m^2/s
# κHs = [5, 15] # m^2/s

mean_days_in_month = umo_ds.mean_days_in_month |> Array
w = Weights(mean_days_in_month)

mlotst = dropdims(mean(readcubedata(mlotst_ds.mlotst), w; dims=:month), dims=:month)
umo = dropdims(mean(readcubedata(umo_ds.umo), w; dims=:month), dims=:month)
vmo = dropdims(mean(readcubedata(vmo_ds.vmo), w; dims=:month), dims=:month)

ψᵢGM = dropdims(mean(readcubedata(ψᵢGM_ds.tx_trans_gm), w; dims=:month), dims=:month)
ψⱼGM = dropdims(mean(readcubedata(ψⱼGM_ds.ty_trans_gm), w; dims=:month), dims=:month)
ψᵢsubmeso = dropdims(mean(readcubedata(ψᵢsubmeso_ds.tx_trans_submeso), w; dims=:month), dims=:month)
ψⱼsubmeso = dropdims(mean(readcubedata(ψⱼsubmeso_ds.ty_trans_submeso), w; dims=:month), dims=:month)

# Replace missing values and convert to arrays
# I think latest YAXArrays converts _FillValues to missing
ψᵢGM = replace(ψᵢGM |> Array, missing => 0) .|> Float64
ψⱼGM = replace(ψⱼGM |> Array, missing => 0) .|> Float64
ψᵢsubmeso = replace(ψᵢsubmeso |> Array, missing => 0) .|> Float64
ψⱼsubmeso = replace(ψⱼsubmeso |> Array, missing => 0) .|> Float64

# Also remove missings in umo and vmo
umo = replace(umo, missing => 0)
vmo = replace(vmo, missing => 0)

# Take the vertical diff of zonal/meridional transport diagnostics to get their mass transport
(nx, ny, _) = size(ψᵢGM)
ϕᵢGM = diff([fill(0.0, nx, ny, 1);;; ψᵢGM |> Array], dims=3)
ϕⱼGM = diff([fill(0.0, nx, ny, 1);;; ψⱼGM |> Array], dims=3)
ϕᵢsubmeso = diff([fill(0.0, nx, ny, 1);;; ψᵢsubmeso |> Array], dims=3)
ϕⱼsubmeso = diff([fill(0.0, nx, ny, 1);;; ψⱼsubmeso |> Array], dims=3)

# TODO fix incompatible dimensions betwewen umo and ϕᵢGM/ϕᵢsubmeso Dim{:i} and Dim{:xu_ocean}
ϕ = let umo = umo + ϕᵢGM + ϕᵢsubmeso, vmo = vmo + ϕⱼGM + ϕⱼsubmeso
    facefluxesfrommasstransport(; umo, vmo, gridmetrics, indices)
end


(; T) = transportmatrix(; ϕ, mlotst, gridmetrics, indices, ρ, κH, κVML, κVdeep)

# Save cyclo matrix only (don't save all the metadata in case IO is a bottleneck)
κVdeep_str = "kVdeep" * format(κVdeep, conversion="e")
κH_str = "kH" * format(κH, conversion="d")
κVML_str = "kVML" * format(κVML, conversion="e")
outputfile = joinpath(cycloinputdir, "cyclo_matrix_$(κVdeep_str)_$(κH_str)_$(κVML_str)_meanflow.jld2")
@info "Saving matrix as $outputfile"
save(outputfile,
    Dict(
        "T" => T,
        "note" => "Test matrix built from averaging the all the transport variables umo and vmo + GM + subeso + mlotst."
    )
)





