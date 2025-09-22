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
inputdir = "/scratch/y99/TMIP/data/ACCESS-OM2-025/omip2/r1i1p1f1/Jan0200-Dec0209/"

# Load areacello and volcello for grid geometry
umo_ds = open_dataset(joinpath(inputdir, "umo.nc"))
vmo_ds = open_dataset(joinpath(inputdir, "vmo.nc"))
ψᵢGM_ds = open_dataset(joinpath(inputdir, "tx_trans_gm.nc"))
ψⱼGM_ds = open_dataset(joinpath(inputdir, "ty_trans_gm.nc"))
ψᵢsubmeso_ds = open_dataset(joinpath(inputdir, "tx_trans_submeso.nc"))
ψⱼsubmeso_ds = open_dataset(joinpath(inputdir, "ty_trans_submeso.nc"))
mlotst_ds = open_dataset(joinpath(inputdir, "mlotst.nc"))
areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))
volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))


# Load fixed variables in memory
areacello = readcubedata(areacello_ds.areacello)
volcello = readcubedata(volcello_ds.volcello)
lon = readcubedata(volcello_ds.lon)
lat = readcubedata(volcello_ds.lat)
lev = volcello_ds.lev
# Identify the vertices keys (vary across CMIPs / models)
volcello_keys = propertynames(volcello_ds)
lon_vertices_key = volcello_keys[findfirst(x -> occursin("lon", x) & occursin("vert", x), string.(volcello_keys))]
lat_vertices_key = volcello_keys[findfirst(x -> occursin("lat", x) & occursin("vert", x), string.(volcello_keys))]
lon_vertices = readcubedata(getproperty(volcello_ds, lon_vertices_key))
lat_vertices = readcubedata(getproperty(volcello_ds, lat_vertices_key))

# Make makegridmetrics
gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
(; lon_vertices, lat_vertices) = gridmetrics

# Make indices
indices = makeindices(gridmetrics.v3D)


# Some parameter values
ρ = 1035.0    # kg/m^3
# ACCESS-ESM1.5 preferred values
upwind = false
κVdeep = 3.0e-5 # m^2/s
κVML = 1.0      # m^2/s
κH = 300.0      # m^2/s


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
outputfile = joinpath(inputdir, "cyclo_matrix_$(κVdeep_str)_$(κH_str)_$(κVML_str)_meanflow.jld2")
@info "Saving matrix as $outputfile"
save(outputfile,
    Dict(
        "T" => T,
        "note" => "Test 0.25-degree transport matrix built from averaging the all the transport variables umo and vmo (no GM nor subeso terms) + mlotst."
    )
)





