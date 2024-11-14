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
using NaNStatistics
using Format
using Dates
using FileIO

# Making monthly matrices for every MONTH from 1990s and averaging them into a single matrix, and computing the age

model = "ACCESS-ESM1-5"
member = "r1i1p1f1"
experiment = "historical"
time_window = "Jan1990-Dec1999"
lumpby = "season"


# Gadi directory for input files
inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)"
cycloinputdir = joinpath(inputdir, "cyclo$season")
umo_ds = open_dataset(joinpath(cycloinputdir, "umo.nc"))
vmo_ds = open_dataset(joinpath(cycloinputdir, "vmo.nc"))
ψᵢGM_ds = open_dataset(joinpath(cycloinputdir, "tx_trans_gm.nc"))
ψⱼGM_ds = open_dataset(joinpath(cycloinputdir, "ty_trans_gm.nc"))
ψᵢsubmeso_ds = open_dataset(joinpath(cycloinputdir, "tx_trans_submeso.nc"))
ψⱼsubmeso_ds = open_dataset(joinpath(cycloinputdir, "ty_trans_submeso.nc"))
mlotst_ds = open_dataset(joinpath(cycloinputdir, "mlotst.nc"))
areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))
volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))



# Load fixed variables in memory
areacello = readcubedata(areacello_ds.areacello)
volcello = readcubedata(volcello_ds.volcello)
lon = readcubedata(volcello_ds.longitude)
lat = readcubedata(volcello_ds.latitude)
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
κH = 500.0    # m^2/s
κVML = 0.1    # m^2/s
κVdeep = 1e-5 # m^2/s

seasons = ("DJF", "MAM", "JJA", "SON")

for season in seasons

    # Load variables in memory
    mlotst = readcubedata(mlotst_ds.mlotst[Ti=At(season)])
    umo = readcubedata(umo_ds.umo[Ti=At(season)])
    vmo = readcubedata(vmo_ds.vmo[Ti=At(season)])

    ψᵢGM = readcubedata(ψᵢGM_ds.tx_trans_gm[Ti=At(season)])
    ψⱼGM = readcubedata(ψⱼGM_ds.ty_trans_gm[Ti=At(season)])
    ψᵢsubmeso = readcubedata(ψᵢsubmeso_ds.tx_trans_submeso[Ti=At(season)])
    ψⱼsubmeso = readcubedata(ψⱼsubmeso_ds.ty_trans_submeso[Ti=At(season)])

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
    outputfile = joinpath(cycloinputdir, "cyclo_matrix_$season.jld2")
    @info "Saving matrix as $outputfile"
    save(outputfile,
        Dict(
            "T" => T,
        )
    )

end



