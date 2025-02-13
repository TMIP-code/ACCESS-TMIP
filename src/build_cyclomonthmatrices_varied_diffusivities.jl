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
using FileIO

# Making monthly matrices for every MONTH from 1990s and averaging them into a single matrix, and computing the age

# script options
@show model = "ACCESS-ESM1-5"
if isempty(ARGS)
    member = "r20i1p1f1"
    experiment = "historical"
    time_window = "Jan1850-Dec1859"
    # time_window = "Jan1990-Dec1999"
    # experiment = "ssp370"
    # time_window = "Jan2030-Dec2039"
    # time_window = "Jan2090-Dec2099"
else
    experiment, member, time_window = ARGS
end
@show experiment
@show member
@show time_window

lumpby = "month"
months = 1:12

# Gadi directory for input files
fixedvarsinputdir = "/scratch/xv83/TMIP/data/$model"
# Load areacello and volcello for grid geometry
volcello_ds = open_dataset(joinpath(fixedvarsinputdir, "volcello.nc"))
areacello_ds = open_dataset(joinpath(fixedvarsinputdir, "areacello.nc"))
inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)"
cycloinputdir = joinpath(inputdir, "cyclo$lumpby")
umo_ds = open_dataset(joinpath(cycloinputdir, "umo.nc"))
vmo_ds = open_dataset(joinpath(cycloinputdir, "vmo.nc"))
ψᵢGM_ds = open_dataset(joinpath(cycloinputdir, "tx_trans_gm.nc"))
ψⱼGM_ds = open_dataset(joinpath(cycloinputdir, "ty_trans_gm.nc"))
ψᵢsubmeso_ds = open_dataset(joinpath(cycloinputdir, "tx_trans_submeso.nc"))
ψⱼsubmeso_ds = open_dataset(joinpath(cycloinputdir, "ty_trans_submeso.nc"))
mlotst_ds = open_dataset(joinpath(cycloinputdir, "mlotst.nc"))


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


for month in months

    mlotst = readcubedata(mlotst_ds.mlotst[month=At(month)])
    umo = readcubedata(umo_ds.umo[month=At(month)])
    vmo = readcubedata(vmo_ds.vmo[month=At(month)])

    mean_days_in_month = umo_ds.mean_days_in_month[month=At(month)] |> Array |> only

    ψᵢGM = readcubedata(ψᵢGM_ds.tx_trans_gm[month=At(month)])
    ψⱼGM = readcubedata(ψⱼGM_ds.ty_trans_gm[month=At(month)])
    ψᵢsubmeso = readcubedata(ψᵢsubmeso_ds.tx_trans_submeso[month=At(month)])
    ψⱼsubmeso = readcubedata(ψⱼsubmeso_ds.ty_trans_submeso[month=At(month)])

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
    outputfile = joinpath(cycloinputdir, "cyclo_matrix_$(κVdeep_str)_$(κH_str)_$(κVML_str)_$(month).jld2")
    @info "Saving matrix as $outputfile"
    save(outputfile,
        Dict(
            "T" => T,
            "mean_days_in_month" => mean_days_in_month,
        )
    )


end



