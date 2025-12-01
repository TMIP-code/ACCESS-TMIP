# qsub -I -P y99 -N yearlyTM_OM2-025 -l ncpus=24 -l mem=95GB -l jobfs=4GB -l walltime=1:00:00 -l storage=scratch/gh0+gdata/xp65+scratch/xv83+scratch/y99 -l wd -o output/PBS/ -j oe

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

# Making yearly matrices for ACCESS-OM2-01

model = "ACCESS-OM2-01"
experiment = "01deg_jra55v140_iaf_cycle4"
time_window = "Jan1960-Dec1979"
@show inputdir = "/scratch/y99/TMIP/data/$model/$experiment/$time_window"

@info "Load dzt for grid geometry"
dzt_ds = open_dataset(joinpath(inputdir, "dzt.nc")) # dht = dzt
dzt = readcubedata(dzt_ds.dzt)
lev = dzt_ds.st_ocean

@info "Load supergrid areacello + coordinates + vertices"
# Unfortunately ACCESS-OM2 raw data does not have coordinates of cell vertices
# So instead I go back to the source: the supergrids
include("supergrid.jl")
(; lon, lat, areacello, lon_vertices, lat_vertices) = supergrid(model; dims = dims(dzt_ds.dzt)[1:2])

@info "Compute volcello"
volcello = readcubedata(dzt .* areacello)

@info "Make makegridmetrics"
gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
(; lon_vertices, lat_vertices) = gridmetrics

@info "Make indices"
indices = makeindices(gridmetrics.v3D)

@info "Load transport and mixed layer depth data"
umo_ds = open_dataset(joinpath(inputdir, "tx_trans.nc"))
vmo_ds = open_dataset(joinpath(inputdir, "ty_trans.nc"))

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

# Replace missing values and convert to arrays
# I think latest YAXArrays converts _FillValues to missing
umo = replace(umo, missing => 0)
vmo = replace(vmo, missing => 0)

@info "build face fluxes from mass transport"
# TODO fix incompatible dimensions betwewen umo and ϕᵢGM/ϕᵢsubmeso Dim{:i} and Dim{:xu_ocean}
# ϕ = let umo = umo + ϕᵢGM + ϕᵢsubmeso, vmo = vmo + ϕⱼGM + ϕⱼsubmeso
ϕ = facefluxesfrommasstransport(; umo, vmo, gridmetrics, indices)

# Some parameter values
ρ = 1035.0    # kg/m^3
# ACCESS-ESM1.5 preferred values
upwind = false
κVdeep = 3.0e-5 # m^2/s
κVML = 1.0      # m^2/s
κH = 300.0 / 10 # m^2/s (÷10 as simple grid scaling, from 1° to 0.1°)

@info "Build the matrix"
(; T) = transportmatrix(; ϕ, mlotst, gridmetrics, indices, ρ, κH, κVML, κVdeep, upwind)

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
        "note" => """Test 0.1-degree transport matrix built from averaging
            all the transport variables umo and vmo + mlotst (no GM!).
            """,
    )
)
