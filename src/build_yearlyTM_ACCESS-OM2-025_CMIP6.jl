# qsub -I -P y99 -N yearlyTM_OM2-1 -l ncpus=28 -l mem=120GB -l jobfs=4GB -l walltime=3:00:00 -l storage=scratch/gh0+gdata/xv83+scratch/xv83+scratch/p66+scratch/y99+gdata/cj50+gdata/ik11 -l wd -o output/PBS/ -j oe
# For debugging no need to request that much resources!

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
# FIXME adding plotting functions temporarily for debugging
try
    using CairoMakie
catch
    using CairoMakie
end
using GeoMakie
using Interpolations
using OceanBasins
using Statistics
using NaNStatistics
using StatsBase
using FileIO
using Contour
using GeometryBasics
using GeometryOps
using LibGEOS
include("plotting_functions.jl")

# Making yearly matrices for ACCESS-OM2-025

# script options
inputdir = "/scratch/y99/TMIP/data/ACCESS-OM2-025/omip2/r1i1p1f1/Jan0200-Dec0209/"

# Load areacello and volcello for grid geometry
umo_ds = open_dataset(joinpath(inputdir, "umo.nc"))
vmo_ds = open_dataset(joinpath(inputdir, "vmo.nc"))
# FIXME No GM or submeso terms for now
# ψᵢGM_ds = open_dataset(joinpath(inputdir, "tx_trans_gm.nc"))
# ψⱼGM_ds = open_dataset(joinpath(inputdir, "ty_trans_gm.nc"))
# ψᵢsubmeso_ds = open_dataset(joinpath(inputdir, "tx_trans_submeso.nc"))
# ψⱼsubmeso_ds = open_dataset(joinpath(inputdir, "ty_trans_submeso.nc"))
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
# FIXME using test CMORized data for vertices. Hopefully they match!
CMORtestfile = "/scratch/p66/yz9299/OM2_CMORised/umo_Omon_ACCESS-OM2-025_omip1_r1i1p1f1_gn_190001-190112.nc"
CMORtest_ds = open_dataset(CMORtestfile)
volcello_keys = propertynames(CMORtest_ds)
lon_vertices_key = volcello_keys[findfirst(x -> occursin("lon", x) & occursin("vert", x), string.(volcello_keys))]
lat_vertices_key = volcello_keys[findfirst(x -> occursin("lat", x) & occursin("vert", x), string.(volcello_keys))]
lon_vertices = readcubedata(getproperty(CMORtest_ds, lon_vertices_key))
lat_vertices = readcubedata(getproperty(CMORtest_ds, lat_vertices_key))

# FIXME checking vertices
# FIXME

# using CairoMakie
# using Makie.GeometryBasics


# # volume (3D)
# FillValue = volcello.properties["_FillValue"]
# v3D = volcello |> Array{Union{Missing, Float64}}
# v3D = replace(v3D, missing => NaN, 0 => NaN, FillValue => NaN)

# # area (2D)
# FillValue = areacello.properties["_FillValue"]
# area2D = areacello |> Array{Union{Missing, Float64}}
# area2D = replace(area2D, missing => NaN, 0 => NaN, FillValue => NaN)

# # depth and cell height (3D)
# thkcello = v3D ./ area2D
# ZBOT3D = cumsum(thkcello, dims = 3)
# Z3D = ZBOT3D - 0.5 * thkcello
# zt = lev |> Array

# lat = lat |> Array
# lon = lon |> Array

# # same with lon_vertices
# lon_vertices = lon_vertices |> Array{Float64}
# lat_vertices = lat_vertices |> Array{Float64}

# i = j = 1
# i, j = 500, 1000

# fig = Figure()
# ax = Axis(fig[1, 1]; title = "ACCESS-OM2-025 neighboring grid cell vertices (v2)", xlabel = "Longitude", ylabel = "Latitude", aspect = DataAspect())
# plotgridcell!(ax, lon_vertices[:, i, j], lat_vertices[:, i, j]; color = (:blue, 0.3), strokewidth = 2)
# text!(ax, mean(lon_vertices[:, i, j]), mean(lat_vertices[:, i, j]); text = "(i=$i,j=$j)", align = (:center, :center))
# plotgridcell!(ax, lon_vertices[:, i+1, j], lat_vertices[:, i+1, j]; color = (:red, 0.3), strokewidth = 2)
# text!(ax, mean(lon_vertices[:, i+1, j]), mean(lat_vertices[:, i+1, j]); text = "(i+1,j)", align = (:center, :center))
# plotgridcell!(ax, lon_vertices[:, i, j+1], lat_vertices[:, i, j+1]; color = (:purple, 0.3), strokewidth = 2)
# text!(ax, mean(lon_vertices[:, i, j+1]), mean(lat_vertices[:, i, j+1]); text = "(i,j+1)", align = (:center, :center))
# plotgridcell!(ax, lon_vertices[:, i+1, j+1], lat_vertices[:, i+1, j+1]; color = (:orange, 0.3), strokewidth = 2)
# text!(ax, mean(lon_vertices[:, i+1, j+1]), mean(lat_vertices[:, i+1, j+1]); text = "(i+1,j+1)", align = (:center, :center))

# save(joinpath(inputdir, "gridcell_vertices_check_$(i)_$(j)_v2.png"), fig)
# println(joinpath(inputdir, "gridcell_vertices_check_$(i)_$(j)_v2.png"))
# FIXME end
# FIXME end


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
κH = 300.0 / 4  # m^2/s (grid-scaling by sqrt(area))


# FIXME no need to average over time for now as the input is already averaged
# FIXME but in future need to try to first build monthly climatologies and then average here
# mean_days_in_month = umo_ds.mean_days_in_month |> Array
# w = Weights(mean_days_in_month)

# mlotst = dropdims(mean(readcubedata(mlotst_ds.mlotst), w; dims=:month), dims=:month)
# umo = dropdims(mean(readcubedata(umo_ds.umo), w; dims=:month), dims=:month)
# vmo = dropdims(mean(readcubedata(vmo_ds.vmo), w; dims=:month), dims=:month)
umo = readcubedata(umo_ds.umo)
vmo = readcubedata(vmo_ds.vmo)
mlotst = readcubedata(mlotst_ds.mlotst)

# FIXME No GM or submeso terms for now
# ψᵢGM = dropdims(mean(readcubedata(ψᵢGM_ds.tx_trans_gm), w; dims=:month), dims=:month)
# ψⱼGM = dropdims(mean(readcubedata(ψⱼGM_ds.ty_trans_gm), w; dims=:month), dims=:month)
# ψᵢsubmeso = dropdims(mean(readcubedata(ψᵢsubmeso_ds.tx_trans_submeso), w; dims=:month), dims=:month)
# ψⱼsubmeso = dropdims(mean(readcubedata(ψⱼsubmeso_ds.ty_trans_submeso), w; dims=:month), dims=:month)
#
# # Replace missing values and convert to arrays
# # I think latest YAXArrays converts _FillValues to missing
# ψᵢGM = replace(ψᵢGM |> Array, missing => 0) .|> Float64
# ψⱼGM = replace(ψⱼGM |> Array, missing => 0) .|> Float64
# ψᵢsubmeso = replace(ψᵢsubmeso |> Array, missing => 0) .|> Float64
# ψⱼsubmeso = replace(ψⱼsubmeso |> Array, missing => 0) .|> Float64

# Also remove missings in umo and vmo
umo = replace(umo, missing => 0)
vmo = replace(vmo, missing => 0)

# # Take the vertical diff of zonal/meridional transport diagnostics to get their mass transport
# (nx, ny, _) = size(ψᵢGM)
# ϕᵢGM = diff([fill(0.0, nx, ny, 1);;; ψᵢGM |> Array], dims=3)
# ϕⱼGM = diff([fill(0.0, nx, ny, 1);;; ψⱼGM |> Array], dims=3)
# ϕᵢsubmeso = diff([fill(0.0, nx, ny, 1);;; ψᵢsubmeso |> Array], dims=3)
# ϕⱼsubmeso = diff([fill(0.0, nx, ny, 1);;; ψⱼsubmeso |> Array], dims=3)

# TODO fix incompatible dimensions betwewen umo and ϕᵢGM/ϕᵢsubmeso Dim{:i} and Dim{:xu_ocean}
# ϕ = let umo = umo + ϕᵢGM + ϕᵢsubmeso, vmo = vmo + ϕⱼGM + ϕⱼsubmeso
#     facefluxesfrommasstransport(; umo, vmo, gridmetrics, indices)
# end
ϕ = facefluxesfrommasstransport(; umo, vmo, gridmetrics, indices)


(; T) = transportmatrix(; ϕ, mlotst, gridmetrics, indices, ρ, κH, κVML, κVdeep)

# Save cyclo matrix only (don't save all the metadata in case IO is a bottleneck)
κVdeep_str = "kVdeep" * format(κVdeep, conversion="e")
κH_str = "kH" * format(κH, conversion="d")
κVML_str = "kVML" * format(κVML, conversion="e")
# outputfile = joinpath(inputdir, "cyclo_matrix_$(κVdeep_str)_$(κH_str)_$(κVML_str)_meanflow.jld2")
outputfile = joinpath(inputdir, "yearly_matrix_$(κVdeep_str)_$(κH_str)_$(κVML_str).jld2")
@info "Saving matrix as $outputfile"
save(outputfile,
    Dict(
        "T" => T,
        "note" => "Test 0.25-degree transport matrix built from averaging the all the transport variables umo and vmo (no GM nor subeso terms) + mlotst."
    )
)





