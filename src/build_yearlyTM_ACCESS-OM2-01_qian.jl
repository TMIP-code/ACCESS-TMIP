# qsub -I -P y99 -N yrlyTM_OM2-01-qian -q hugemem -l ncpus=24 -l mem=735GB -l jobfs=4GB -l walltime=3:00:00 -l storage=scratch/gh0+gdata/xv83+scratch/xv83+scratch/p66+scratch/y99+gdata/cj50+gdata/ik11 -l wd -o output/PBS/ -j oe

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
model = "ACCESS-OM2-01"

# script options
inputdir = "/scratch/y99/TMIP/data/$model/01deg_jra55v13_ryf9091_qian_wthmp/Jan2150-Dec2159/"

# Load transport
tx_trans_ds = open_dataset(joinpath(inputdir, "tx_trans.nc"))
ty_trans_ds = open_dataset(joinpath(inputdir, "ty_trans.nc"))

# Load MLD
mld_ds = open_dataset(joinpath(inputdir, "mld.nc"))

# Load areacello and volcello for grid geometry
# But not available in raw output, so must rebuild from area_t and dzt
area_t_ds = open_dataset(joinpath(inputdir, "area_t.nc"))
dzt_ds = open_dataset(joinpath(inputdir, "dzt.nc"))
areacello = readcubedata(area_t_ds.area_t)
dzt = readcubedata(dzt_ds.dzt)
volcello = areacello .* dzt


# Load fixed variables in memory
# Identify the vertices keys (vary across CMIPs / models)
# FIXME using test CMORized data for vertices. Hopefully they match!
CMORtestfile = "/scratch/p66/yz9299/OM2_CMORised/umo_Omon_$(model)_omip1_r1i1p1f1_gn_214501-214503.nc"
CMORtest_ds = open_dataset(CMORtestfile)
lon = readcubedata(CMORtest_ds.longitude)
lat = readcubedata(CMORtest_ds.latitude)
lev = CMORtest_ds.lev
volcello_keys = propertynames(CMORtest_ds)
lon_vertices_key = volcello_keys[findfirst(x -> occursin("lon", x) & occursin("vert", x), string.(volcello_keys))]
lat_vertices_key = volcello_keys[findfirst(x -> occursin("lat", x) & occursin("vert", x), string.(volcello_keys))]
lon_vertices = readcubedata(getproperty(CMORtest_ds, lon_vertices_key))
lat_vertices = readcubedata(getproperty(CMORtest_ds, lat_vertices_key))

# FIXME checking vertices
# FIXME

using CairoMakie
using Makie.GeometryBasics


# volume (3D)
FillValue = dzt.properties["_FillValue"]
v3D = volcello |> Array{Union{Missing, Float64}}
v3D = replace(v3D, missing => NaN, 0 => NaN, FillValue => NaN)

# area (2D)
FillValue = areacello.properties["_FillValue"]
area2D = areacello |> Array{Union{Missing, Float64}}
area2D = replace(area2D, missing => NaN, 0 => NaN, FillValue => NaN)

# depth and cell height (3D)
thkcello = v3D ./ area2D
ZBOT3D = cumsum(thkcello, dims = 3)
Z3D = ZBOT3D - 0.5 * thkcello
zt = lev |> Array

lat = lat |> Array
lon = lon |> Array

# same with lon_vertices
lon_vertices = lon_vertices |> Array{Float64}
lat_vertices = lat_vertices |> Array{Float64}

# i = j = 1
i, j = 1025, 2500

fig = Figure()
ax = Axis(fig[1, 1]; title = "$(model) neighboring grid cell vertices (v2)", xlabel = "Longitude", ylabel = "Latitude", aspect = DataAspect())
plotgridcell!(ax, lon_vertices[:, i, j], lat_vertices[:, i, j]; color = (:blue, 0.3), strokewidth = 2)
text!(ax, mean(lon_vertices[:, i, j]), mean(lat_vertices[:, i, j]); text = "(i=$i,j=$j)", align = (:center, :center))
plotgridcell!(ax, lon_vertices[:, i+1, j], lat_vertices[:, i+1, j]; color = (:red, 0.3), strokewidth = 2)
text!(ax, mean(lon_vertices[:, i+1, j]), mean(lat_vertices[:, i+1, j]); text = "(i+1,j)", align = (:center, :center))
plotgridcell!(ax, lon_vertices[:, i, j+1], lat_vertices[:, i, j+1]; color = (:purple, 0.3), strokewidth = 2)
text!(ax, mean(lon_vertices[:, i, j+1]), mean(lat_vertices[:, i, j+1]); text = "(i,j+1)", align = (:center, :center))
plotgridcell!(ax, lon_vertices[:, i+1, j+1], lat_vertices[:, i+1, j+1]; color = (:orange, 0.3), strokewidth = 2)
text!(ax, mean(lon_vertices[:, i+1, j+1]), mean(lat_vertices[:, i+1, j+1]); text = "(i+1,j+1)", align = (:center, :center))

save(joinpath(inputdir, "$(model)_gridcell_vertices_check_$(i)_$(j)_v2.png"), fig)
println(joinpath(inputdir, "$(model)_gridcell_vertices_check_$(i)_$(j)_v2.png"))
# FIXME end
# FIXME end

# Add _FillValue to volcello (as OceanTransportMatrixBuilder expects it)
volcello.properties["_FillValue"] = FillValue

# DEBUGGING
OceanTransportMatrixBuilder.getgridtopology(lon_vertices, lat_vertices, lev) # unknowngrid topology
    nx = size(lon_vertices, 2)
    ny = size(lon_vertices, 3)
    nz = length(lev)
    # North pole vertices
    NPlon = @view lon_vertices[3:4,:,end]
    NPlat = @view lat_vertices[3:4,:,end]


# Make makegridmetrics
gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
(; lon_vertices, lat_vertices) = gridmetrics

# Make indices
indices = makeindices(gridmetrics.v3D)


# Some parameter values
ρ = 1035.0    # kg/m^3
# ACCESS-ESM1.5 preferred values
# FIXME should I make background diffusivities much smaller?
upwind = false
κVdeep = 3.0e-5 # m^2/s
κVML = 1.0      # m^2/s
κH = 300.0 / 10 # m^2/s (grid-scaling by sqrt(area))


umo = readcubedata(tx_trans_ds.tx_trans)
vmo = readcubedata(ty_trans_ds.ty_trans)
mlotst = readcubedata(mld_ds.mld)

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
        "note" => "Test 0.1-degree transport matrix built from averaging the transport variables tx_trans, ty_trans, and mld."
    )
)





