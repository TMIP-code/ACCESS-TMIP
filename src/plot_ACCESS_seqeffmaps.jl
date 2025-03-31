# qsub -I -P xv83 -l mem=180GB -l storage=scratch/gh0+scratch/xv83 -l walltime=01:00:00 -l ncpus=48

using Pkg
Pkg.activate(".")
Pkg.instantiate()

using OceanTransportMatrixBuilder
using NetCDF
using YAXArrays
using DataFrames
using DimensionalData
# using SparseArrays
# using LinearAlgebra
using Unitful
using Unitful: s, yr
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
using Format
# using LaTeXStrings

model = "ACCESS-ESM1-5"

time_window = "Jan2030-Dec2039"
# time_window = "Jan2090-Dec2099"
experiment = parse(Int, time_window[4:7]) ≤ 2010 ? "historical" : "ssp370"

members = ["r$(r)i1p1f1" for r in 1:40]
# members = ["r$(r)i1p1f1" for r in 1:3]

# Gadi directory for input files
# inputdirfun(member) = "/scratch/xv83/TMIP/data/$model/$experiment/all members/$(time_window)"
inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/all_members/$(time_window)/cyclomonth"
outputdir = inputdir
mkpath(inputdir)



mlotst_files = ["/scratch/xv83/TMIP/data/$model/$experiment/$member/$time_window/cyclomonth/mlotst.nc" for member in members]
mlotst_ds = open_mfdataset(DimArray(mlotst_files, Dim{:member}(members)))
mlotst = readcubedata(mlotst_ds.mlotst)

mlotst_yearmax_ensemblemean = dropdims(maximum(maximum(mlotst, dims=:month), dims=:member), dims=(:month, :member)) .-
    dropdims(minimum(maximum(mlotst, dims=:month), dims=:member), dims=(:month, :member))

# Load areacello and volcello for grid geometry
fixedvarsinputdir = "/scratch/xv83/TMIP/data/$model"
volcello_ds = open_dataset(joinpath(fixedvarsinputdir, "volcello.nc"))
areacello_ds = open_dataset(joinpath(fixedvarsinputdir, "areacello.nc"))

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
(; lon_vertices, lat_vertices, lon, lat, zt, v3D, thkcello, Z3D) = gridmetrics
lev = zt
# Make indices
indices = makeindices(gridmetrics.v3D)
(; wet3D, N) = indices

# Load \Gamma out
κVdeep = 3.0e-5 # m^2/s
κVML = 1.0      # m^2/s
κH = 300.0      # m^2/s
κVdeep_str = "kVdeep" * format(κVdeep, conversion="e")
κVML_str = "kVML" * format(κVML, conversion="e")
κH_str = "kH" * format(κH, conversion="d")
upwind = false
upwind_str = upwind ? "" : "_centered"
upwind_str2 = upwind ? "upwind" : "centered"
yearly = true
yearly_str = yearly ? "_yearly" : ""
yearly_str2 = yearly ? "(yearly)" : ""

# TODO Adapt idea of just grabbing τs before from diff plot
if yearly
    ℰ_files = ["/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)/seqeff$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str).nc" for member in members]
    ℰ_ds = open_mfdataset(DimArray(ℰ_files, Dim{:member}(members)))
    ℰ = readcubedata(ℰ_ds.seqeff)
else
    ℰ_files = ["/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)/calE$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str).nc" for member in members]
    ℰ_ds = open_mfdataset(DimArray(ℰ_files, Dim{:member}(members)))
    ℰ = readcubedata(ℰ_ds.calE)
end



years = ℰ_ds.Ti |> Array

ℰ_ensemblemean = dropdims(mean(ℰ, dims = 4), dims = 4)
ℰ_ensemblemax = dropdims(maximum(ℰ, dims = 4), dims = 4)
ℰ_ensemblemin = dropdims(minimum(ℰ, dims = 4), dims = 4)
ℰ_ensemblerange = ℰ_ensemblemax - ℰ_ensemblemin


include("plotting_functions.jl")

usecontourf = false

axs = Array{Any,2}(undef, (3, 2))
contours = Array{Any,2}(undef, (3, 2))
nrows, ncols = size(axs)

fig = Figure(size = (ncols * 500, nrows * 250 + 100), fontsize = 18)

yticks = -60:30:60
xticks = -120:60:120 + 360

for (irow, year) in enumerate([100, 300, 1000])

    iyear = year + 1

    # Plot ensemble mean
    icol = 1
    levels = 0:10:100
    colormap = cgrad(:Zissou1Continuous, length(levels) - 1, categorical = true, rev = true)
    # colormap = cgrad(:Hiroshige, 10, categorical = true, rev = true)
    colorrange = extrema(levels)

    axs[irow, icol] = ax = Axis(fig[irow, icol]; yticks, xticks, xtickformat, ytickformat)

    contours[irow, icol] = if usecontourf
        plotcontourfmap!(ax, 100 * ℰ_ensemblemean[:, :, iyear], gridmetrics; levels, colormap)
    else
        plotmap!(ax, 100 * ℰ_ensemblemean[:, :, iyear], gridmetrics; colorrange, colormap) # <- need to fix wrapping longitude for contour levels
    end

    myhidexdecorations!(ax, irow < nrows)
    myhideydecorations!(ax, icol > 1)

    # Plot ensemble range
    icol = 2
    levels = 0:10:60
    colormap = cgrad(:tol_ylorbr, length(levels), categorical = true)
    highclip = colormap[end]
    colormap = cgrad(colormap[1:end-1], categorical = true)
    colorrange = extrema(levels)

    axs[irow, icol] = ax = Axis(fig[irow, icol]; yticks, xticks, xtickformat, ytickformat)

    contours[irow, icol] = if usecontourf
        plotcontourfmap!(ax, 100 * ℰ_ensemblerange[:, :, iyear], gridmetrics; levels, colormap)
    else
        contours[irow, icol] = plotmap!(ax, 100 * ℰ_ensemblerange[:, :, iyear], gridmetrics; colorrange, colormap, highclip) # <- need to fix wrapping longitude for contour levels
    end

    myhidexdecorations!(ax, irow < nrows)
    myhideydecorations!(ax, icol > 1)

    Label(fig[irow, 0]; text = "τ = $year years", rotation = π/2, tellheight = false)

end

# ℰstr = rich("ℰ", superscript("—"ℰ̅ ‾))
ℰstr = rich("ℰ", rich("‾", offset = (-0.5, 0.15)))
# τstr = rich("τ", font = :italic)
𝒓 = rich("r", font = :bold_italic)
ℰfun = rich(ℰstr, "(", 𝒓, ", τ)")
label = rich("ensemble mean ", ℰfun, " (%)")
cb = Colorbar(fig[nrows + 1, 1], contours[1, 1]; label, vertical = false, flipaxis = false, ticks = 0:20:100)
cb.width = Relative(2/3)

label = rich("ensemble range, max ", ℰfun, " − min ", ℰfun, " (%)")
cb = Colorbar(fig[nrows + 1, 2], contours[1, 2]; label, vertical = false, flipaxis = false, ticks = 0:10:60)
cb.width = Relative(2/3)

# column labels
# Label(fig[0, 1]; text = "ensemble mean", tellwidth = false)
# Label(fig[0, 2]; text = "ensemble range (internal variability)", tellwidth = false)

labels = [
    "a" "d"
    "b" "e"
    "c" "f"
]

for (ax, label) in zip(axs, labels)
    txt = text!(ax, 0, 1; text = "$label", labeloptions..., strokecolor = :white, strokewidth = 3)
    translate!(txt, 0, 0, 100)
    txt = text!(ax, 0, 1; text = "$label", labeloptions...)
    translate!(txt, 0, 0, 100)
end

Label(fig[0, 1:2]; text = "$(time_window[4:7])s Seafloor Sequestration Efficiency ($(length(members)) members)$(yearly_str2)", fontsize = 24, tellwidth = false)
rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)

# save plot
suffix = usecontourf ? "_ctrf" : ""

outputfile = joinpath(outputdir, "seqeff$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str)_$(time_window)$(suffix).png")
@info "Saving seqeff on sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)
outputfile = joinpath(outputdir, "seqeff$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str)_$(time_window)$(suffix).pdf")
@info "Saving seqeff on sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)

