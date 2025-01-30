# qsub -I -P xv83 -l mem=64GB -l storage=scratch/gh0+scratch/xv83 -l walltime=02:00:00 -l ncpus=12

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

model = "ACCESS-ESM1-5"

time_window1 = "Jan2030-Dec2039"
time_window2 = "Jan2090-Dec2099"
experiment1 = parse(Int, time_window1[4:7]) ≤ 2010 ? "historical" : "ssp370"
experiment2 = parse(Int, time_window2[4:7]) ≤ 2010 ? "historical" : "ssp370"

# members = ["r$(r)i1p1f1" for r in 1:40]
members = ["r$(r)i1p1f1" for r in 1:3]

# Gadi directory for input files
# inputdirfun(member) = "/scratch/xv83/TMIP/data/$model/$experiment/all members/$(time_window)"
inputdir1 = "/scratch/xv83/TMIP/data/$model/$experiment1/all_members/$(time_window1)/cyclomonth"
inputdir2 = "/scratch/xv83/TMIP/data/$model/$experiment2/all_members/$(time_window2)/cyclomonth"
outputdir = inputdir2
mkpath(inputdir2)



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



ℰ_files1 = ["/scratch/xv83/TMIP/data/$model/$experiment1/$member/$(time_window1)/calE.nc" for member in members]
ℰ_files2 = ["/scratch/xv83/TMIP/data/$model/$experiment2/$member/$(time_window2)/calE.nc" for member in members]
ℰ_ds1 = open_mfdataset(DimArray(ℰ_files1, Dim{:member}(members)))
ℰ_ds2 = open_mfdataset(DimArray(ℰ_files2, Dim{:member}(members)))
ℰ1 = readcubedata(ℰ_ds1.calE)
ℰ2 = readcubedata(ℰ_ds2.calE)



years = ℰ_ds1.Ti |> Array

# little helper function to get the year given a quantile
function yearatquantile(ℰ, ℰlevel)
    isnan(ℰ[1]) && return NaN
    out = findfirst(ℰ .< ℰlevel)
    isnothing(out) ? maximum(years) : Float64(out)
end
ℰ10_1 = map(
    ts -> yearatquantile(ts, 0.9),
    view(ℰ1, i, j, :, m) for i in 1:size(ℰ1,1), j in 1:size(ℰ1,2), m in 1:size(ℰ1,4)
)
ℰ10_1_ensemblemean = dropdims(mean(ℰ10_1, dims = 3), dims = 3)
ℰ10_2 = map(
    ts -> yearatquantile(ts, 0.9),
    view(ℰ2, i, j, :, m) for i in 1:size(ℰ2,1), j in 1:size(ℰ2,2), m in 1:size(ℰ2,4)
)
ℰ10_2_ensemblemean = dropdims(mean(ℰ10_2, dims = 3), dims = 3)
ℰ10_diff = ℰ10_2_ensemblemean - ℰ10_1_ensemblemean

ℰ50_1 = map(
    ts -> yearatquantile(ts, 0.5),
    view(ℰ1, i, j, :, m) for i in 1:size(ℰ1,1), j in 1:size(ℰ1,2), m in 1:size(ℰ1,4)
)
ℰ50_1_ensemblemean = dropdims(mean(ℰ50_1, dims = 3), dims = 3)
ℰ50_2 = map(
    ts -> yearatquantile(ts, 0.5),
    view(ℰ2, i, j, :, m) for i in 1:size(ℰ2,1), j in 1:size(ℰ2,2), m in 1:size(ℰ2,4)
)
ℰ50_2_ensemblemean = dropdims(mean(ℰ50_2, dims = 3), dims = 3)
ℰ50_diff = ℰ50_2_ensemblemean - ℰ50_1_ensemblemean



# Load \Gamma out
Gammaoutfile1 = "/scratch/xv83/TMIP/data/$model/$experiment1/all_members/$(time_window1)/cyclomonth/adjointage_timemean.nc"
Γoutyr3D1_timemean = readcubedata(open_dataset(Gammaoutfile1).adjointage_timemean)
Γoutyr3D1_ensemblemean = dropdims(mean(Γoutyr3D1_timemean, dims = 4), dims = 4)
Gammaoutfile2 = "/scratch/xv83/TMIP/data/$model/$experiment2/all_members/$(time_window2)/cyclomonth/adjointage_timemean.nc"
Γoutyr3D2_timemean = readcubedata(open_dataset(Gammaoutfile2).adjointage_timemean)
Γoutyr3D2_ensemblemean = dropdims(mean(Γoutyr3D2_timemean, dims = 4), dims = 4)

include("plotting_functions.jl")

Γout1_ensemblemean = seafloorvalue(Γoutyr3D1_ensemblemean, wet3D, gridmetrics)
Γout2_ensemblemean = seafloorvalue(Γoutyr3D2_ensemblemean, wet3D, gridmetrics)
Γout_ensemblemean_diff = Γout2_ensemblemean - Γout1_ensemblemean

usecontourf = false

axs = Array{Any,2}(undef, (3, 2))
contours = Array{Any,2}(undef, (3, 2))
nrows, ncols = size(axs)

fig = Figure(size = (ncols * 500, nrows * 250 + 100), fontsize = 18)

yticks = -60:30:60
xticks = -120:60:120 + 360

datamean = (Γout2_ensemblemean, ℰ50_2_ensemblemean, ℰ10_2_ensemblemean)
datadiff = (Γout_ensemblemean_diff, ℰ50_diff, ℰ10_diff)
Γstr = rich("Γ", superscript("↑"))
Qstr = rich("Q", font = :italic)
Q10 = rich(Qstr, "(0.1)")
Q50 = rich(Qstr, "(0.5)")
rowlabels = (rich("Mean, ", Γstr), rich("Median, ", Q50), rich("10th percentile, ", Q10))

for (irow, (x2Dmean, x2Ddiff, text)) in enumerate(zip(datamean, datadiff, rowlabels))

    # Plot ensemble mean
    icol = 1
    levels = 0:200:2000
    colormap = cgrad(:viridis, 11, categorical = true)
    highclip = colormap[end]
    colormap = cgrad(colormap[1:end-1], categorical = true)
    colorrange = extrema(levels)

    axs[irow, icol] = ax = Axis(fig[irow, icol]; yticks, xticks, xtickformat, ytickformat)

    contours[irow, icol] = if usecontourf
        plotcontourfmap!(ax, x2Dmean, gridmetrics; levels, colormap)
    else
        plotmap!(ax, x2Dmean, gridmetrics; colorrange, colormap, highclip) # <- need to fix wrapping longitude for contour levels
    end

    myhidexdecorations!(ax, irow < nrows)
    myhideydecorations!(ax, icol > 1)

    # Plot ensemble diff
    icol = 2
    levels = -300:50:300
    colormap = cgrad(cgrad(:balance, 13, categorical = true)[[1:end÷2+1; end÷2+1:end]], categorical = true)
    highclip = colormap[end]
    lowclip = colormap[1]
    colormap = cgrad(colormap[2:end-1], categorical = true)
    colorrange = extrema(levels)

    axs[irow, icol] = ax = Axis(fig[irow, icol]; yticks, xticks, xtickformat, ytickformat)

    contours[irow, icol] = if usecontourf
        plotcontourfmap!(ax, x2Ddiff, gridmetrics; levels, colormap)
    else
        plotmap!(ax, x2Ddiff, gridmetrics; colorrange, colormap, lowclip, highclip) # <- need to fix wrapping longitude for contour levels
    end

    myhidexdecorations!(ax, irow < nrows)
    myhideydecorations!(ax, icol > 1)

    Label(fig[irow, 0]; text, rotation = π/2, tellheight = false)

end


label = rich("ensemble mean (years)")
cb = Colorbar(fig[nrows + 1, 1], contours[1, 1]; label, vertical = false, flipaxis = false, ticks = 0:400:2000)
cb.width = Relative(2/3)

label = rich("mean 2090s − mean 2030s (years)")
cb = Colorbar(fig[nrows + 1, 2], contours[1, 2]; label, vertical = false, flipaxis = false, ticks = -300:200:300, tickformat = divergingcbarticklabelformat)
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

Label(fig[0, 1:2]; text = "Climate Change Effect on Seafloor Reemergence Time ($(length(members)) members)", fontsize = 24, tellwidth = false)
rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)

# save plot
suffix = usecontourf ? "_ctrf" : ""
outputfile = joinpath(outputdir, "reemergencetime_diff_$(time_window2)$(suffix).png")
@info "Saving reemergencetime on sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)
outputfile = joinpath(outputdir, "reemergencetime_diff_$(time_window2)$(suffix).pdf")
@info "Saving reemergencetime on sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)

