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
using Format

model = "ACCESS-ESM1-5"

time_window = "Jan2030-Dec2039"
experiment = parse(Int, time_window[4:7]) ≤ 2010 ? "historical" : "ssp370"

members_mean = ["r$(r)i1p1f1" for r in 1:40]
members = ["r$(r)i1p1f1" for r in 1:3]

# Gadi directory for input files
# inputdirfun(member) = "/scratch/xv83/TMIP/data/$model/$experiment/all members/$(time_window)"
inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/all_members/$(time_window)/cyclomonth"
outputdir = inputdir
mkpath(inputdir)



# mlotst_files = ["/scratch/xv83/TMIP/data/$model/$experiment/$member/$time_window/cyclomonth/mlotst.nc" for member in members]
# mlotst_ds = open_mfdataset(DimArray(mlotst_files, Dim{:member}(members)))
# mlotst = readcubedata(mlotst_ds.mlotst)

# mlotst_yearmax_ensemblemean = dropdims(maximum(maximum(mlotst, dims=:month), dims=:member), dims=(:month, :member)) .-
#     dropdims(minimum(maximum(mlotst, dims=:month), dims=:member), dims=(:month, :member))

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



ℰ_files = ["/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)/calE.nc" for member in members]
ℰ_ds = open_mfdataset(DimArray(ℰ_files, Dim{:member}(members)))
ℰ = readcubedata(ℰ_ds.calE)



years = ℰ_ds.Ti |> Array

# little helper function to get the year given a quantile
function yearatquantile(ℰ, ℰlevel)
    isnan(ℰ[1]) && return NaN
    out = findfirst(ℰ .< ℰlevel)
    isnothing(out) ? maximum(years) : Float64(out)
end
ℰ10 = map(
    ts -> yearatquantile(ts, 0.9),
    view(ℰ, i, j, :, m) for i in 1:size(ℰ,1), j in 1:size(ℰ,2), m in 1:size(ℰ,4)
)
ℰ10_ensemblemean = dropdims(mean(ℰ10, dims = 3), dims = 3)
ℰ10_ensemblemax = dropdims(maximum(ℰ10, dims = 3), dims = 3)
ℰ10_ensemblemin = dropdims(minimum(ℰ10, dims = 3), dims = 3)
ℰ10_ensemblerange = ℰ10_ensemblemax - ℰ10_ensemblemin

ℰ50 = map(
    ts -> yearatquantile(ts, 0.5),
    view(ℰ, i, j, :, m) for i in 1:size(ℰ,1), j in 1:size(ℰ,2), m in 1:size(ℰ,4)
)
ℰ50_ensemblemean = dropdims(mean(ℰ50, dims = 3), dims = 3)
ℰ50_ensemblemax = dropdims(maximum(ℰ50, dims = 3), dims = 3)
ℰ50_ensemblemin = dropdims(minimum(ℰ50, dims = 3), dims = 3)
ℰ50_ensemblerange = ℰ50_ensemblemax - ℰ50_ensemblemin


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

Gammaout_files = ["/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)/cyclomonth/reemergence_time$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str).nc" for member in members_mean]
Γoutyr3D_ds = open_mfdataset(DimArray(Gammaout_files, Dim{:member}(members_mean)))
Γoutyr3D = readcubedata(Γoutyr3D_ds.adjointage)
Γoutyr3D_timemean = dropdims(mean(Γoutyr3D, dims = :month), dims = :month) # TODO use monthly weights
Γoutyr3D_ensemblemean = dropdims(mean(Γoutyr3D_timemean, dims = 4), dims = 4)
Γoutyr3D_ensemblemax = dropdims(maximum(Γoutyr3D_timemean, dims = 4), dims = 4)
Γoutyr3D_ensemblemin = dropdims(minimum(Γoutyr3D_timemean, dims = 4), dims = 4)
Γoutyr3D_ensemblerange = Γoutyr3D_ensemblemax - Γoutyr3D_ensemblemin

include("plotting_functions.jl")

Γout_ensemblemean = seafloorvalue(Γoutyr3D_ensemblemean, wet3D, gridmetrics)
Γout_ensemblerange = seafloorvalue(Γoutyr3D_ensemblerange, wet3D, gridmetrics)

usecontourf = false

axs = Array{Any,2}(undef, (3, 2))
contours = Array{Any,2}(undef, (3, 2))
nrows, ncols = size(axs)

fig = Figure(size = (ncols * 500, nrows * 250 + 100), fontsize = 18)

yticks = -60:30:60
xticks = -120:60:120 + 360

datamean = (Γout_ensemblemean, ℰ50_ensemblemean, ℰ10_ensemblemean)
datarange = (Γout_ensemblerange, ℰ50_ensemblerange, ℰ10_ensemblerange)
𝒓 = rich("r", font = :bold_italic)
Γstr = rich("Γ", superscript("†"), rich("‾", offset = (-0.55, 0.25)), rich("‾", offset = (-0.85, 0.25)))
Γfun = rich(Γstr, rich("(", 𝒓, ")", offset = (0.4, 0)))
# Qstr = rich("Q", font = :italic)
ℰstr = rich("ℰ", rich("‾", offset = (-0.5, 0.15)))
# Q10 = rich(Qstr, "(0.1)")
# Q50 = rich(Qstr, "(0.5)")
rowlabels = (rich("Mean time, ", Γfun), rich("Median, ", ℰstr, " = 50 %"), rich("10th percentile, ", ℰstr, " = 90 %"))


for (irow, (x2Dmean, x2Drange, text)) in enumerate(zip(datamean, datarange, rowlabels))

    # Plot ensemble mean
    icol = 1
    levels = 0:200:3000
    colormap = cgrad(:viridis, length(levels), categorical = true)
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

    # Plot ensemble range
    icol = 2
    levels = 0:100:1000
    colormap = cgrad(:amp, length(levels), categorical = true)
    highclip = colormap[end]
    colormap = cgrad(colormap[1:end-1], categorical = true)
    colorrange = extrema(levels)

    axs[irow, icol] = ax = Axis(fig[irow, icol]; yticks, xticks, xtickformat, ytickformat)

    contours[irow, icol] = if usecontourf
        plotcontourfmap!(ax, x2Drange, gridmetrics; levels, colormap)
    else
        plotmap!(ax, x2Drange, gridmetrics; colorrange, colormap, highclip) # <- need to fix wrapping longitude for contour levels
    end

    myhidexdecorations!(ax, irow < nrows)
    myhideydecorations!(ax, icol > 1)

    Label(fig[irow, 0]; text, rotation = π/2, tellheight = false)

end


label = rich("ensemble mean (years)")
cb = Colorbar(fig[nrows + 1, 1], contours[1, 1]; label, vertical = false, flipaxis = false, ticks = 0:1000:3000)
cb.width = Relative(2/3)

label = rich("ensemble range, max − min (years)")
cb = Colorbar(fig[nrows + 1, 2], contours[1, 2]; label, vertical = false, flipaxis = false, ticks = 0:200:1000)
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

# Label(fig[0, 1:2]; text = "$(time_window[4:7])s Seafloor Reemergence Time ($(length(members)) members)", fontsize = 24, tellwidth = false)
Label(fig[0, 1:2]; text = "$(time_window[4:7])s Characteristic Timescales of Reemergence ($(length(members)) members)", fontsize = 24, tellwidth = false)
rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)

# save plot
suffix = usecontourf ? "_ctrf" : ""
outputfile = joinpath(outputdir, "reemergencetime_$(time_window)$(suffix).png")
@info "Saving reemergencetime on sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)
outputfile = joinpath(outputdir, "reemergencetime_$(time_window)$(suffix).pdf")
@info "Saving reemergencetime on sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)



# datamean = (Γout_ensemblemean, ℰ50_ensemblemean, ℰ10_ensemblemean)
ikeep = .!isnan.(Γout_ensemblemean) .& (seafloorvalue(Z3D, wet3D) .> 3000)
data = ℰ50_ensemblemean[ikeep] ./ Γout_ensemblemean[ikeep]
weights = Weights(gridmetrics.area2D[ikeep])
fig, ax, plt = hist(data; weights, bins = 0:0.05:3)
outputfile = joinpath(outputdir, "median_over_mean_histogram_$(time_window)$(suffix).png")
@info "Saving reemergencetime on sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)
mean(data, weights)
quantile(data, weights, 0:0.1:1)


ikeep = .!isnan.(Γout_ensemblemean) .& (seafloorvalue(Z3D, wet3D) .> 3000)
data = ℰ10_ensemblemean[ikeep] ./ ℰ50_ensemblemean[ikeep] .- 1
weights = Weights(gridmetrics.area2D[ikeep])
fig, ax, plt = hist(data; weights, bins = -1:0.05:1)
outputfile = joinpath(outputdir, "tenthpercentile_over_median_histogram_$(time_window)$(suffix).png")
@info "Saving reemergencetime on sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)
mean(data, weights)
quantile(data, weights, 0:0.1:1)


ikeep = .!isnan.(Γout_ensemblemean) .& (seafloorvalue(Z3D, wet3D) .> 3000)
data = ℰ50_ensemblemean[ikeep]
weights = Weights(gridmetrics.area2D[ikeep])
fig, ax, plt = hist(data; weights, bins = 0:100:3000)
outputfile = joinpath(outputdir, "median_histogram_$(time_window)$(suffix).png")
@info "Saving reemergencetime on sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)
mean(data, weights)
quantile(data, weights, 0:0.1:1)