# qsub -I -P xv83 -q express -l mem=47GB -l storage=scratch/gh0+scratch/xv83 -l walltime=01:00:00 -l ncpus=12

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

time_window1 = "Jan2030-Dec2039"
time_window2 = "Jan2090-Dec2099"
experiment1 = parse(Int, time_window1[4:7]) ≤ 2010 ? "historical" : "ssp370"
experiment2 = parse(Int, time_window2[4:7]) ≤ 2010 ? "historical" : "ssp370"

members_mean = ["r$(r)i1p1f1" for r in 1:40]
members = ["r$(r)i1p1f1" for r in 1:40]

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



# Preferred diffusivities
κVdeep = 3.0e-5 # m^2/s
κVML = 1.0      # m^2/s
κH = 300.0      # m^2/s
κVdeep_str = "kVdeep" * format(κVdeep, conversion="e")
κVML_str = "kVML" * format(κVML, conversion="e")
κH_str = "kH" * format(κH, conversion="d")
upwind = false
upwind_str = upwind ? "" : "_centered"
upwind_str2 = upwind ? "upwind" : "centered"

# Use yearly time-stepped simulations or monthly ones?
yearly = true
yearly_str = yearly ? "_yearly" : ""
yearly_str2 = yearly ? "(yearly)" : ""

# little helper function to get the year given a quantile
function yearatquantile(ℰ, ℰlevel)
    isnan(ℰ[1]) && return NaN
    out = findfirst(ℰ .< ℰlevel)
    isnothing(out) ? maximum(years) : Float64(out)
end

# To avoid loading and carrying 100s of GB of data around,
# preprocess each member before and only save the 2D data needed for plots.
if yearly
    ℰ_file0 = "/scratch/xv83/TMIP/data/$model/$experiment1/$(first(members))/$(time_window1)/seqeff$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str).nc"
    ℰ_ds0 = open_dataset(ℰ_file0)
    years = ℰ_ds0.Ti |> Array
    ℰ1050_ensemble1 = reduce((x,y) -> cat(x, y, dims = 4), map(members) do member
        @info "loading $member ℰ"
        ℰ_file = "/scratch/xv83/TMIP/data/$model/$experiment1/$member/$(time_window1)/seqeff$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str).nc"
        ℰ_ds = open_dataset(ℰ_file)
        ℰ = readcubedata(ℰ_ds.seqeff)
        ℰ10 = map(
            ts -> yearatquantile(ts, 0.9),
            view(ℰ, i, j, :) for i in 1:size(ℰ,1), j in 1:size(ℰ,2)
        )
        ℰ50 = map(
            ts -> yearatquantile(ts, 0.5),
            view(ℰ, i, j, :) for i in 1:size(ℰ,1), j in 1:size(ℰ,2)
        )
        [ℰ10;;; ℰ50]
    end)
    ℰ1050_ensemble2 = reduce((x,y) -> cat(x, y, dims = 4), map(members) do member
        @info "loading $member ℰ"
        ℰ_file = "/scratch/xv83/TMIP/data/$model/$experiment2/$member/$(time_window2)/seqeff$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str).nc"
        ℰ_ds = open_dataset(ℰ_file)
        ℰ = readcubedata(ℰ_ds.seqeff)
        ℰ10 = map(
            ts -> yearatquantile(ts, 0.9),
            view(ℰ, i, j, :) for i in 1:size(ℰ,1), j in 1:size(ℰ,2)
        )
        ℰ50 = map(
            ts -> yearatquantile(ts, 0.5),
            view(ℰ, i, j, :) for i in 1:size(ℰ,1), j in 1:size(ℰ,2)
        )
        [ℰ10;;; ℰ50]
    end)
    ℰ1050_ensemblemean1 = dropdims(mean(ℰ1050_ensemble1, dims = 4), dims = 4)
    ℰ1050_ensemblemean2 = dropdims(mean(ℰ1050_ensemble2, dims = 4), dims = 4)
    ℰ1050_ensemblemean_diff = ℰ1050_ensemblemean2 - ℰ1050_ensemblemean1


else
    # TODO (don't forget to deal with months)
end

include("plotting_functions.jl") # load seafloorvalue function

Γout_ensemble1 = reduce((x, y) -> cat(x, y, dims = 3), map(members) do member
    @info "loading $member Γ†"
    Γout_file = "/scratch/xv83/TMIP/data/$model/$experiment1/$member/$(time_window1)/cyclomonth/reemergence_time$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str).nc"
    Γoutyr4D_ds = open_dataset(Γout_file)
    Γoutyr3D = dropdims(mean(readcubedata(Γoutyr4D_ds.adjointage), dims = Ti), dims = Ti)
    seafloorvalue(Γoutyr3D, wet3D, gridmetrics)
end)
Γout_ensemble2 = reduce((x, y) -> cat(x, y, dims = 3), map(members) do member
    @info "loading $member Γ†"
    Γout_file = "/scratch/xv83/TMIP/data/$model/$experiment2/$member/$(time_window2)/cyclomonth/reemergence_time$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str).nc"
    Γoutyr4D_ds = open_dataset(Γout_file)
    Γoutyr3D = dropdims(mean(readcubedata(Γoutyr4D_ds.adjointage), dims = Ti), dims = Ti)
    seafloorvalue(Γoutyr3D, wet3D, gridmetrics)
end)

Γout_ensemblemean1 = dropdims(mean(Γout_ensemble1, dims = 3), dims = 3)
Γout_ensemblemean2 = dropdims(mean(Γout_ensemble2, dims = 3), dims = 3)
Γout_ensemblemean_diff = Γout_ensemblemean2 - Γout_ensemblemean1




usecontourf = false

axs = Array{Any,2}(undef, (3, 2))
contours = Array{Any,2}(undef, (3, 2))
nrows, ncols = size(axs)

fig = Figure(size = (ncols * 500, nrows * 250 + 100), fontsize = 18)

yticks = -60:30:60
xticks = -120:60:120 + 360

datamean = (Γout_ensemblemean2, ℰ1050_ensemblemean2[:,:,2], ℰ1050_ensemblemean2[:,:,1])
datadiff = (Γout_ensemblemean_diff, ℰ1050_ensemblemean_diff[:,:,2], ℰ1050_ensemblemean_diff[:,:,1])
𝒓 = rich("r", font = :bold_italic)
Γstr = rich("Γ", superscript("†"), rich("‾", offset = (-0.55, 0.25)), rich("‾", offset = (-0.85, 0.25)))
Γfun = rich(Γstr, rich("(", 𝒓, ")", offset = (0.4, 0)))
# Qstr = rich("Q", font = :italic)
ℰstr = rich("ℰ", rich("‾", offset = (-0.5, 0.15)))
# Q10 = rich(Qstr, "(0.1)")
# Q50 = rich(Qstr, "(0.5)")
rowlabels = (rich("Mean time, ", Γfun), rich("Median, ", ℰstr, " = 50 %"), rich("10th percentile, ", ℰstr, " = 90 %"))


for (irow, (x2Dmean, x2Ddiff, text)) in enumerate(zip(datamean, datadiff, rowlabels))

    # Plot ensemble mean
    icol = 1
    levels = 0:200:4000
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

    # Plot ensemble diff
    icol = 2
    levels = -1200:100:1200
    colormap = cgrad(cgrad(:balance, length(levels), categorical = true)[[1:end÷2+1; end÷2+1:end]], categorical = true)
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
cb = Colorbar(fig[nrows + 1, 1], contours[1, 1]; label, vertical = false, flipaxis = false, ticks = 0:1000:4000)
cb.width = Relative(2/3)

label = rich("ensemble mean 2090s − 2030s (years)")
cb = Colorbar(fig[nrows + 1, 2], contours[1, 2]; label, vertical = false, flipaxis = false, ticks = -1200:600:1200, tickformat = divergingcbarticklabelformat)
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

Label(fig[0, 1:2]; text = "Climate Change Effect on Seafloor Reemergence Time ($(length(members)) members)$(yearly_str2)", fontsize = 24, tellwidth = false)
rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)

# save plot
suffix = usecontourf ? "_ctrf" : ""
outputfile = joinpath(outputdir, "reemergencetime_diff$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str)_$(time_window2)$(suffix).png")
@info "Saving reemergencetime on sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)
outputfile = joinpath(outputdir, "reemergencetime_diff$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str)_$(time_window2)$(suffix).pdf")
@info "Saving reemergencetime on sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)









ikeep = .!isnan.(datamean[1]) .& (seafloorvalue(Z3D, wet3D) .> 3000)
data = Γout_ensemblemean_diff[ikeep] ./ Γout_ensemblemean1[ikeep]
weights = Weights(gridmetrics.area2D[ikeep])
fig, ax, plt = hist(data; weights, bins = -1:0.1:3)
outputfile = joinpath(outputdir, "diffmeanage_over_meanage_$(time_window1)-$(time_window2).png")
@info "Saving as image file:\n  $(outputfile)"
save(outputfile, fig)
@show mean(data, weights)
@show std(data, weights)
@show quantile(data, weights, 0:0.1:1)




ikeep = .!isnan.(datamean[1]) .& (seafloorvalue(Z3D, wet3D) .> 3000)
data = ℰ1050_ensemblemean_diff[:,:,2][ikeep] ./ ℰ1050_ensemblemean1[:,:,2][ikeep]
weights = Weights(gridmetrics.area2D[ikeep])
fig, ax, plt = hist(data; weights, bins = -1:0.1:3)
outputfile = joinpath(outputdir, "diffE50_over_E50_$(time_window1)-$(time_window2).png")
@info "Saving as image file:\n  $(outputfile)"
save(outputfile, fig)
@show mean(data, weights)
@show std(data, weights)
@show quantile(data, weights, 0:0.1:1)



ikeep = .!isnan.(datamean[1]) .& (seafloorvalue(Z3D, wet3D) .> 3000)
data = ℰ1050_ensemblemean_diff[:,:,1][ikeep] ./ ℰ1050_ensemblemean1[:,:,1][ikeep]
weights = Weights(gridmetrics.area2D[ikeep])
fig, ax, plt = hist(data; weights, bins = -1:0.1:3)
outputfile = joinpath(outputdir, "diffE90_over_E90_$(time_window1)-$(time_window2).png")
@info "Saving as image file:\n  $(outputfile)"
save(outputfile, fig)
@show mean(data, weights)
@show std(data, weights)
@show quantile(data, weights, 0:0.1:1)

