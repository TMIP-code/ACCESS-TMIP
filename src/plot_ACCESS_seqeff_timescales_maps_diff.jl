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
yearly = false
yearly_str = yearly ? "_yearly" : ""
yearly_str2 = yearly ? "(yearly)" : ""

τs = [100, 300, 1000]

# little helper function to get the year given a quantile
function yearatquantile(ℰ, ℰlevel)
    isnan(ℰ[1]) && return NaN
    out = findfirst(ℰ .< ℰlevel)
    isnothing(out) ? maximum(years) : Float64(out)
end

# To avoid loading and carrying 100s of GB of data around,
# preprocess each member before and only save the 2D data needed for plots.
varname = yearly ? "seqeff" : "calE"
ℰ_file0 = "/scratch/xv83/TMIP/data/$model/$experiment1/$(first(members))/$(time_window1)/$(varname)$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str).nc"
ℰ_ds0 = open_dataset(ℰ_file0)
years = ℰ_ds0.Ti |> Array
τℰ1050_ensemble1 = reduce((x,y) -> cat(x, y, dims = 4), map(members) do member
    @info "loading $member ℰ"
    ℰ_file = "/scratch/xv83/TMIP/data/$model/$experiment1/$member/$(time_window1)/$(varname)$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str).nc"
    ℰ_ds = open_dataset(ℰ_file)
    ℰ = readcubedata(ℰ_ds[varname])
    τℰ10 = map(
        ts -> yearatquantile(ts, 0.9),
        view(ℰ, i, j, :) for i in 1:size(ℰ,1), j in 1:size(ℰ,2)
    )
    τℰ50 = map(
        ts -> yearatquantile(ts, 0.5),
        view(ℰ, i, j, :) for i in 1:size(ℰ,1), j in 1:size(ℰ,2)
    )
    [τℰ10;;; τℰ50]
end)
τℰ1050_ensemble2 = reduce((x,y) -> cat(x, y, dims = 4), map(members) do member
    @info "loading $member ℰ"
    ℰ_file = "/scratch/xv83/TMIP/data/$model/$experiment2/$member/$(time_window2)/$(varname)$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str).nc"
    ℰ_ds = open_dataset(ℰ_file)
    ℰ = readcubedata(ℰ_ds[varname])
    τℰ10 = map(
        ts -> yearatquantile(ts, 0.9),
        view(ℰ, i, j, :) for i in 1:size(ℰ,1), j in 1:size(ℰ,2)
    )
    τℰ50 = map(
        ts -> yearatquantile(ts, 0.5),
        view(ℰ, i, j, :) for i in 1:size(ℰ,1), j in 1:size(ℰ,2)
    )
    [τℰ10;;; τℰ50]
end)
τℰ1050_ensemblemean1 = dropdims(mean(τℰ1050_ensemble1, dims = 4), dims = 4)
τℰ1050_ensemblemean2 = dropdims(mean(τℰ1050_ensemble2, dims = 4), dims = 4)
τℰ1050_ensemblemean_diff = τℰ1050_ensemblemean2 - τℰ1050_ensemblemean1


# Read 2030s files again for seq eff values now
ℰ_files = ["/scratch/xv83/TMIP/data/$model/$experiment1/$member/$(time_window1)/$(varname)$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str).nc" for member in members]
ℰ_ds = open_mfdataset(DimArray(ℰ_files, Dim{:member}(members)))
# Load the data for the selected timescales only (to avoid using too much memory)
ℰ = readcubedata(ℰ_ds[varname][Ti = At(τs)])
years = ℰ_ds.Ti |> Array
ℰ_ensemblemean1 = dropdims(mean(ℰ, dims = :member), dims = :member)
# Read 2090s files
ℰ_files = ["/scratch/xv83/TMIP/data/$model/$experiment2/$member/$(time_window2)/$(varname)$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str).nc" for member in members]
ℰ_ds = open_mfdataset(DimArray(ℰ_files, Dim{:member}(members)))
ℰ = readcubedata(ℰ_ds[varname][Ti = At(τs)])
ℰ_ensemblemean2 = dropdims(mean(ℰ, dims = :member), dims = :member)
ℰ_diff = ℰ_ensemblemean2 - ℰ_ensemblemean1


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

axs = Array{Any,2}(undef, (2, 2))
contours = Array{Any,2}(undef, (2, 2))
nrows, ncols = size(axs)

fig = Figure(size = (ncols * 500, nrows * 250 + 200), fontsize = 18)

yticks = -60:30:60
xticks = -120:60:120 + 360

# datamean = (Γout_ensemblemean2, τℰ1050_ensemblemean2[:,:,2], τℰ1050_ensemblemean2[:,:,1])
# datadiff = (Γout_ensemblemean_diff, τℰ1050_ensemblemean_diff[:,:,2], τℰ1050_ensemblemean_diff[:,:,1])
datamean = (Γout_ensemblemean2, 100 * ℰ_ensemblemean2[:, :, 3])
datadiff = (Γout_ensemblemean_diff, ℰ_diff[:,:,3])

𝒓 = rich("r", font = :bold_italic)
Γstr = rich("Γ", superscript("†"), rich("‾", offset = (-0.55, 0.25)), rich("‾", offset = (-0.85, 0.25)))
Γfun = rich(Γstr, rich("(", 𝒓, ")", offset = (0.4, 0)))
# Qstr = rich("Q", font = :italic)
ℰstr = rich("ℰ", rich("‾", offset = (-0.5, 0.15)))
# τstr = rich("τ", font = :italic)
# ℰfun = rich(ℰstr, "(", 𝒓, ", τ)")
ℰfun = rich(ℰstr, "(τ)")
# Q10 = rich(Qstr, "(0.1)")
# Q50 = rich(Qstr, "(0.5)")
# rowlabels = (rich("Mean time, ", Γfun), rich(ℰstr, "(τ = $year years)"))
# levels = (0:200:4000, 0:10:100)
# levels_diff = (-1200:100:1200, -50:10:50)
# cmaps = (:viridis, :balance)



# a 2090s seqeff
irow, icol = 1, 1
levels = 0:10:100
colormap = cgrad(:Zissou1Continuous, length(levels) - 1, categorical = true, rev = true)
highclip = colormap[end]
colorrange = extrema(levels)
axs[irow, icol] = ax = Axis(fig[irow + 1, icol]; yticks, xticks, xtickformat, ytickformat, aspect = DataAspect())
co = plotmap!(ax, 100 * ℰ_ensemblemean2[:, :, 3], gridmetrics; colorrange, colormap, highclip)
myhidexdecorations!(ax, irow < nrows)
myhideydecorations!(ax, icol > 1)
year = τs[3]
Label(fig[irow + 1, 0]; text = rich(ℰstr, "(τ = $year years)"), rotation = π/2, tellheight = false)
label = rich("ensemble mean 2090s ", ℰfun, " (%)")
cb = Colorbar(fig[irow, icol], co; label, vertical = false, flipaxis = true, ticks = 0:20:100)
cb.width = Relative(2/3)

# b 2090s ensemble mean reemergence time
irow, icol = 2, 1
levels = 0:200:4000
colormap = cgrad(:viridis, length(levels), categorical = true)
highclip = colormap[end]
colormap = cgrad(colormap[1:end-1], categorical = true)
colorrange = extrema(levels)
axs[irow, icol] = ax = Axis(fig[irow + 1, icol]; yticks, xticks, xtickformat, ytickformat, aspect = DataAspect())
co = plotmap!(ax, Γout_ensemblemean2, gridmetrics; colorrange, colormap, highclip)
myhidexdecorations!(ax, irow < nrows)
myhideydecorations!(ax, icol > 1)
Label(fig[irow + 1, 0]; text = rich("Mean time, ", Γstr), rotation = π/2, tellheight = false)
label = rich("ensemble mean 2090s ", Γstr, "  (years)")
cb = Colorbar(fig[irow + 2, icol], co; label, vertical = false, flipaxis = false, ticks = 0:1000:4000)
cb.width = Relative(2/3)



# c seqeff diff
irow, icol = 1, 2
levels = -50:10:50
colormap = cgrad(cgrad(:tol_bu_rd, length(levels), categorical = true)[[1:end÷2+1; end÷2+1:end]], categorical = true)
highclip = colormap[end]
lowclip = colormap[1]
colormap = cgrad(colormap[2:end-1], categorical = true)
colorrange = extrema(levels)
axs[irow, icol] = ax = Axis(fig[irow + 1, icol]; yticks, xticks, xtickformat, ytickformat, aspect = DataAspect())
co = plotmap!(ax, 100 * ℰ_diff[:,:,3], gridmetrics; colorrange, colormap, highclip, lowclip)
myhidexdecorations!(ax, irow < nrows)
myhideydecorations!(ax, icol > 1)
label = rich("ensemble mean 2090s − 2030s ", ℰfun, " (%)")
cb = Colorbar(fig[irow, icol], co; label, vertical = false, flipaxis = true, ticks = -50:20:50, tickformat = divergingcbarticklabelformat)
cb.width = Relative(2/3)


# d ensemble mean reemergence time diff
irow, icol = 2, 2
levels = -1200:100:1200
colormap = cgrad(cgrad(:balance, length(levels), categorical = true)[[1:end÷2+1; end÷2+1:end]], categorical = true)
highclip = colormap[end]
lowclip = colormap[1]
colormap = cgrad(colormap[2:end-1], categorical = true)
colorrange = extrema(levels)
axs[irow, icol] = ax = Axis(fig[irow + 1, icol]; yticks, xticks, xtickformat, ytickformat, aspect = DataAspect())
co = plotmap!(ax, Γout_ensemblemean_diff, gridmetrics; colorrange, colormap, highclip, lowclip)
myhidexdecorations!(ax, irow < nrows)
myhideydecorations!(ax, icol > 1)
label = rich("ensemble mean 2090s − 2030s ", Γstr, "  (years)")
cb = Colorbar(fig[irow + 2, icol], co; label, vertical = false, flipaxis = false, ticks = -1200:600:1200, tickformat = divergingcbarticklabelformat)
cb.width = Relative(2/3)





# column labels
# Label(fig[0, 1]; text = "ensemble mean", tellwidth = false)
# Label(fig[0, 2]; text = "ensemble range (internal variability)", tellwidth = false)

labels = [
    "a" "c"
    "b" "d"
]

labeloptions = (
    font = :bold,
    align = (:left, :bottom),
    offset = (5, 2),
    space = :relative,
    fontsize = 24
)

for (ax, label) in zip(axs, labels)
    txt = text!(ax, 0, 0; text = "$label", labeloptions..., strokecolor = :white, strokewidth = 3)
    translate!(txt, 0, 0, 100)
    txt = text!(ax, 0, 0; text = "$label", labeloptions...)
    translate!(txt, 0, 0, 100)
end

# Label(fig[0, 1:2]; text = "Climate Change Effect on Seafloor Reemergence Time ($(length(members)) members)$(yearly_str2)", fontsize = 24, tellwidth = false)
rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)

colsize!(fig.layout, 1, Aspect(2, 2.0))
colsize!(fig.layout, 2, Aspect(2, 2.0))
resize_to_layout!(fig)


# save plot
suffix = usecontourf ? "_ctrf" : ""
outputfile = joinpath(outputdir, "seqeff_timescales_diff$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str)_$(time_window2)$(suffix).png")
@info "Saving as image file:\n  $(outputfile)"
save(outputfile, fig)
outputfile = joinpath(outputdir, "seqeff_timescales_diff$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str)_$(time_window2)$(suffix).pdf")
@info "Saving as image file:\n  $(outputfile)"
save(outputfile, fig)












# Now same plot but for the SI figure

axs = Array{Any,2}(undef, (4, 2))
contours = Array{Any,2}(undef, (4, 2))
nrows, ncols = size(axs)

fig = Figure(size = (ncols * 500, nrows * 250 + 200), fontsize = 18)

yticks = -60:30:60
xticks = -120:60:120 + 360

# datamean = (Γout_ensemblemean2, τℰ1050_ensemblemean2[:,:,2], τℰ1050_ensemblemean2[:,:,1])
# datadiff = (Γout_ensemblemean_diff, τℰ1050_ensemblemean_diff[:,:,2], τℰ1050_ensemblemean_diff[:,:,1])
datamean = (Γout_ensemblemean2, 100 * ℰ_ensemblemean2[:, :, 3])
datadiff = (Γout_ensemblemean_diff, ℰ_diff[:,:,3])

𝒓 = rich("r", font = :bold_italic)
Γstr = rich("Γ", superscript("†"), rich("‾", offset = (-0.55, 0.25)), rich("‾", offset = (-0.85, 0.25)))
Γfun = rich(Γstr, rich("(", 𝒓, ")", offset = (0.4, 0)))
# Qstr = rich("Q", font = :italic)
ℰstr = rich("ℰ", rich("‾", offset = (-0.5, 0.15)))
# τstr = rich("τ", font = :italic)
# ℰfun = rich(ℰstr, "(", 𝒓, ", τ)")
ℰfun = rich(ℰstr, "(τ)")
# Q10 = rich(Qstr, "(0.1)")
# Q50 = rich(Qstr, "(0.5)")
# rowlabels = (rich("Mean time, ", Γfun), rich(ℰstr, "(τ = $year years)"))
# levels = (0:200:4000, 0:10:100)
# levels_diff = (-1200:100:1200, -50:10:50)
# cmaps = (:viridis, :balance)


# a,b Plot 2090s seqeff for 100 and 300 years
levels = 0:10:100
colormap = cgrad(:Zissou1Continuous, length(levels) - 1, categorical = true, rev = true)
highclip = colormap[end]
colorrange = extrema(levels)

irow, icol = 1, 1
axs[irow, icol] = ax = Axis(fig[irow + 1, icol]; yticks, xticks, xtickformat, ytickformat, aspect = DataAspect())
co = plotmap!(ax, 100 * ℰ_ensemblemean2[:, :, 2], gridmetrics; colorrange, colormap, highclip)
myhidexdecorations!(ax, irow < nrows)
myhideydecorations!(ax, icol > 1)
year = τs[2]
Label(fig[irow + 1, 0]; text = rich(ℰstr, "(τ = $year years)"), rotation = π/2, tellheight = false)

label = rich("ensemble mean 2090s ", ℰfun, " (%)")
cb = Colorbar(fig[irow, icol], co; label, vertical = false, flipaxis = true, ticks = 0:20:100)
cb.width = Relative(2/3)

irow, icol = 2, 1
axs[irow, icol] = ax = Axis(fig[irow + 1, icol]; yticks, xticks, xtickformat, ytickformat, aspect = DataAspect())
co = plotmap!(ax, 100 * ℰ_ensemblemean2[:, :, 1], gridmetrics; colorrange, colormap, highclip)
myhidexdecorations!(ax, irow < nrows)
myhideydecorations!(ax, icol > 1)
year = τs[1]
Label(fig[irow + 1, 0]; text = rich(ℰstr, "(τ = $year years)"), rotation = π/2, tellheight = false)


# c,d Plots 2090s ensemble mean median and 10th percentile time
levels = 0:200:4000
colormap = cgrad(:viridis, length(levels), categorical = true)
highclip = colormap[end]
colormap = cgrad(colormap[1:end-1], categorical = true)
colorrange = extrema(levels)

irow, icol = 3, 1
axs[irow, icol] = ax = Axis(fig[irow + 1, icol]; yticks, xticks, xtickformat, ytickformat, aspect = DataAspect())
co = plotmap!(ax, τℰ1050_ensemblemean2[:,:,2], gridmetrics; colorrange, colormap, highclip)
myhidexdecorations!(ax, irow < nrows)
myhideydecorations!(ax, icol > 1)
Label(fig[irow + 1, 0]; text = rich("Median time (", ℰstr, " = 50 %)"), rotation = π/2, tellheight = false)


irow, icol = 4, 1
axs[irow, icol] = ax = Axis(fig[irow + 1, icol]; yticks, xticks, xtickformat, ytickformat, aspect = DataAspect())
co = plotmap!(ax, τℰ1050_ensemblemean2[:,:,1], gridmetrics; colorrange, colormap, highclip)
myhidexdecorations!(ax, irow < nrows)
myhideydecorations!(ax, icol > 1)
Label(fig[irow + 1, 0]; text = rich("10th %ile time (", ℰstr, " = 90 %)"), rotation = π/2, tellheight = false)

label = rich("ensemble mean 2090s characteristic timescales (years)")
cb = Colorbar(fig[irow + 2, icol], co; label, vertical = false, flipaxis = false, ticks = 0:1000:4000)
cb.width = Relative(2/3)


# e,f Plot seqeff diff
levels = -50:10:50
colormap = cgrad(cgrad(:tol_bu_rd, length(levels), categorical = true)[[1:end÷2+1; end÷2+1:end]], categorical = true)
highclip = colormap[end]
lowclip = colormap[1]
colormap = cgrad(colormap[2:end-1], categorical = true)
colorrange = extrema(levels)

irow, icol = 1, 2
axs[irow, icol] = ax = Axis(fig[irow + 1, icol]; yticks, xticks, xtickformat, ytickformat, aspect = DataAspect())
co = plotmap!(ax, 100 * ℰ_diff[:,:,2], gridmetrics; colorrange, colormap, highclip, lowclip)
myhidexdecorations!(ax, irow < nrows)
myhideydecorations!(ax, icol > 1)

label = rich("ensemble mean 2090s − 2030s ", ℰfun, " (%)")
cb = Colorbar(fig[irow, icol], co; label, vertical = false, flipaxis = true, ticks = -50:20:50, tickformat = divergingcbarticklabelformat)
cb.width = Relative(2/3)

irow, icol = 2, 2
axs[irow, icol] = ax = Axis(fig[irow + 1, icol]; yticks, xticks, xtickformat, ytickformat, aspect = DataAspect())
co = plotmap!(ax, 100 * ℰ_diff[:,:,1], gridmetrics; colorrange, colormap, highclip, lowclip)
myhidexdecorations!(ax, irow < nrows)
myhideydecorations!(ax, icol > 1)


# g,h Plot ensemble mean median and 10th percentile time diff
levels = -1200:100:1200
colormap = cgrad(cgrad(:balance, length(levels), categorical = true)[[1:end÷2+1; end÷2+1:end]], categorical = true)
highclip = colormap[end]
lowclip = colormap[1]
colormap = cgrad(colormap[2:end-1], categorical = true)
colorrange = extrema(levels)

irow, icol = 3, 2
axs[irow, icol] = ax = Axis(fig[irow + 1, icol]; yticks, xticks, xtickformat, ytickformat, aspect = DataAspect())
co = plotmap!(ax, τℰ1050_ensemblemean_diff[:, :, 2], gridmetrics; colorrange, colormap, highclip, lowclip)
myhidexdecorations!(ax, irow < nrows)
myhideydecorations!(ax, icol > 1)

irow, icol = 4, 2
axs[irow, icol] = ax = Axis(fig[irow + 1, icol]; yticks, xticks, xtickformat, ytickformat, aspect = DataAspect())
co = plotmap!(ax, τℰ1050_ensemblemean_diff[:, :, 1], gridmetrics; colorrange, colormap, highclip, lowclip)
myhidexdecorations!(ax, irow < nrows)
myhideydecorations!(ax, icol > 1)

label = rich("ensemble mean 2090s − 2030s  (years)")
cb = Colorbar(fig[irow + 2, icol], co; label, vertical = false, flipaxis = false, ticks = -1200:600:1200, tickformat = divergingcbarticklabelformat)
cb.width = Relative(2/3)



# column labels
# Label(fig[0, 1]; text = "ensemble mean", tellwidth = false)
# Label(fig[0, 2]; text = "ensemble range (internal variability)", tellwidth = false)

labels = [
    "a" "e"
    "b" "f"
    "c" "g"
    "d" "h"
]

labeloptions = (
    font = :bold,
    align = (:left, :bottom),
    offset = (5, 2),
    space = :relative,
    fontsize = 24
)

for (ax, label) in zip(axs, labels)
    txt = text!(ax, 0, 0; text = "$label", labeloptions..., strokecolor = :white, strokewidth = 3)
    translate!(txt, 0, 0, 100)
    txt = text!(ax, 0, 0; text = "$label", labeloptions...)
    translate!(txt, 0, 0, 100)
end

# Label(fig[0, 1:2]; text = "Climate Change Effect on Seafloor Reemergence Time ($(length(members)) members)$(yearly_str2)", fontsize = 24, tellwidth = false)
rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)
rowgap!(fig.layout, 3, 20)
colgap!(fig.layout, 2, 20)

colsize!(fig.layout, 1, Aspect(2, 2.0))
colsize!(fig.layout, 2, Aspect(2, 2.0))
resize_to_layout!(fig)


# save plot
suffix = usecontourf ? "_ctrf" : ""
outputfile = joinpath(outputdir, "SI_seqeff_timescales_diff$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str)_$(time_window2)$(suffix).png")
@info "Saving as image file:\n  $(outputfile)"
save(outputfile, fig)
outputfile = joinpath(outputdir, "SI_seqeff_timescales_diff$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str)_$(time_window2)$(suffix).pdf")
@info "Saving as image file:\n  $(outputfile)"
save(outputfile, fig)









# Save the data to be uploaded with paper


axlist = (
    dims(areacello_ds["areacello"])...,
    dims(DimArray(ones(2), Dim{:climatology}(["2030s", "2090s"])))[1],
)

metadata = Dict(
    "description" => "Sequestration efficiency",
    "unit" => "",
    "Ti unit" => "years",
)
seqeff3D = DimensionalData.rebuild(areacello_ds["areacello"];
    data = permutedims([ℰ_ensemblemean1;;;; ℰ_ensemblemean2], (1, 2, 4, 3)) .|> Float64,
    dims = (axlist..., dims(ℰ_ensemblemean1)[3]),
    metadata = metadata,
)

metadata = Dict(
    "description" => "Characteristic timescales",
    "unit" => "years",
)
timescales3D = DimensionalData.rebuild(areacello_ds["areacello"];
    data = permutedims([Γout_ensemblemean1;;; τℰ1050_ensemblemean1;;;; Γout_ensemblemean2;;; τℰ1050_ensemblemean2], (1, 2, 4, 3)) .|> Float64,
    dims = (axlist..., dims(DimArray(ones(3), Dim{:timescale}(["mean", "median", "10th percentile"])))[1]),
    metadata = metadata,
)

properties = Dict(
    "description" => "Sequestration efficiencies and characteristic timescales as plotted in Figs. 4 and S6 in Pasquier et al. (2025)",
    "model" => model,
    "experiment" => experiment1,
)
arrays = Dict(:seqeff => seqeff3D, :timescales => timescales3D, :lat => readcubedata(volcello_ds.lat), :lon => readcubedata(volcello_ds.lon))
ds = Dataset(; properties = properties, arrays...)

# Save to netCDF file
outputfile = joinpath(outputdir, "Pasquier_etal_GRL_2025_Fig4_data.nc")
@info "Saving climate-change figure data as netCDF file:\n  $(outputfile)"
# ds_chunked = setchunks(ds, (x = 60, y = 60, Ti = length(ds.Ti)))
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)

