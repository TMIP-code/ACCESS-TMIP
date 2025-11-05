# qsub -I -P xv83 -q express -l mem=47GB -l storage=scratch/gh0+scratch/xv83 -l walltime=00:30:00 -l ncpus=12

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
# using LaTeXStrings

# Load some stuff from my ACCESS runs (for the grid, lat, lon, depth etc.)
model = "ACCESS-ESM1-5"

time_window = "Jan2030-Dec2039"
# time_window = "Jan2090-Dec2099"
experiment = parse(Int, time_window[4:7]) ≤ 2010 ? "historical" : "ssp370"

# Gadi directory for input files
# inputdirfun(member) = "/scratch/xv83/TMIP/data/$model/$experiment/all members/$(time_window)"
inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/all_members/$(time_window)/cyclomonth"
outputdir = "/scratch/xv83/TMIP/AndersonAcceleration"
mkpath(outputdir)


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

# Directory where the Anderson Accelerated age fields are
age_dir = "/scratch/xv83/bp3051/access-esm/archive/andersonacceleration_test-n10-5415f621/age_output"
age_files = [joinpath(age_dir, f) for f in readdir(age_dir) if startswith(f, "ocean_age.res")]
age_files = sort(age_files)
age_idx = [parse(Int, f[(end - 6):(end - 3)]) for f in age_files]
age_ds = open_mfdataset(DimArray(age_files, Dim{:Decade}(age_idx)))
age = readcubedata(age_ds.age_global)
age = dropdims(age, dims = 4)
age4D = age.data
age4D[.!wet3D, :] .= NaN;

# Load ages at start of cycles for drift plot
cyclestart_age_dir = "/scratch/xv83/bp3051/access-esm/archive/andersonacceleration_test-n10-5415f621/age_output"
cyclestart_age_files = [joinpath(cyclestart_age_dir, f) for f in readdir(cyclestart_age_dir) if startswith(f, "ocean_age.cyclestart")]
cyclestart_age_files = sort(cyclestart_age_files)
cyclestart_idx = [parse(Int, f[(end - 6):(end - 3)]) for f in cyclestart_age_files] # CHECK shift?
cyclestart_age_ds = open_mfdataset(DimArray(cyclestart_age_files, Dim{:Decade}(cyclestart_idx)))
cyclestart_age = readcubedata(cyclestart_age_ds.age_global)
cyclestart_age = dropdims(cyclestart_age, dims = 4)
cyclestart_age4D = cyclestart_age.data
cyclestart_age4D[.!wet3D, :] .= NaN;


maxage = dropdims(nanmaximum(age4D, dims = (1, 2, 3)), dims = (1, 2, 3))
i2000m = argmin(abs.(zt .- 2000))
volumeintegralage2000m = dropdims(nansum((age4D .* v3D)[:, :, i2000m, :], dims = (1, 2)), dims = (1, 2))
volumeintergal2000m = dropdims(nansum(v3D[:, :, i2000m], dims = (1, 2)), dims = (1, 2))
meanage2000m = volumeintegralage2000m ./ volumeintergal2000m
i4000m = argmin(abs.(zt .- 4000))
volumeintegralage4000m = dropdims(nansum((age4D .* v3D)[:, :, i4000m, :], dims = (1, 2)), dims = (1, 2))
volumeintergal4000m = dropdims(nansum(v3D[:, :, i4000m], dims = (1, 2)), dims = (1, 2))
meanage4000m = volumeintegralage4000m ./ volumeintergal4000m

include("plotting_functions.jl")


fig = Figure()
axisoptions = (
    limits = (0, nothing, nothing, nothing),
    xlabel = "ACCESS-ESM1.5 simulation years",
)
ax1 = Axis(fig[1, 1]; axisoptions..., ylabel = "age (years)")
x = 10 * age_idx
lines!(ax1, x, maxage; label = "maximum")
lines!(ax1, x, meanage2000m; label = "2000m mean")
lines!(ax1, x, meanage4000m; label = "4000m mean")
ylims!(ax1, 0, nothing)
axislegend(ax1; position = :lt, framevisible = false)
hidexdecorations!(ax1, grid = false)

ax2 = Axis(fig[2, 1]; axisoptions..., ylabel = "10-yr cycle Δage (years)")
x = 10 * age_idx[2:end]
lines!(ax2, x, diff(maxage); label = "maximum")
lines!(ax2, x, diff(meanage2000m); label = "2000m mean")
lines!(ax2, x, diff(meanage4000m); label = "4000m mean")

linkxaxes!(ax2, ax1)

outputfile = joinpath(outputdir, "AA_age_timeseries_postmaxagefix.png")
@info "Saving image file:\n  $(outputfile)"
save(outputfile, fig)


# Rebuild dims from volcello
finalage = rebuild(
    age;
    data = age4D[:, :, :, end],
    dims = dims(volcello),
)
finalage2000m = finalage[lev = Near(2000)]
finalage4000m = finalage[lev = Near(4000)]

fig = Figure()
Ncycles = length(age_idx) - 1
Label(fig[0, :], "Anderson-Acceleration age after $(10Ncycles) simulation years ($(Ncycles) cycles)", tellwidth = false)
yticks = -60:30:60
xticks = -120:60:(120 + 360)
levels = 0:100:2500
colormap = cgrad(:viridis, length(levels), categorical = true)
highclip = colormap[end]
colormap = cgrad(colormap[1:(end - 1)], categorical = true)
lowclip = :red
colorrange = extrema(levels)


# 2000m age
ax = Axis(fig[1, 1]; yticks, xticks, xtickformat, ytickformat)
ctrf = plotmap!(ax, finalage2000m, gridmetrics; colorrange, colormap, highclip, lowclip) # <- need to fix wrapping longitude for contour levels
myhidexdecorations!(ax, true)
𝑧 = rich("z", font = :italic)
Label(fig[1, 0]; text = rich(𝑧, " = 2000 m"), rotation = π / 2, tellheight = false)

# 4000m age
ax = Axis(fig[2, 1]; yticks, xticks, xtickformat, ytickformat)
ctrf = plotmap!(ax, finalage4000m, gridmetrics; colorrange, colormap, highclip, lowclip) # <- need to fix wrapping longitude for contour levels
Label(fig[2, 0]; text = rich(𝑧, " = 4000 m"), rotation = π / 2, tellheight = false)

# Colorbar
cb = Colorbar(fig[1:2, 2], ctrf; label = "age (years)", vertical = true, flipaxis = true, ticks = levels[1:2:end])
cb.height = Relative(2 / 3)

outputfile = joinpath(outputdir, "AA_age_maps_$(Ncycles)_postmaxagefix.png")
@info "Saving image file:\n  $(outputfile)"
save(outputfile, fig)


# Below I try to plot the drift but it does not work with the AA data,
# since its drift is large as it is accelerated. I should
# use the ACCESS-ESM1.5 output instead, but that data has been erased.
# I could potentially retsart a single 1-yr run from each AA age to plot it.
# Fraction of ocean volume with an ideal age drift of <10 years per 1,000 years
# Δage = diff(age4D, dims = 4)
x = cyclestart_idx
startidx, endidx = extrema(x) # Cannot use the last index of start of cycle (model has not finished corresponding run)
Δage = (age[Decade = startidx .. endidx] .- cyclestart_age[Decade = startidx .. endidx]).data
Δage[.!wet3D, :] .= NaN;
fracdrift = abs.(Δage) / 10
issmalldrift = fracdrift .< 0.01
volumeintegralissmalldrift = dropdims(nansum(issmalldrift .* v3D, dims = (1, 2, 3)), dims = (1, 2, 3))
fracsmalldrift = volumeintegralissmalldrift ./ nansum(v3D)

volumeintegralfracdrift = dropdims(nansum(fracdrift .* v3D, dims = (1, 2, 3)), dims = (1, 2, 3))
meanfracdrift = volumeintegralfracdrift ./ nansum(v3D)
volumeintegralfracdrift2000 = dropdims(nansum(fracdrift[:, :, i2000m, :] .* v3D[:, :, i2000m], dims = (1, 2)), dims = (1, 2))
meanfracdrift2000 = volumeintegralfracdrift2000 ./ nansum(v3D[:, :, i2000m])
volumeintegralfracdrift4000 = dropdims(nansum(fracdrift[:, :, i4000m, :] .* v3D[:, :, i4000m], dims = (1, 2)), dims = (1, 2))
meanfracdrift4000 = volumeintegralfracdrift4000 ./ nansum(v3D[:, :, i4000m])
maxfracdrift = dropdims(nanmaximum(Δage / 10, dims = (1, 2, 3)), dims = (1, 2, 3))
minfracdrift = dropdims(nanminimum(Δage / 10, dims = (1, 2, 3)), dims = (1, 2, 3))

fig = Figure(size = (800, 1000))
axisoptions = (
    xlabel = "ACCESS-ESM1.5 simulation years",
    # yscale = log10,
    ylabel = "small drift (< 0.01 yr/yr) frac. vol.",
    yminorticksvisible = true,
)
ax = Axis(fig[1, 1]; axisoptions...)
lines!(ax, 10 * x, fracsmalldrift)
xlims!(ax, 0, nothing)
ylims!(ax, 0, 1)
hidexdecorations!(ax, grid = false)

ax = Axis(fig[2, 1]; axisoptions..., ylabel = "drift (yr/yr)", yscale = log10)
lines!(ax, 10 * x, meanfracdrift2000, label = "mean 2000 m")
lines!(ax, 10 * x, meanfracdrift4000, label = "mean 4000 m")
lines!(ax, 10 * x, meanfracdrift, color = :black, label = "mean global")
xlims!(ax, 0, nothing)
ylims!(ax, 1.0e-5, 1)
hidexdecorations!(ax, grid = false)
axislegend(ax, position = :lt)

ax = Axis(fig[3, 1]; axisoptions..., ylabel = "drift (yr/yr)")
band!(ax, 10 * x, minfracdrift, maxfracdrift; color = :lightgray, label = "drift range global")
xlims!(ax, 0, nothing)
ylims!(ax, nothing, nothing)
axislegend(ax, position = :lt)

outputfile = joinpath(outputdir, "AA_fraction_of_small_drift_postmaxagefix.png")
@info "Saving image file:\n  $(outputfile)"
save(outputfile, fig)
