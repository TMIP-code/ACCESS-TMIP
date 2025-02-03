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

# test single file
# age_ds = open_dataset(joinpath(age_dir, "ocean_age.res_0000.nc"))

# Open all files into single xarray
Ndecades = length(age_files)
age_ds = open_mfdataset(DimArray(age_files, Dim{:Decade}(1:Ndecades)))
age = readcubedata(age_ds.age_global)
age = dropdims(age, dims = 4)
age4D = age.data
age4D[.!wet3D, :] .= NaN;

maxage = dropdims(nanmaximum(age4D, dims = (1,2,3)), dims = (1,2,3))
i2000m = argmin(abs.(zt .- 2000))
volumeintegralage2000m = dropdims(nansum((age4D .* v3D)[:,:,i2000m,:], dims = (1,2)), dims = (1,2))
volumeintergal2000m = dropdims(nansum(v3D[:,:,i2000m], dims = (1,2)), dims = (1,2))
meanage2000m = volumeintegralage2000m ./ volumeintergal2000m
i4000m = argmin(abs.(zt .- 4000))
volumeintegralage4000m = dropdims(nansum((age4D .* v3D)[:,:,i4000m,:], dims = (1,2)), dims = (1,2))
volumeintergal4000m = dropdims(nansum(v3D[:,:,i4000m], dims = (1,2)), dims = (1,2))
meanage4000m = volumeintegralage4000m ./ volumeintergal4000m

include("plotting_functions.jl")


fig = Figure()
axisoptions = (
    limits = (0, 10Ndecades, nothing, nothing),
    # xticks = 0:50:10Ndecades,
    xlabel = "ACCESS-ESM1.5 simulation years"
)
ax1 = Axis(fig[1, 1]; axisoptions..., ylabel = "age (years)")
x = range(start = 0, step = 10, length = Ndecades)
lines!(ax1, x, maxage; label = "maximum")
lines!(ax1, x, meanage2000m; label = "2000m mean")
lines!(ax1, x, meanage4000m; label = "4000m mean")
ylims!(ax1, 0, nothing)
axislegend(ax1; position = :lt, framevisible = false)
hidexdecorations!(ax1, grid = false)

ax2 = Axis(fig[2, 1]; axisoptions..., ylabel = "10-yr cycle Δage (years)")
x = range(start = 5, step = 10, length = Ndecades - 1)
lines!(ax2, x, diff(maxage); label = "maximum")
lines!(ax2, x, diff(meanage2000m); label = "2000m mean")
lines!(ax2, x, diff(meanage4000m); label = "4000m mean")

linkxaxes!(ax2, ax1)

outputfile = joinpath(outputdir, "AA_age_timeseries.png")
@info "Saving seqeff on sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)





# Rebuild dims from volcello
finalage = rebuild(age;
    data = age4D[:,:,:,end],
    dims = dims(volcello),
)
finalage2000m = finalage[lev = Near(2000)]
finalage4000m = finalage[lev = Near(4000)]

fig = Figure()
Label(fig[0, :], "Anderson-Acceleration age after $(10Ndecades) simulation years ($(Ndecades) cycles)", tellwidth = false)
yticks = -60:30:60
xticks = -120:60:120 + 360
levels = 0:100:1800
colormap = cgrad(:viridis, 19, categorical = true)
highclip = colormap[end]
colormap = cgrad(colormap[1:end-1], categorical = true)
colorrange = extrema(levels)


# 2000m age
ax = Axis(fig[1, 1]; yticks, xticks, xtickformat, ytickformat)
ctrf = plotmap!(ax, finalage2000m, gridmetrics; colorrange, colormap, highclip) # <- need to fix wrapping longitude for contour levels
myhidexdecorations!(ax, true)
𝑧 = rich("z", font = :italic)
Label(fig[1, 0]; text = rich(𝑧, " = 2000 m"), rotation = π/2, tellheight = false)

# 4000m age
ax = Axis(fig[2, 1]; yticks, xticks, xtickformat, ytickformat)
ctrf = plotmap!(ax, finalage4000m, gridmetrics; colorrange, colormap, highclip) # <- need to fix wrapping longitude for contour levels
Label(fig[2, 0]; text = rich(𝑧, " = 4000 m"), rotation = π/2, tellheight = false)

# Colorbar
cb = Colorbar(fig[1:2, 2], ctrf; label = "age (years)", vertical = true, flipaxis = true, ticks = 0:200:2000)
cb.height = Relative(2/3)

outputfile = joinpath(outputdir, "AA_age_maps_$(Ndecades).png")
@info "Saving seqeff on sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)




# # Below I try to plot the drift but it does not work with the AA data,
# # since its drift is large as it is accelerated. I should
# # use the ACCESS-ESM1.5 output instead, but that data has been erased.
# # I could potentially retsart a single 1-yr run from each AA age to plot it.
# # Fraction of ocean volume with an ideal age drift of <10 years per 1,000 years
# Δage = diff(age4D, dims = 4)
# issmalldrift = abs.(Δage) .< 1
# volumeintegralissmalldrift = dropdims(nansum(issmalldrift .* v3D, dims = (1,2,3)), dims = (1,2,3))
# fracsmalldrift = volumeintegralissmalldrift ./ nansum(v3D)

# fracdrift = abs.(Δage) / 10
# volumeintegralfracdrift = dropdims(nansum(fracdrift .* v3D, dims = (1,2,3)), dims = (1,2,3))
# meanfracdrift = volumeintegralfracdrift ./ nansum(v3D)

# fig = Figure()
# axisoptions = (
#     # xticks = 0:50:10Ndecades,
#     yticks = (exp10.(-3:0), ["0.001", "0.01", "0.1", "1"]),
#     xlabel = "ACCESS-ESM1.5 simulation years",
#     yscale = log10,
#     ylabel = "small drift frac. vol.",
#     yminorticksvisible = true,
#     )
# ax = Axis(fig[1, 1]; axisoptions...)
# x = range(start = 5, step = 10, length = Ndecades - 1)
# lines!(ax, x, fracsmalldrift)
# xlims!(ax, 0, 10Ndecades)
# ylims!(ax, 0.001, 1)

# ax = Axis(fig[2, 1]; axisoptions..., ylabel = "mean drift")
# lines!(ax, x, meanfracdrift)
# xlims!(ax, 0, 10Ndecades)
# ylims!(ax, 0.001, 1)

# outputfile = joinpath(outputdir, "AA_fraction_of_small_drift.png")
# @info "Saving seqeff on sea floor as image file:\n  $(outputfile)"
# save(outputfile, fig)