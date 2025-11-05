# # qsub -I -P xv83 -q hugemem -l mem=360GB -l storage=scratch/gh0+scratch/xv83 -l walltime=02:00:00 -l ncpus=48

# using Pkg
# Pkg.activate(".")
# Pkg.instantiate()


# ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
# using OceanTransportMatrixBuilder
# using NetCDF
# using YAXArrays
# using DataFrames
# using DimensionalData
# using SparseArrays
# using LinearAlgebra
# using Unitful
# using Unitful: s, yr, d
# using Statistics
# using Format
# using Dates
# using FileIO
# using LinearSolve
# import Pardiso # import Pardiso instead of using (to avoid name clash?)
# using NonlinearSolve
# using ProgressMeter
# try
#     using CairoMakie
# catch
#     using CairoMakie
# end
# using GeoMakie
# using OceanBasins
# using NaNStatistics

# include("plotting_functions.jl")


# # Load matrix and grid metrics
# # @show model = "ACCESS-ESM1-5"
# # @show experiment = ARGS[1]
# # @show member = ARGS[2]
# # @show time_window = ARGS[3]
# @show model = "ACCESS-ESM1-5"
# @show experiment = "ssp370"
# @show member = "r1i1p1f1"
# @show time_window = "Jan2030-Dec2039"


# # Gadi directory for input files
# inputdir = "/scratch/xv83/TMIP/data/$model"
# # Load areacello and volcello for grid geometry
# volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
# areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))

# # Load fixed variables in memory
# areacello = readcubedata(areacello_ds.areacello)
# volcello = readcubedata(volcello_ds.volcello)
# lon = readcubedata(volcello_ds.lon)
# lat = readcubedata(volcello_ds.lat)
# lev = volcello_ds.lev
# # Identify the vertices keys (vary across CMIPs / models)
# volcello_keys = propertynames(volcello_ds)
# lon_vertices_key = volcello_keys[findfirst(x -> occursin("lon", x) & occursin("vert", x), string.(volcello_keys))]
# lat_vertices_key = volcello_keys[findfirst(x -> occursin("lat", x) & occursin("vert", x), string.(volcello_keys))]
# lon_vertices = readcubedata(getproperty(volcello_ds, lon_vertices_key))
# lat_vertices = readcubedata(getproperty(volcello_ds, lat_vertices_key))

# # Make makegridmetrics
# gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
# (; lon_vertices, lat_vertices, v3D, ) = gridmetrics

# # Make indices
# indices = makeindices(v3D)
# (; N, wet3D) = indices


# members = map(m -> "r$(m)i1p1f1", 1:40)

# data = Dict{String, Any}()
# for member in members
#     @show datainputdir = joinpath(fixedvarsinputdir, experiment, member, time_window)
#     outputfile = joinpath(datainputdir, "timeseries_injected.jld2")
#     isfile(outputfile) || continue
#     @info "Loading injection time series file:\n  $(outputfile)"
#     data[member] = load(outputfile)
# end


colors = Makie.wong_colors()[[1, 3, 6]]
fig = Figure(size = (400, 300))
ax = Axis(
    fig[1, 1];
    yticks = (-60:20:0, ["$(-i)°S" for i in -60:20:0]),
    xticks = (100:20:180, ["$(i)°E" for i in 100:20:180]),
)

depth2D = nansum(gridmetrics.thkcello; dim = 3)
depth2D[.!wet3D[:, :, 1]] .= NaN
# plotmap!(ax, depth2D, gridmetrics; colormap = :deep, colorscale = log10)
hm = plotmap!(ax, depth2D, gridmetrics; colormap = :GnBu, colorscale = log10)
src_P = (110, -15) # <- Choose (lon,lat) of source here
# TODO: read src_P somehow (file name or variable inside file)
for (ksrc, (src_name, src_P)) in enumerate(pairs(first(values(data))["src_Ps"]))
    sc = scatter!(ax, src_P; marker = :star5, markersize = 10, color = colors[ksrc], strokecolor = :black, strokewidth = 1)
    translate!(sc, 0, 0, 100)
end
xlims!(ax, (100, 180))
ylims!(ax, (-60, 0))
cb = Colorbar(fog[1, 2], hm)

# save plot
outputfile = joinpath(inputdir, "injected_tracer_location.png")
@info "Saving injection location as image file:\n  $(outputfile)"
save(outputfile, fig)
