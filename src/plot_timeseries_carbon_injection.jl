# # qsub -I -P xv83 -q hugemem -l mem=360GB -l storage=gdata/gh0+scratch/gh0+scratch/xv83 -l walltime=02:00:00 -l ncpus=48

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
# @show experiment = "historical"
# @show member = "r1i1p1f1"
# @show time_window = "Jan1990-Dec1999"

# lumpby = "month"
# steps = 1:12
# Nsteps = length(steps)
# δt = ustrip(s, 1yr / Nsteps) # TODO maybe use exact mean number of days (more important for monthly because Feb)?


# # Gadi directory for input files
# inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)"
# cycloinputdir = joinpath(inputdir, "cyclo$lumpby")
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


# outputfile = joinpath(inputdir, "injected_tracer_timeseries.jld2")
# @info "Saving injection time series file:\n  $(outputfile)"
# data = load(outputfile)


fig = Figure(size=(500, 300))
ax = Axis(fig[1, 1];
    title="Seafloor injected carbon sequestration",
    ytrimspine = true,
    xtrimspine = true,
    xgridvisible = false,
    ygridvisible = false,
    rightspinevisible = false,
    topspinevisible = false,
    xticks = 2000:100:10_000,
    yticks = 0:20:100,
    ylabel = "Fraction sequestered (%)",
    xlabel = rich("year (repeating 2020s climatology)"),
)
x = 2020 .+ (0:length(data["umass"]))
y = 100 * [0; data["umass"]] / data["src_mass"]
ylevs = 0:10:100

color = Makie.wong_colors()[6]
xmin, xmax = 2000, 2500
for ylev in ylevs
    ix = findlast(y .> ylev)
    (isnothing(ix) || x[ix] > xmax) && continue
    hspan!(ax, 0, ylev; color = (:black, 0.05))
    lines!(ax, [x[ix], x[ix], NaN, x[ix], x[ix]], [0, 10, NaN, 30, ylev]; color = (:black, 0.2), linestyle = :dash)
    text!(ax, x[ix], 20; text = "$ylev%", color = (:black, 0.2), rotation = π/2, align = (:center, :center))
end
ibnd = vspan!(ax, 2020, 2030; color = (color, 0.1))
text!(ax, 2010, 50; text = "injection during 2020s", rotation = π/2, align = (:center, :center), color)
ln = lines!(ax, x, y; color, linewidth=2)
xlims!(ax, (2000, 2500))
ylims!(ax, (0, 102))

# Add insert of injection location
# The inset axis
inset_ax = Axis(fig[1, 1],
    width=Relative(0.4),
    height=Relative(0.4),
    halign=0.15,
    valign=0.4,
    backgroundcolor=:lightgray)

hidedecorations!(inset_ax)
depth2D = nansum(gridmetrics.thkcello; dim = 3)
depth2D[.!wet3D[:,:,1]] .= NaN
plotmap!(inset_ax, depth2D, gridmetrics; colormap = :dense)
# src_P = (110, -15) # <- Choose (lon,lat) of source here
# TODO: read src_P somehow (file name or variable inside file)
sc = scatter!(inset_ax, src_P; marker=:star5, markersize=15, color, strokecolor=:black, strokewidth=1)
translate!(sc, 0, 0, 100)

# save plot
outputfile = joinpath(inputdir, "injected_tracer_timeseries.png")
@info "Saving injection location as image file:\n  $(outputfile)"
save(outputfile, fig)