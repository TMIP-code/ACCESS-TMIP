# # qsub -I -P xv83 -l mem=32GB -l storage=scratch/gh0+scratch/xv83 -l walltime=02:00:00 -l ncpus=6

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
# # @show experiment = "historical"
# @show experiment = "ssp370"
# # @show time_window = "Jan1850-Dec1859"
# # @show time_window = "Jan1990-Dec1999"
# @show time_window = "Jan2030-Dec2039"

# lumpby = "month"
# steps = 1:12
# Nsteps = length(steps)
# δt = ustrip(s, 1yr / Nsteps) # TODO maybe use exact mean number of days (more important for monthly because Feb)?


# # Gadi directory for input files
# fixedvarsinputdir = "/scratch/xv83/TMIP/data/$model"
# # Load areacello and volcello for grid geometry
# volcello_ds = open_dataset(joinpath(fixedvarsinputdir, "volcello.nc"))
# areacello_ds = open_dataset(joinpath(fixedvarsinputdir, "areacello.nc"))

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

validmembers = keys(data)

src_names = keys(first(values(data))["src_Ps"])

# ymean = mean([data[m]["umass"] for m in keys(data)], dim) # cannot work if saved umass have different time span

Cseq = NamedTuple(map(enumerate(src_names)) do (ksrc, src_name)
    src_name => reduce(hcat, 100 * [0; data[m]["umass"]][:, ksrc] ./ data[m]["src_mass"][ksrc] for m in keys(data))
end)
# ymean = mean(100 * [0; data[m]["umass"]] ./ data[m]["src_mass"] for m in keys(data))
# ymax = maximum(100 * [0; data[m]["umass"]] ./ data[m]["src_mass"] for m in keys(data))
# ymin = minimum(100 * [0; data[m]["umass"]] ./ data[m]["src_mass"] for m in keys(data))

Nyearsplus1, Nvalidmembers = size(first(Cseq))
Nsrc = length(src_names)

year_start = parse(Int, time_window[4:7])
fig = Figure(size=(500, 300Nsrc))
axisoptions = (
    ytrimspine = true,
    xtrimspine = false,
    xgridvisible = false,
    ygridvisible = false,
    rightspinevisible = false,
    topspinevisible = false,
    xticks = 1800:100:3000,
    yticks = 0:20:100,
    ylabel = "Fraction sequestered (%)",
    xlabel = rich("year (repeating $(time_window) climatology)"),
)
ylevs = 0:10:100
color = Makie.wong_colors()[6]
xmin, xmax = year_start - 20, year_start + 500
for (ksrc, src_name) = enumerate(src_names)
    Cseqksrc = Cseq[ksrc]
    Cseqmean = dropdims(mean(Cseqksrc, dims = 2), dims = 2)
    Cseqmin = dropdims(minimum(Cseqksrc, dims = 2), dims = 2)
    Cseqmax = dropdims(maximum(Cseqksrc, dims = 2), dims = 2)
    ax = Axis(fig[ksrc, 1]; axisoptions...)
    x = year_start .+ (0:Nyearsplus1 - 1)
    for ylev in ylevs
        ix = findlast(Cseqmean .> ylev)
        (isnothing(ix) || (ix == length(Cseqmean)) || x[ix] > xmax) && continue
        hspan!(ax, 0, ylev; color = (:black, 0.05))
        lines!(ax, [x[ix], x[ix], NaN, x[ix], x[ix]], [0, 10, NaN, 30, ylev]; color = (:black, 0.2), linestyle = :dash)
        text!(ax, x[ix], 20; text = "$ylev%", color = (:black, 0.2), rotation = π/2, align = (:center, :center))
    end
    # TODO use saved duration of injection in output file
    ibnd = vspan!(ax, year_start, year_start + 10; color = (color, 0.1))
    text!(ax, year_start - 10, 50; text = "injection during $time_window", rotation = π/2, align = (:center, :center), color)
    bd = band!(ax, x, Cseqmin, Cseqmax; color=(color, 0.3))
    for Cseqksrc_m in eachslice(Cseqksrc, dims = 2)
        # cannot work if saved umass have different time span
        # y = 100 * [0; data[m]["umass"][:, ksrc]] / data[m]["src_mass"][ksrc]
        # TODO use a ribbon instead of plotting each trajectory
        lines!(ax, x, Cseqksrc_m; color = :gray, linewidth=1)
    end
    ln = lines!(ax, x, Cseqmean; color, linewidth=2)
    xlims!(ax, (xmin, xmax))
    ylims!(ax, (0, 102))
    myhidexdecorations!(ax, ksrc < Nsrc)
end




# # Add insert of injection location
# # The inset axis
# inset_ax = Axis(fig[1, 1],
#     width=Relative(0.4),
#     height=Relative(0.4),
#     halign=0.15,
#     valign=0.4,
#     backgroundcolor=:lightgray)

# hidedecorations!(inset_ax)
# depth2D = nansum(gridmetrics.thkcello; dim = 3)
# depth2D[.!wet3D[:,:,1]] .= NaN
# plotmap!(inset_ax, depth2D, gridmetrics; colormap = :dense)
# # src_P = (110, -15) # <- Choose (lon,lat) of source here
# # TODO: read src_P somehow (file name or variable inside file)
# sc = scatter!(inset_ax, src_P; marker=:star5, markersize=15, color, strokecolor=:black, strokewidth=1)
# translate!(sc, 0, 0, 100)

# save plot
outputdir = joinpath(fixedvarsinputdir, "all_members")
mkpath(outputdir)
outputfile = joinpath(outputdir, "injected_tracer_timeseries.png")
@info "Saving injection location as image file:\n  $(outputfile)"
save(outputfile, fig)