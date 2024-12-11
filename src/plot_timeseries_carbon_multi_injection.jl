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
# using NaturalEarth
# using GeometryOps

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

# function sourcelocation(srcname)
#     if srcname == "Karratha"
#         return (115.45849390000001,-16.56466979999999) # Carnarvon Basin?" North West of Australia
#     elseif srcname == "Portland"
#         return (141.73529860000008,-38.93477809999996) # Otway Basin" South West of Melbourne (West of Tas)
#     elseif srcname == "Marlo"
#         return (149.05333500000006, -38.25798499999996) # "Shark 1" Gippsland Basin" South East (East of Tas)
#     else
#         error("No source name matchin $srcname")
#     end
# end

# members = map(m -> "r$(m)i1p1f1", 1:40)
# srcnames = ["Karratha", "Portland", "Marlo"]
# year_start = parse(Int, time_window[4:7])
# Nyears = 501
# times = range(Date("$year_start-01-01"); step=Year(1), length=Nyears)
# axlist = (
#     Dim{:time}(times),
#     Dim{:source}(srcnames),
#     Dim{:member}(members),
# )
# Nmembers = length(members)
# Nsrc = length(srcnames)
# data = fill(NaN, (Nyears, Nsrc, Nmembers))
# yaxdata = YAXArray(axlist, data)

# for member in members
#     for srcname in srcnames
#         datainputdir = joinpath(fixedvarsinputdir, experiment, member, time_window)
#         outputfile = joinpath(datainputdir, "$(srcname)_timeseries_injected.jld2")
#         isfile(outputfile) || continue
#         @info "Loading injection time series file:\n  $(outputfile)"
#         umass = load(outputfile, "umass")
#         src_mass = load(outputfile, "src_mass")
#         yaxdata[member=At(member), source=At(srcname)] .= 100 * [0.0; umass] / src_mass
#         # TODO: Check that src_P match (since I have changed these along the way)
#         src_P = load(outputfile, "src_P")
#         src_P == sourcelocation(srcname) || @info "issue $member $srcname $scr_P ≠ $(sourcelocation(srcname))"
#     end
# end

depths = [10000, 9000, 8000, 7000, 6000, 5000, 4000, 3000, 2000, 1000, 200, 0]
maxdepth = 6000
@time "simplepolygons" simplepolygons = [GeometryOps.simplify(VisvalingamWhyatt(tol=0.2), NaturalEarth.bathymetry(z).geometry) for z in depths[depths .≤ maxdepth]]
# @time "simplepolygons" simplepolygons = [GeometryOps.simplify(NaturalEarth.bathymetry(z).geometry, tol=0.2) for z in depths[depths .≤ maxdepth]]
# @time "simplepolygons" simplepolygons = [GeometryOps.simplify(NaturalEarth.bathymetry(z).geometry, ratio=0.05) for z in depths[depths .≤ maxdepth]]

fig = Figure(size=(900, 400))

# # Add location as first axis
# ax = Axis(fig[1, 1];
#     yticks = (-60:20:0, ["$(-i)°S" for i in -60:20:0]),
#     xticks = (100:20:180, ["$(i)°E" for i in 100:20:180]),
# )
# Try GeoMakie
limits = ((107, 158), (-48, -7))
xticks = -180:10:360
yticks = -90:10:90
xticks = (xticks, lonticklabel.(xticks))
yticks = (yticks, latticklabel.(yticks))
# ga = GeoAxis(fig[1, 1];
ga = Axis(fig[1, 1];
    backgroundcolor = :black,
    limits,
    # dest = "+proj=longlat +datum=WGS84",
    xgridvisible = true,
    ygridvisible = true,
    # yaxisposition = :right,
    xticks,
    yticks,
)
# xlims!(ga, limits[1])
# ylims!(ga, limits[2])
# image!(ga, -180..180, -90..90, GeoMakie.earth() |> rotr90; interpolate = false)
colors = cgrad(:jblue, length(depths), categorical=true)
# colors = cgrad(:blues, length(depths), categorical=true)
# colors = cgrad(:grays, length(depths), categorical=true, rev=true)

isinlimits(P) = !ismissing(P[1]) && !ismissing(P[2]) && (limits[1][1] ≤ P[1] ≤ limits[1][2]) && (limits[2][1] ≤ P[2] ≤ limits[2][2])
isinlimits(vecP::Vector) = any(isinlimits, vecP)
for (simplepolygon, color) in zip(reverse(simplepolygons), colors)
    poly!(ga, simplepolygon; color)
end
# poly!(ga, GeoMakie.land(); color = :white, strokecolor = :black, strokewidth = 1)
poly!(ga, GeoMakie.land(); color = :lightyellow, strokecolor = :black, strokewidth = 1)
# depth2D = nansum(gridmetrics.thkcello; dim = 3)
# depth2D[.!wet3D[:,:,1]] .= NaN
# # plotmap!(ax, depth2D, gridmetrics; colormap = :deep, colorscale = log10)
# hm = plotmap!(ax, depth2D, gridmetrics; colormap = :GnBu, colorscale = log10)
# TODO: read src_P somehow (file name or variable inside file)
colors = cgrad(:Egypt, categorical=true)[[3, 4, 1]]
# colors = Makie.wong_colors()[[1, 3, 6]]
offsets = map(x -> x.* 2, [(-2, 1), (-2, -1), (2, -1)])
# aligns = [(:right, :bottom), (:right, :top), (:left, :top)]
aligns = [(:right, :center), (:right, :center), (:left, :center)]
texts = ["A", "B", "C"]
function kinkline(A, B, c=2/3)
    # A -- C
    #       \
    #        B
    return [(0.9A[1] + 0.1B[1], A[2]), (c * B[1] + (1-c) * A[1], A[2]), B]
end
for (ksrc, (srcname, offset, align, color, text)) in enumerate(zip(srcnames, offsets, aligns, colors, texts))
    src_P = sourcelocation(srcname)
    # sc = scatter!(ga, src_P; marker=:star5, markersize=20, color=colors[ksrc], strokecolor=:black, strokewidth=1)
    # sc1 = scatter!(ga, src_P; marker=:circle, markersize=10, color=(:black, 0), strokecolor=:black, strokewidth=3)
    sc2 = scatter!(ga, src_P; marker=:circle, markersize=10, color=(:black, 0), strokecolor=:black, strokewidth=4)
    sc2 = scatter!(ga, src_P; marker=:circle, markersize=10, color=(:black, 0), strokecolor=color, strokewidth=2)
    # lines!(ga, [src_P, src_P .+ offset]; color=:white)
    lines!(ga, kinkline(src_P .+ offset, src_P); color=:black)
    # lines!(ga, kinkline(src_P .+ offset, src_P); color=:black, linewidth=3)
    # lines!(ga, kinkline(src_P .+ offset, src_P); color)
    text!(ga, src_P .+ offset; text, align, color=:black, strokecolor=:black)
    # text!(ga, src_P .+ offset; text, align, color=:black, font=:bold, fontsize=18, strokecolor=:black, strokewidth=2)
    # text!(ga, src_P .+ offset; text, align, color, font=:bold, fontsize=18)
    # translate!(sc1, 0, 0, 99)
    translate!(sc2, 0, 0, 100)
end

# xlims!(ax, (100, 180))
# ylims!(ax, (-60, 0))
# cb = Colorbar(fog[1, 2], hm)

axisoptions = (
    ytrimspine = (false, true),
    xtrimspine = (false, true),
    xgridvisible = true,
    ygridvisible = false,
    rightspinevisible = false,
    # leftspinevisible = true,
    topspinevisible = false,
    # yaxisposition = :right,
    xticks = -100:100:500,
    yticks = 0:20:100,
    ylabel = "Fraction sequestered (%)",
    xlabel = rich("year after injection (repeating $(year_start)s)"),
)
xmin, xmax = -30, 500
# xmin, xmax = year_start, 2100
ax = Axis(fig[1, 2]; axisoptions...)
x = -10:(length(times) - 11)
# BG of discrete gradient?
# for ylev in 0:10:100
#     hspan!(ax, 0, ylev; color = (:black, 0.025))
# end
for ylev in 20:40:80
    hspan!(ax, ylev, ylev+20; color = (:black, 0.025))
end
ibnd = vspan!(ax, -10, 0; color = (:black, 0.1))
text!(ax, 10, 50; text = "10-year injection", rotation = π/2, align = (:center, :center))
for (ksrc, srcname) = enumerate(srcnames)
    Cseqksrc = yaxdata[source = At(srcname)]
    Cseqmean = dropdims(mean(Cseqksrc, dims=:member), dims=:member).data
    Cseqmin = dropdims(minimum(Cseqksrc, dims=:member), dims=:member).data
    Cseqmax = dropdims(maximum(Cseqksrc, dims=:member), dims=:member).data
    color = colors[ksrc]
    # Cseqmin = dropdims(minimum(Cseqksrc, dims = 2), dims = 2)
    # Cseqmax = dropdims(maximum(Cseqksrc, dims = 2), dims = 2)
    # TODO use saved duration of injection in output file
    # hide injection and pre injection
    inan = x .≤ 0
    Cseqmean[inan] .= NaN
    Cseqmin[inan] .= NaN
    Cseqmax[inan] .= NaN
    # bd = band!(ax, x, Cseqmin, Cseqmax; color=(:black, 0.3))
    bd = band!(ax, x, Cseqmin, Cseqmax; color=(color, 0.3))
    # for Cseqksrc_m in eachslice(Cseqksrc, dims = 2)
    #     # cannot work if saved umass have different time span
    #     # y = 100 * [0; data[m]["umass"][:, ksrc]] / data[m]["src_mass"][ksrc]
    #     # TODO use a ribbon instead of plotting each trajectory
    #     lines!(ax, x, Cseqksrc_m; color = :black, linewidth=0.1)
    # end
    # ln = lines!(ax, x, Cseqmean; color=:black, linewidth=2, linecap=:round, joinstyle=:round)
    ln = lines!(ax, x, Cseqmean; color, linewidth=2, linecap=:round, joinstyle=:round)
    i = 250 - 10ksrc
    text = texts[ksrc] * (ksrc == 1 ? " injection location" : "")
    text!(ax, x[i], Cseqmean[i]; text, offset = (10, 1), align = (:left, :bottom), color=:black)
    # text!(ax, x[i], Cseqmean[i]; text, offset = (10, 1), align = (:left, :bottom), color, fontsize=18)
    # Add a bar chart in the middle?
    # for xbar in 2100:100:2500
    #     ix = only(findall(x .== xbar))
    #     categories = fill(xbar, Nmembers)
    #     values = Cseqksrc[time = ix].data
    #     boxplot!(ax, categories, values; color, markersize = 0.1, width = 20)
    # end
    xlims!(ax, (xmin, xmax))
    # ylims!(ax, (0, 102))
    ylims!(ax, (0, nothing))
    # ylims!(ax, (70, 102))
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

# colgap!(fig.layout, 1, -50)

for (ax, label) in zip([ga, ax], ["a", "b"])
    text!(
        ax, 0, 0,
        text = label,
        font = :bold,
        align = (:left, :bottom),
        offset = (5, 2),
        space = :relative,
        fontsize = 24
    )
end


# save plot
outputdir = joinpath(fixedvarsinputdir, "all_members")
mkpath(outputdir)
# outputfile = joinpath(outputdir, "injected_tracer_timeseries.png")
outputfile = joinpath(outputdir, "injected_tracer_timeseries.pdf")
@info "Saving injection location as image file:\n  $(outputfile)"
save(outputfile, fig)


# fig = Figure(size=(500, 300Nsrc))
# axisoptions = (
#     ytrimspine = true,
#     xtrimspine = false,
#     xgridvisible = true,
#     ygridvisible = false,
#     rightspinevisible = false,
#     topspinevisible = false,
#     xticks = 100 * floor(Int, year_start / 100):100:100 * floor(Int, (year_start + Nyears) / 100),
#     yticks = 0:20:100,
#     ylabel = "Fraction sequestered (%)",
#     xlabel = rich("year (repeating $(time_window) climatology)"),
# )
# ylevs = 0:10:100
# colors = Makie.wong_colors()[[1, 3, 6]]
# xmin, xmax = extrema(axisoptions.xticks)
# for (ksrc, srcname) = enumerate(srcnames)
#     Cseqksrc = yaxdata[source = At(srcname)].data
#     Cseqmean = dropdims(nanmean(Cseqksrc, dims=2), dims=2)
#     Cseqmin = dropdims(nanminimum(Cseqksrc, dims=2), dims=2)
#     Cseqmax = dropdims(nanmaximum(Cseqksrc, dims=2), dims=2)
#     color = colors[ksrc]
#     # Cseqmin = dropdims(minimum(Cseqksrc, dims = 2), dims = 2)
#     # Cseqmax = dropdims(maximum(Cseqksrc, dims = 2), dims = 2)
#     ax = Axis(fig[ksrc, 1]; axisoptions...)
#     x = year.(times)
#     for ylev in ylevs
#         ix = findlast(Cseqmean .> ylev)
#         (isnothing(ix) || (ix == length(Cseqmean)) || x[ix] > xmax) && continue
#         hspan!(ax, 0, ylev; color = (:black, 0.025))
#         # lines!(ax, [x[ix], x[ix], NaN, x[ix], x[ix]], [0, 10, NaN, 30, ylev]; color = (:black, 0.2), linestyle = :dash)
#         # text!(ax, x[ix], 20; text = "$ylev%", color = (:black, 0.2), rotation = π/2, align = (:center, :center))
#     end
#     # TODO use saved duration of injection in output file
#     ibnd = vspan!(ax, year_start, year_start + 10; color = (:black, 0.1))
#     if ksrc == 1
#         text!(ax, year_start - 10, 50; text = "injection during $time_window", rotation = π/2, align = (:center, :center), color)
#     end
#     bd = band!(ax, x, Cseqmin, Cseqmax; color=(color, 0.3))
#     # for Cseqksrc_m in eachslice(Cseqksrc, dims = 2)
#     #     # cannot work if saved umass have different time span
#     #     # y = 100 * [0; data[m]["umass"][:, ksrc]] / data[m]["src_mass"][ksrc]
#     #     # TODO use a ribbon instead of plotting each trajectory
#     #     lines!(ax, x, Cseqksrc_m; color = :gray, linewidth=1)
#     # end
#     ln = lines!(ax, x, Cseqmean; color, linewidth=2)
#     xlims!(ax, (xmin, xmax))
#     ylims!(ax, (0, 120))
#     myhidexdecorations!(ax, ksrc < Nsrc)
# end

# rowgap!(fig.layout, 1, 0.0)
# rowgap!(fig.layout, 2, 0.0)


# # # Add insert of injection location
# # # The inset axis
# # inset_ax = Axis(fig[1, 1],
# #     width=Relative(0.4),
# #     height=Relative(0.4),
# #     halign=0.15,
# #     valign=0.4,
# #     backgroundcolor=:lightgray)

# # hidedecorations!(inset_ax)
# # depth2D = nansum(gridmetrics.thkcello; dim = 3)
# # depth2D[.!wet3D[:,:,1]] .= NaN
# # plotmap!(inset_ax, depth2D, gridmetrics; colormap = :dense)
# # # src_P = (110, -15) # <- Choose (lon,lat) of source here
# # # TODO: read src_P somehow (file name or variable inside file)
# # sc = scatter!(inset_ax, src_P; marker=:star5, markersize=15, color, strokecolor=:black, strokewidth=1)
# # translate!(sc, 0, 0, 100)

# # save plot
# outputdir = joinpath(fixedvarsinputdir, "all_members")
# mkpath(outputdir)
# outputfile = joinpath(outputdir, "injected_tracer_timeseries.png")
# @info "Saving injection location as image file:\n  $(outputfile)"
# save(outputfile, fig)