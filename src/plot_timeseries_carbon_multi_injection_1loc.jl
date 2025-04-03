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
# using GeometryBasics
# using LibGEOS

# include("plotting_functions.jl")





# # Load matrix and grid metrics
# # @show model = "ACCESS-ESM1-5"
# # @show experiment = ARGS[1]
# # @show member = ARGS[2]
# # @show time_window = ARGS[3]
# @show model = "ACCESS-ESM1-5"
# # @show experiment = "historical"
# @show experiment = "ssp370"
# @show experiment2 = "ssp370"
# # @show time_window = "Jan1850-Dec1859"
# # @show time_window = "Jan1990-Dec1999"
# @show time_window = "Jan2030-Dec2039"
# @show time_window2 = "Jan2090-Dec2099"

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
# srcnames = ["Karratha"]
# year_start = parse(Int, time_window[4:7])
# # Nyears = 501
# # Nyears = 2001
# Nyears = 1001
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
#         @show size(umass)
#         yaxdata[member=At(member), source=At(srcname)] .= 100 * [0.0; umass] / src_mass
#         # TODO: Check that src_P match (since I have changed these along the way)
#         src_P = load(outputfile, "src_P")
#         src_P == sourcelocation(srcname) || @info "issue $member $srcname $scr_P ≠ $(sourcelocation(srcname))"
#     end
# end

# @info "Grab mean reemergence time at injection locations for all members"
# all_members_inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/all_members/$(time_window)/cyclomonth"
# Γoutyr3D_timemean_ds = open_dataset(joinpath(all_members_inputdir, "adjointage_timemean.nc")).adjointage_timemean
# Γouts = map(srcnames) do srcname
#     src_P = sourcelocation(srcname)
#     src_i, src_j = Tuple(argmin(map(P -> norm(P .- src_P), zip(lon, lat))))
#     src_k = findlast(wet3D[src_i, src_j,:])
#     readcubedata(Γoutyr3D_timemean_ds[src_i, src_j, src_k, :]).data
# end

# depths = [10000, 9000, 8000, 7000, 6000, 5000, 4000, 3000, 2000, 1000, 200, 0]
# # maxdepth = 6000
# maxdepth = 5000
# @time "simplepolygons" simplepolygons = [GeometryOps.simplify(VisvalingamWhyatt(tol=0.2), NaturalEarth.bathymetry(z).geometry) for z in depths[depths .≤ maxdepth]]
# # @time "simplepolygons" simplepolygons = [GeometryOps.simplify(NaturalEarth.bathymetry(z).geometry, tol=0.2) for z in depths[depths .≤ maxdepth]]
# # @time "simplepolygons" simplepolygons = [GeometryOps.simplify(NaturalEarth.bathymetry(z).geometry, ratio=0.05) for z in depths[depths .≤ maxdepth]]





# @info "Grab TTDs"
# members = ["r$(r)i1p1f1" for r in 1:3]
# # Gadi directory for input files
# inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/all_members/$(time_window)/cyclomonth"
# ℰ_files = ["/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)/calE.nc" for member in members]
# ℰ_ds = open_mfdataset(DimArray(ℰ_files, Dim{:member}(members)))
# ℰ = readcubedata(ℰ_ds.calE)
# ℰs = map(srcnames) do srcname
#     src_P = sourcelocation(srcname)
#     src_i, src_j = Tuple(argmin(map(P -> norm(P .- src_P), zip(lon, lat))))
#     ℰ[src_i, src_j, :, :]
# end
# # TTDs for 2090s circulation
# ℰ_files2 = ["/scratch/xv83/TMIP/data/$model/$experiment2/$member/$(time_window2)/calE.nc" for member in members]
# ℰ_ds2 = open_mfdataset(DimArray(ℰ_files2, Dim{:member}(members)))
# ℰ2 = readcubedata(ℰ_ds2.calE)
# ℰs2 = map(srcnames) do srcname
#     src_P = sourcelocation(srcname)
#     src_i, src_j = Tuple(argmin(map(P -> norm(P .- src_P), zip(lon, lat))))
#     ℰ2[src_i, src_j, :, :]
# end
# TTD_time = ℰs[1].Ti |> Array




fig = Figure(size=(800, 600))


limits = ((107, 158), (-48, -7))
xticks = -180:10:360
yticks = -90:10:90
xticks = (xticks, lonticklabel.(xticks))
yticks = (yticks, latticklabel.(yticks))
panelh = Axis(fig[2, 2];
    backgroundcolor = :black,
    limits,
    # dest = "+proj=longlat +datum=WGS84",
    xgridvisible = true,
    ygridvisible = true,
    # yaxisposition = :right,
    xticks,
    yticks,
)
panelh2 = Axis(fig[2, 2]; yaxisposition = :right, limits, xticks, yticks)
linkaxes!(panelh2, panelh)
hidespines!(panelh2)
hidexdecorations!(panelh2)
hideydecorations!(panelh)
# xlims!(panelh2, limits[1])
# ylims!(panelh2, limits[2])
# image!(panelh2, -180..180, -90..90, GeoMakie.earth() |> rotr90; interpolate = false)
colors = cgrad(:jblue, length(depths), categorical=true)
# colors = cgrad(:blues, length(depths), categorical=true)
# colors = cgrad(:grays, length(depths), categorical=true, rev=true)

isinlimits(P) = !ismissing(P[1]) && !ismissing(P[2]) && (limits[1][1] ≤ P[1] ≤ limits[1][2]) && (limits[2][1] ≤ P[2] ≤ limits[2][2])
isinlimits(vecP::Vector) = any(isinlimits, vecP)
for (simplepolygon, color) in zip(reverse(simplepolygons), colors)
    p = poly!(panelh2, simplepolygon; color)
    translate!(p, 0, 0, -100)
end
# p = poly!(panelh2, simplepolygons[2]; color = (:red, 0), strokecolor = :black, strokewidth = 1) # should be 5000m
# translate!(p, 0, 0, -90)
# poly!(panelh2, GeoMakie.land(); color = :white, strokecolor = :black, strokewidth = 1)
poly!(panelh2, GeoMakie.land(); color = :lightgray, strokecolor = :black, strokewidth = 1)

# Add depth colorbar inside
cbarlabelformat(x) = isinteger(x) ? string(round(Int, x)) : string(x)
ilow = only(findall(depths .== 0))
ihigh = only(findall(depths .== maxdepth))
Colorbar(fig;
    limits = (1,length(simplepolygons) - 2),
    # ticks = (1:length(simplepolygons) - 2, cbarlabelformat.(1e-3 * reverse(depths[ihigh + 1 : ilow - 1]))),
    ticks = (1:length(simplepolygons) - 2, string.(reverse(depths[ihigh + 1 : ilow - 1]))),
    colormap = cgrad(reverse(colors)[ihigh + 2 : ilow - 1], categorical = true, rev = true),
    highclip = reverse(colors)[ihigh + 1],
    lowclip = reverse(colors)[ilow],
    label = "seafloor depth (m)",
    height = 10,
    width = 140,
    labelsize = 10,
    vertical = false,
    flipaxis = true,
    ticklabelsize = 10,
    bbox = panelh2.scene.viewport,
    alignmode = Outside(10),
    valign = :bottom,
    halign = :left,
    # topspinecolor = :white,
    # bottomspinecolor = :white,
    # rightspinecolor = :white,
    # leftspinecolor = :white,
    # ticklabelcolor = :white,
    # labelcolor = :white,
    # tickcolor = :white
)


# depth2D = nansum(gridmetrics.thkcello; dim = 3)
# depth2D[.!wet3D[:,:,1]] .= NaN
# # plotmap!(ax, depth2D, gridmetrics; colormap = :deep, colorscale = log10)
# hm = plotmap!(ax, depth2D, gridmetrics; colormap = :GnBu, colorscale = log10)
# TODO: read src_P somehow (file name or variable inside file)
# colors = cgrad(:Egypt, categorical=true)[[3, 4, 1]]
colors = cgrad(:Egypt, categorical=true)[[3]]
# colors = Makie.wong_colors()[[1, 3, 6]]
# offsets = map(x -> x.* 2, [(-2, 1), (-2, -1), (2, -1)])
offsets = map(x -> x.* 2, [(-2, 1)])
# aligns = [(:right, :bottom), (:right, :top), (:left, :top)]
# aligns = [(:right, :center), (:right, :center), (:left, :center)]
aligns = [(:right, :center)]
texts = ["A"]

for (ksrc, (srcname, offset, align, color, text)) in enumerate(zip(srcnames, offsets, aligns, colors, texts))
    src_P = sourcelocation(srcname)
    # sc = scatter!(panelh2, src_P; marker=:star5, markersize=20, color=colors[ksrc], strokecolor=:black, strokewidth=1)
    # sc1 = scatter!(panelh2, src_P; marker=:circle, markersize=10, color=(:black, 0), strokecolor=:black, strokewidth=3)
    sc2 = scatter!(panelh2, src_P; marker=:circle, markersize=10, color=(:black, 0), strokecolor=:black, strokewidth=4)
    sc2 = scatter!(panelh2, src_P; marker=:circle, markersize=10, color=(:black, 0), strokecolor=color, strokewidth=2)
    # lines!(panelh2, [src_P, src_P .+ offset]; color=:white)
    lines!(panelh2, kinkline(src_P .+ offset, src_P); color=:black, linewidth = 1)
    # lines!(panelh2, kinkline(src_P .+ offset, src_P); color=:black, linewidth=3)
    # lines!(panelh2, kinkline(src_P .+ offset, src_P); color)
    text!(panelh2, src_P .+ offset; text, align, color=:black, strokecolor=:black)
    # text!(panelh2, src_P .+ offset; text, align, color=:black, font=:bold, fontsize=18, strokecolor=:black, strokewidth=2)
    # text!(panelh2, src_P .+ offset; text, align, color, font=:bold, fontsize=18)
    # translate!(sc1, 0, 0, 99)
    translate!(sc2, 0, 0, 100)
end


# xlims!(ax, (100, 180))
# ylims!(ax, (-60, 0))
# cb = Colorbar(fog[1, 2], hm)

ℰstr = rich("ℰ", rich("‾", offset = (-0.5, 0.15)))
# τstr = rich("τ", font = :italic)
𝒓 = rich("r", font = :bold_italic)
ℰfun = rich(ℰstr, "(", 𝒓, ", τ)")

axisoptions = (
    # ytrimspine = (false, true),
    # xtrimspine = (false, true),
    xgridvisible = true,
    ygridvisible = false,
    # rightspinevisible = false,
    # leftspinevisible = true,
    # topspinevisible = false,
    # bottomspinevisible = false,
    yaxisposition = :left,
    xaxisposition = :top,
    # xticks = -100:100:Nyears,
    xticks = WilkinsonTicks(7),
    # xticks = MultipleTicks(100),
    yticks = 0:20:100,
    ylabel = rich("sequestration efficiency, ", ℰfun, " (%)"),
    xlabel = rich("time after injection, τ (years)"),
)
# xmin, xmax = -30, Nyears - 1
# xmin, xmax = 0, Nyears - 1
xmin, xmax = 0, 2000
# xmin, xmax = 0, 1020
ymin, ymax = 0, 100
# xmin, xmax = year_start, 2100
panela = Axis(fig[1, 1]; axisoptions...)
x = -10:(length(times) - 11)
# BG of discrete gradient?
# for ylev in 0:10:100
#     hspan!(panela, 0, ylev; color = (:black, 0.025))
# end
# for ylev in 20:40:80
#     hspan!(panela, ylev, ylev+20; color = (:black, 0.025))
# end







color = :gray
linestyle = :dash
align = (:center, :center)

lines!(panela, [-1000, 1850 - 80, NaN, 1850 + 80, 3000], fill(90, 5); color, linestyle)
lines!(panela, [-1000, 1850 - 80, NaN, 1850 + 80, 3000], fill(50, 5); color, linestyle)
lines!(panela, fill(100, 5), [-10, 5 - 3, NaN, 5 + 3, 110]; color, linestyle)
lines!(panela, fill(300, 5), [-10, 5 - 3, NaN, 5 + 3, 110]; color, linestyle)
lines!(panela, fill(1000, 5), [-10, 5 - 3, NaN, 5 + 3, 110]; color, linestyle)

text!(panela, 1850, 90; text = "g", align, color)
text!(panela, 1850, 50; text = "f", align, color)
text!(panela, 100, 5; text = "b", align, color)
text!(panela, 300, 5; text = "c", align, color)
text!(panela, 1000, 5; text = "d", align, color)

for ylev in 10:20:90
    hspan!(panela, ylev, ylev+10; color = (:black, 0.025))
end
# Band for injection time window
# ibnd = vspan!(panela, -10, 0; color = (:black, 0.1))
# text!(panela, 10, 50; text = "10-year injection", rotation = π/2, align = (:center, :center))
for (ksrc, (srcname, text)) = enumerate(zip(srcnames, texts))
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
    # bd = band!(panela, x, Cseqmin, Cseqmax; color=(:black, 0.3))
    bd = band!(panela, x, clamp.(Cseqmin, 0, 100), clamp.(Cseqmax, 0, 100); color=(color, 0.3))
    # for Cseqksrc_m in eachslice(Cseqksrc, dims = 2)
    #     # cannot work if saved umass have different time span
    #     # y = 100 * [0; data[m]["umass"][:, ksrc]] / data[m]["src_mass"][ksrc]
    #     # TODO use a ribbon instead of plotting each trajectory
    #     lines!(panela, x, Cseqksrc_m; color = :black, linewidth=0.1)
    # end
    # ln = lines!(panela, x, Cseqmean; color=:black, linewidth=2, linecap=:round, joinstyle=:round)
    ln = lines!(panela, x, clamp.(Cseqmean, 0, 100); color, linewidth=2, linecap=:round, joinstyle=:round)
    # The line below is from the adjoint propagator
    ln2 = lines!(panela, TTD_time, dropdims(mean(100 * ℰs[ksrc].data, dims = 2), dims = 2); color = :black, linewidth=2, linecap=:round, joinstyle=:round)
    ln3 = lines!(panela, TTD_time, dropdims(mean(100 * ℰs2[ksrc].data, dims = 2), dims = 2); color = :black, linewidth=2, linecap=:round, joinstyle=:round, linestyle = :dash)
    i = 250 - 10ksrc
    # (ksrc == 1) && (text = "tracer injected at $text")
    # text!(panela, x[i], Cseqmean[i]; text, offset = (1.5, 1.5), align = (:left, :bottom), color=:black)
    text!(panela, x[i], Cseqmean[i]; text, offset = (-1.5, 1.5), align = (:right, :bottom), color=:black)
    # text!(panela, x[i], Cseqmean[i]; text, offset = (3, 3), align = (:left, :bottom), color, fontsize=18)
    # Add a bar chart in the middle?
    # for xbar in 2100:100:2500
    #     ix = only(findall(x .== xbar))
    #     categories = fill(xbar, Nmembers)
    #     values = Cseqksrc[time = ix].data
    #     boxplot!(panela, categories, values; color, markersize = 0.1, width = 20)
    # end
    # Add vertical bands of mean reemergence times across members
    # vspan!(panela, extrema(Γouts[ksrc])...; color = (color, 0.3))
    # vspan!(panela, quantile(Γouts[ksrc], [0.25, 0.75])...; color = (color, 0.3))
    # vlines!(panela, quantile(Γouts[ksrc], [0.5]); color, linestyle = :dash)
    # vlines!(panela, Γouts[ksrc][1]; color = (:black, 0.3))
    # Add barplot of mean reemergence times across members
    # categories = fill(5ksrc, 40)
    # values = Γouts[ksrc]
    # boxplot!(panela, categories, values; color, markersize = 0.1, width = 5, orientation = :horizontal)
    xlims!(panela, (xmin, xmax))
    # ylims!(panela, (0, 102))
    # ylims!(panela, (0, nothing))
    ylims!(panela, (ymin, ymax))
    # ylims!(panela, (70, 102))
end
# hidexdecorations!(panela, grid=false)


# Bpx plot below

# yticks = (reverse(1:3), fill("", 3))
yticks = (reverse(eachindex(texts)), texts)
Γup = rich("Γ", superscript("↑"))
axisoptions = (
    # ytrimspine = (false, true),
    # xtrimspine = (false, true),
    xgridvisible = true,
    ygridvisible = false,
    yticksvisible = true,
    # rightspinevisible = false,
    # leftspinevisible = false,
    # topspinevisible = false,
    yaxisposition = :left,
    xaxisposition = :bottom,
    # xticks = -100:100:Nyears,
    xticks = WilkinsonTicks(7),
    # xticks = MultipleTicks(100),
    yticks,
    # ylabel = "injection location",
    xlabel = rich("time after injection, τ (years)"),
)


panelsefg = fig[2, 1] = GridLayout()
ylabel = rich("""
    mean time
    """,
    Γup, "(", 𝒓, ")"
)
panele = Axis(panelsefg[1, 1]; axisoptions..., ylabel)
values = reduce(vcat, Γouts)
categories = reduce(vcat, fill(label, 40) for label in texts)
categorypositions = reduce(vcat, fill(ilabel, 40) for ilabel in reverse(eachindex(texts)))
color = reduce(vcat, fill(color, 40) for color in colors)
boxplot!(panele, categorypositions, values; color=color, orientation = :horizontal)
# for (ksrc, text) in enumerate(texts)
#     # (ksrc == 1) && (text = "mean reemergence time at $text")
#     text!(panele, minimum(Γouts[ksrc]), Nsrc + 1 - ksrc; text, offset = (-3,0), align = (:right, :center))
# end
# raincloudoptions = (
#     boxplot_nudge = -0.5,
#     boxplot_width = 0.5,
#     jitter_width = 0.2,
#     clouds = nothing,
#     # hist_bins = 0:5:3000,
#     # cloud_width = 0.3,
#     gap = 0.01,
# )
# rainclouds!(panele, categorypositions, values; raincloudoptions..., color=color, orientation = :horizontal)
hidexdecorations!(panele, grid = false)





ylabel = rich("""
    median time
    """,
    ℰstr, " = 50 %"
)
panelf = Axis(panelsefg[2, 1]; axisoptions..., ylabel)
values = [findfirst(ℰi .< 0.5) for ℰ in ℰs for ℰi in eachcol(ℰ.data)]
categories = reduce(vcat, fill(label, 3) for label in texts)
categorypositions = reduce(vcat, fill(ilabel, 3) for ilabel in reverse(eachindex(texts)))
color = reduce(vcat, fill(color, 3) for color in colors)
boxplot!(panelf, categorypositions, values; color=color, orientation = :horizontal)
linkxaxes!(panelf, panele)
xlims!(panelf, (xmin, xmax))
hidexdecorations!(panelf, grid = false)

ylabel = rich("""
    time to 10 % leak
    """,
    ℰstr, " = 90 %"
)
panelg = Axis(panelsefg[3, 1]; axisoptions..., ylabel)
values = [findfirst(ℰi .< 0.9) for ℰ in ℰs for ℰi in eachcol(ℰ.data)]
categories = reduce(vcat, fill(label, 3) for label in texts)
categorypositions = reduce(vcat, fill(ilabel, 3) for ilabel in reverse(eachindex(texts)))
color = reduce(vcat, fill(color, 3) for color in colors)
boxplot!(panelg, categorypositions, values; color=color, orientation = :horizontal)
linkxaxes!(panelg, panelf)
xlims!(panelg, (xmin, xmax))


panelsbcd = fig[1, 2] = GridLayout()

axisoptions = (
    # ytrimspine = (false, true),
    # xtrimspine = (false, true),
    xgridvisible = false,
    ygridvisible = true,
    xticksvisible = true,
    # rightspinevisible = false,
    # leftspinevisible = false,
    # topspinevisible = false,
    yaxisposition = :right,
    xaxisposition = :top,
    # xticks = -100:100:Nyears,
    yticks = WilkinsonTicks(7),
    # xticks = MultipleTicks(100),
    xticks = (reverse(eachindex(texts)), texts),
    # ylabel = "injection location",
    ylabel = rich("sequestration efficiency, ", ℰfun, " (%)"),
)

panelb = Axis(panelsbcd[1, 1]; axisoptions..., xlabel = "τ = 100 years")
values = [100 * ℰi[100 + 1] for ℰ in ℰs for ℰi in eachcol(ℰ.data)]
categories = reduce(vcat, fill(label, 3) for label in texts)
categorypositions = reduce(vcat, fill(ilabel, 3) for ilabel in reverse(eachindex(texts)))
color = reduce(vcat, fill(color, 3) for color in colors)
boxplot!(panelb, categorypositions, values; color=color, orientation = :vertical)
linkyaxes!(panelb, panela)
ylims!(panelb, (ymin, ymax))
hideydecorations!(panelb)


panelc = Axis(panelsbcd[1, 2]; axisoptions..., xlabel = "τ = 300 years")
values = [100 * ℰi[300 + 1] for ℰ in ℰs for ℰi in eachcol(ℰ.data)]
categories = reduce(vcat, fill(label, 3) for label in texts)
categorypositions = reduce(vcat, fill(ilabel, 3) for ilabel in reverse(eachindex(texts)))
color = reduce(vcat, fill(color, 3) for color in colors)
boxplot!(panelc, categorypositions, values; color=color, orientation = :vertical)
linkyaxes!(panelc, panela)
ylims!(panelc, (ymin, ymax))
hideydecorations!(panelc)


paneld = Axis(panelsbcd[1, 3]; axisoptions..., xlabel = "τ = 1000 years")
values = [100 * ℰi[1000 + 1] for ℰ in ℰs for ℰi in eachcol(ℰ.data)]
categories = reduce(vcat, fill(label, 3) for label in texts)
categorypositions = reduce(vcat, fill(ilabel, 3) for ilabel in reverse(eachindex(texts)))
color = reduce(vcat, fill(color, 3) for color in colors)
boxplot!(paneld, categorypositions, values; color=color, orientation = :vertical)
linkyaxes!(paneld, panela)
ylims!(paneld, (ymin, ymax))
hideydecorations!(paneld, label = false, ticklabels = false, ticks = false)


for panel in (panelb, panelc, paneld)
    for ylev in 10:20:90
        bg = hspan!(panel, ylev, ylev+10; color = (:black, 0.025))
        translate!(bg, 0, 0, -100)
    end
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

rowsize!(fig.layout, 1, Relative(2/3))
colsize!(fig.layout, 1, Relative(2/3))
rowgap!(panelsefg, 10)
colgap!(panelsefg, 10)
rowgap!(panelsbcd, 10)
colgap!(panelsbcd, 10)

# colgap!(fig.layout, 1, -50)

labeloptions = (
    font = :bold,
    align = (:left, :top),
    offset = (5, -2),
    space = :relative,
    fontsize = 24
)
for (ax, label) in zip([panela, panelb, panelc, paneld, panele, panelf, panelg, panelh], string.('a':'h'))
    text!(ax, 0, 1; text = label, labeloptions..., strokecolor = :white, strokewidth = 3)
    text!(ax, 0, 1; text = label, labeloptions...)
end


# save plot
outputdir = joinpath(fixedvarsinputdir, "all_members")
mkpath(outputdir)
outputfile = joinpath(outputdir, "injected_tracer_timeseries_v2_1loc.png")
@info "Saving injection location as image file:\n  $(outputfile)"
save(outputfile, fig)
outputfile = joinpath(outputdir, "injected_tracer_timeseries_v2_1loc.pdf")
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
#     ylabel = "fraction sequestered (%)",
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