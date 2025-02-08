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
# using Unitful: s, yr, d, kyr
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
# all_members_inputdir2 = "/scratch/xv83/TMIP/data/$model/$experiment2/all_members/$(time_window2)/cyclomonth"
# Γoutyr3D_timemean_ds2 = open_dataset(joinpath(all_members_inputdir2, "adjointage_timemean.nc")).adjointage_timemean
# Γouts2 = map(srcnames) do srcname
#     src_P = sourcelocation(srcname)
#     src_i, src_j = Tuple(argmin(map(P -> norm(P .- src_P), zip(lon, lat))))
#     src_k = findlast(wet3D[src_i, src_j,:])
#     readcubedata(Γoutyr3D_timemean_ds2[src_i, src_j, src_k, :]).data
# end

# depths = [10000, 9000, 8000, 7000, 6000, 5000, 4000, 3000, 2000, 1000, 200, 0]
# # maxdepth = 6000
# maxdepth = 5000
# @time "simplepolygons" simplepolygons = [GeometryOps.simplify(VisvalingamWhyatt(tol=0.2), NaturalEarth.bathymetry(z).geometry) for z in depths[depths .≤ maxdepth]]
# # @time "simplepolygons" simplepolygons = [GeometryOps.simplify(NaturalEarth.bathymetry(z).geometry, tol=0.2) for z in depths[depths .≤ maxdepth]]
# # @time "simplepolygons" simplepolygons = [GeometryOps.simplify(NaturalEarth.bathymetry(z).geometry, ratio=0.05) for z in depths[depths .≤ maxdepth]]





# @info "Grab TTDs at injection location"
# members = ["r$(r)i1p1f1" for r in 1:3]
# # Gadi directory for input files
# inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/all_members/$(time_window)/cyclomonth"
# 𝒢_files = ["/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)/calgtilde.nc" for member in members]
# 𝒢_ds = open_mfdataset(DimArray(𝒢_files, Dim{:member}(members)))
# 𝒢 = readcubedata(𝒢_ds.calgtilde)
# 𝒢s = map(srcnames) do srcname
#     src_P = sourcelocation(srcname)
#     src_i, src_j = Tuple(argmin(map(P -> norm(P .- src_P), zip(lon, lat))))
#     𝒢[src_i, src_j, :, :]
# end
# # TTDs for 2090s circulation
# 𝒢_files2 = ["/scratch/xv83/TMIP/data/$model/$experiment2/$member/$(time_window2)/calgtilde.nc" for member in members]
# 𝒢_ds2 = open_mfdataset(DimArray(𝒢_files2, Dim{:member}(members)))
# 𝒢2 = readcubedata(𝒢_ds2.calgtilde)
# 𝒢s2 = map(srcnames) do srcname
#     src_P = sourcelocation(srcname)
#     src_i, src_j = Tuple(argmin(map(P -> norm(P .- src_P), zip(lon, lat))))
#     𝒢2[src_i, src_j, :, :]
# end
# TTD_time = 𝒢s[1].Ti |> Array

# @info "Grab sequestration efficiency at injection location"
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



fig = Figure(size=(800, 600))



#######
# MAP #
#######
limits = ((107, 158), (-48, -7))
xticks = -180:10:360
yticks = -90:10:90
xticks = (xticks, lonticklabel.(xticks))
yticks = (yticks, latticklabel.(yticks))
panela = Axis(fig[1, 1];
    backgroundcolor = :black,
    limits,
    # dest = "+proj=longlat +datum=WGS84",
    xgridvisible = true,
    ygridvisible = true,
    # yaxisposition = :right,
    xaxisposition = :top,
    xticks,
    yticks,
)
# panela2 = Axis(fig[1, 1]; xaxisposition = :top, limits, xticks, yticks)
# linkaxes!(panela2, panela)
# hidespines!(panela2)
# hidexdecorations!(panela2)
# hideydecorations!(panela)
# xlims!(panela2, limits[1])
# ylims!(panela2, limits[2])
# image!(panela2, -180..180, -90..90, GeoMakie.earth() |> rotr90; interpolate = false)
colors = cgrad(:jblue, length(depths), categorical=true)
# colors = cgrad(:blues, length(depths), categorical=true)
# colors = cgrad(:grays, length(depths), categorical=true, rev=true)

isinlimits(P) = !ismissing(P[1]) && !ismissing(P[2]) && (limits[1][1] ≤ P[1] ≤ limits[1][2]) && (limits[2][1] ≤ P[2] ≤ limits[2][2])
isinlimits(vecP::Vector) = any(isinlimits, vecP)
for (simplepolygon, color) in zip(reverse(simplepolygons), colors)
    p = poly!(panela, simplepolygon; color)
    translate!(p, 0, 0, -100)
end
# p = poly!(panela2, simplepolygons[2]; color = (:red, 0), strokecolor = :black, strokewidth = 1) # should be 5000m
# translate!(p, 0, 0, -90)
# poly!(panela2, GeoMakie.land(); color = :white, strokecolor = :black, strokewidth = 1)
poly!(panela, GeoMakie.land(); color = :lightgray, strokecolor = :black, strokewidth = 1)

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
    bbox = panela.scene.px_area,
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
# colors = cgrad(:Egypt, categorical=true)[[3]]
colors = cgrad(:galah, categorical=true)[[3]]
# color2 = cgrad(:Egypt, categorical=true)[4]
colors2 = cgrad(:galah, categorical=true)[[6]]
bdcolors = [(c + 2 * RGBA(1,1,1,1)) / 3 for c in colors]
bdcolors2 = [(c + 2 * RGBA(1,1,1,1)) / 3 for c in colors2]
# colors = Makie.wong_colors()[[1, 3, 6]]
# offsets = map(x -> x.* 2, [(-2, 1), (-2, -1), (2, -1)])
# offsets = map(x -> x.* 2, [(-2, 1)])
offsets = map(x -> x.* 2, [(6, -3)])
# aligns = [(:right, :bottom), (:right, :top), (:left, :top)]
# aligns = [(:right, :center), (:right, :center), (:left, :center)]
aligns = [(:right, :center)]
aligns = [(:left, :center)]
# texts = ["A"]
texts = ["injection\nlocation"]

for (ksrc, (srcname, offset, align, color, text)) in enumerate(zip(srcnames, offsets, aligns, colors, texts))
    src_P = sourcelocation(srcname)
    # sc1 = scatter!(panela2, src_P; marker=:circle, markersize=10, color=(:black, 0), strokecolor=:black, strokewidth=3)
    # sc2 = scatter!(panela, src_P; marker=:circle, markersize=10, color=(:black, 0), strokecolor=:black, strokewidth=4)
    # sc2 = scatter!(panela, src_P; marker=:circle, markersize=10, color=(:black, 0), strokecolor=color, strokewidth=2)
    # lines!(panela2, [src_P, src_P .+ offset]; color=:white)
    lines!(panela, kinkline(src_P .+ offset, src_P); color=:black)
    sc = scatter!(panela, src_P; marker=:star5, markersize=20, color=colors[ksrc], strokecolor=:black, strokewidth=1)
    # lines!(panela2, kinkline(src_P .+ offset, src_P); color=:black, linewidth=3)
    # lines!(panela2, kinkline(src_P .+ offset, src_P); color)
    text!(panela, src_P .+ offset; text, align, color=:black, strokecolor=:black)
    # text!(panela2, src_P .+ offset; text, align, color=:black, font=:bold, fontsize=18, strokecolor=:black, strokewidth=2)
    # text!(panela2, src_P .+ offset; text, align, color, font=:bold, fontsize=18)
    # translate!(sc1, 0, 0, 99)
    # translate!(sc2, 0, 0, 100)
end


# xlims!(ax, (100, 180))
# ylims!(ax, (-60, 0))
# cb = Colorbar(fog[1, 2], hm)

𝒓 = rich("r", font = :bold_italic)

ℰstr = rich("ℰ", rich("‾", offset = (-0.5, 0.15)))
ℰfun = rich(ℰstr, "(", 𝒓, ", τ)")
Γstr = rich("Γ", superscript("†"), rich("‾", offset = (-0.55, 0.25)), rich("‾", offset = (-0.85, 0.25)))
Γfun = rich(Γstr, rich("(", 𝒓, ")", offset = (0.4, 0)))
𝒢str = rich("𝒢", superscript("†"), rich("‾", offset = (-0.55, 0.25)), rich("‾", offset = (-0.85, 0.25)))
𝒢fun = rich(𝒢str, rich("(", 𝒓, ", τ)", offset = (0.4, 0)))







#######
# TTD #
#######
xmin, xmax = 0, 2000
ymin, ymax = 0, 1.2 # It's a pain but this needs to be set for my boxplots to match the size
axisoptions = (
    # ytrimspine = (false, true),
    # xtrimspine = (false, true),
    xgridvisible = true,
    ygridvisible = true,
    # rightspinevisible = false,
    # leftspinevisible = true,
    # topspinevisible = false,
    # bottomspinevisible = false,
    yaxisposition = :right,
    xaxisposition = :top,
    # xticks = -100:100:Nyears,
    xticks = WilkinsonTicks(7),
    # xticks = MultipleTicks(100),
    # yticks = 0:20:100,
    ylabel = rich("transit time distribution, ", 𝒢str, rich(" (kyr", superscript("−1"), ")", offset = (0.5, 0))),
    # ylabel = rich("adjoint propagator, ", 𝒢str, rich(" (kyr", superscript("−1"), ")", offset = (0.5, 0))),
    # xlabel = rich("time after injection, τ (years)"),
    xlabel = rich("time after injection (years)"),
    limits = (xmin, xmax, ymin, ymax),
)
panelb = Axis(fig[1, 2]; axisoptions...)
maxall𝒢s = ustrip.(kyr^-1, maximum(maximum(𝒢s[ksrc].data) for ksrc in eachindex(srcnames)) * s^-1)
for (ksrc, (srcname, text, color, bdcolor, color2, bdcolor2)) = enumerate(zip(srcnames, texts, colors, bdcolors, colors2, bdcolors2))
    vspan!(panelb, minimum(Γouts[ksrc]), maximum(Γouts[ksrc]), color = bdcolor)
    vspan!(panelb, minimum(Γouts2[ksrc]), maximum(Γouts2[ksrc]), color = bdcolor2)
    𝒢2030smin = ustrip.(kyr^-1, dropdims(minimum(𝒢s[ksrc].data, dims = 2), dims = 2) * s^-1)
    𝒢2030smax = ustrip.(kyr^-1, dropdims(maximum(𝒢s[ksrc].data, dims = 2), dims = 2) * s^-1)
    𝒢2090smin = ustrip.(kyr^-1, dropdims(minimum(𝒢s2[ksrc].data, dims = 2), dims = 2) * s^-1)
    𝒢2090smax = ustrip.(kyr^-1, dropdims(maximum(𝒢s2[ksrc].data, dims = 2), dims = 2) * s^-1)
    𝒢2030s = ustrip.(kyr^-1, dropdims(mean(𝒢s[ksrc].data, dims = 2), dims = 2) * s^-1)
    𝒢2090s = ustrip.(kyr^-1, dropdims(mean(𝒢s2[ksrc].data, dims = 2), dims = 2) * s^-1)
    bd2030s = band!(panelb, TTD_time, 𝒢2030smin, 𝒢2030smax; color = bdcolor)
    bd2090s = band!(panelb, TTD_time, 𝒢2090smin, 𝒢2090smax; color = bdcolor2)
    ln2030s = lines!(panelb, TTD_time, 𝒢2030s; color = color, linewidth=2, linecap=:round, joinstyle=:round)
    ln2090s = lines!(panelb, TTD_time, 𝒢2090s; color = color2, linewidth=2, linecap=:round, joinstyle=:round, linestyle = :dash)
    text!(panelb, 350, 𝒢2030smax[350]; text = "2030s", align = (:left, :bottom), offset = (3, 3), color = color)
    text!(panelb, 500, 𝒢2090smin[500]; text = "2090s", align = (:right, :top), offset = (-3, -3), color = color2)
end




############################
# sequestration efficiency #
############################
xmin, xmax = 0, 350
ymin, ymax = 70, 100
axisoptions = (
    # ytrimspine = (false, true),
    # xtrimspine = (false, true),
    xgridvisible = true,
    ygridvisible = true,
    # rightspinevisible = false,
    # leftspinevisible = true,
    # topspinevisible = false,
    # bottomspinevisible = false,
    yaxisposition = :left,
    xaxisposition = :bottom,
    # xticks = -100:100:Nyears,
    xticks = 0:100:xmax,
    # xticks = MultipleTicks(100),
    # yticks = 0:10:100,
    # yticks = [70, 90, 100],
    ylabel = rich("sequestration efficiency, ", ℰstr, " (%)"),
    # xlabel = rich("time after injection, τ (years)"),
    xlabel = rich("time after injection (years)"),
    limits = (xmin, xmax, ymin, ymax),
    backgroundcolor = (:white, 0),
)
panelc = Axis(fig[2, 1]; axisoptions...)
# dashed lines slicing ℰ
color = :gray
linestyle = :dash
align = (:center, :center)

# hlines!(panelc, [50, 90]; color, linestyle)
# vlines!(panelc, [100, 300, 1000]; color, linestyle)
# alternating bg bands instead of grid
# for ylev in 10:20:90
#     hspan!(panelc, ylev, ylev+10; color = (:black, 0.025))
# end
# Band for injection time window
# ibnd = vspan!(panelc, -10, 0; color = (:black, 0.1))
# text!(panelc, 10, 50; text = "10-year injection", rotation = π/2, align = (:center, :center))
for (ksrc, (srcname, text, color, bdcolor, color2, bdcolor2)) = enumerate(zip(srcnames, texts, colors, bdcolors, colors2, bdcolors2))
    ℰ2030smin = dropdims(minimum(100 * ℰs[ksrc].data, dims = 2), dims = 2)
    ℰ2030smax = dropdims(maximum(100 * ℰs[ksrc].data, dims = 2), dims = 2)
    ℰ2090smin = dropdims(minimum(100 * ℰs2[ksrc].data, dims = 2), dims = 2)
    ℰ2090smax = dropdims(maximum(100 * ℰs2[ksrc].data, dims = 2), dims = 2)
    ℰ2030s = dropdims(mean(100 * ℰs[ksrc].data, dims = 2), dims = 2)
    ℰ2090s = dropdims(mean(100 * ℰs2[ksrc].data, dims = 2), dims = 2)
    bd2030s = band!(panelc, TTD_time, ℰ2030smin, ℰ2030smax; color = bdcolor)
    bd2090s = band!(panelc, TTD_time, ℰ2090smin, ℰ2090smax; color = bdcolor2)
    ln2030s = lines!(panelc, TTD_time, ℰ2030s; color, linewidth=2, linecap=:round, joinstyle=:round)
    ln2090s = lines!(panelc, TTD_time, ℰ2090s; color = color2, linewidth=2, linecap=:round, joinstyle=:round, linestyle = :dash)
end
ylims!(panelc, ymax, ymin)





# Repeat for unzoomed time series
xmin, xmax = 0, 2000
ymin, ymax = 0, 100
axisoptions = (
    # ytrimspine = (false, true),
    # xtrimspine = (false, true),
    xgridvisible = true,
    ygridvisible = true,
    # rightspinevisible = false,
    # leftspinevisible = true,
    # topspinevisible = false,
    # bottomspinevisible = false,
    yaxisposition = :right,
    xaxisposition = :bottom,
    # xticks = -100:100:Nyears,
    xticks = WilkinsonTicks(7),
    # xticks = [0, 100, 300, 500, 1000, 1500, 2000],
    # xticks = MultipleTicks(100),
    # yticks = (0:10:100, [(mod(x, 20) == 0 ? string(x) : "") for x in 0:10:100]),
    # yticks = [0, 50, 90, 100],
    ylabel = rich("sequestration efficiency, ", ℰstr, " (%)"),
    # xlabel = rich("time after injection, τ (years)"),
    xlabel = rich("time after injection (years)"),
    limits = (xmin, xmax, ymin, ymax),
    backgroundcolor = (:white, 0),
)
paneld = Axis(fig[2, 2]; axisoptions...)
# dashed lines slicing ℰ
color = :gray
linestyle = :dash
align = (:center, :center)
# hlines!(paneld, [50, 90]; color, linestyle)
# vlines!(paneld, [100, 300, 1000]; color, linestyle)
for (ksrc, (srcname, text, color, bdcolor, color2, bdcolor2)) = enumerate(zip(srcnames, texts, colors, bdcolors, colors2, bdcolors2))
    vspan!(paneld, minimum(Γouts[ksrc]), maximum(Γouts[ksrc]), color = bdcolor)
    vspan!(paneld, minimum(Γouts2[ksrc]), maximum(Γouts2[ksrc]), color = bdcolor2)
    ℰ2030smin = dropdims(minimum(100 * ℰs[ksrc].data, dims = 2), dims = 2)
    ℰ2030smax = dropdims(maximum(100 * ℰs[ksrc].data, dims = 2), dims = 2)
    ℰ2090smin = dropdims(minimum(100 * ℰs2[ksrc].data, dims = 2), dims = 2)
    ℰ2090smax = dropdims(maximum(100 * ℰs2[ksrc].data, dims = 2), dims = 2)
    ℰ2030s = dropdims(mean(100 * ℰs[ksrc].data, dims = 2), dims = 2)
    ℰ2090s = dropdims(mean(100 * ℰs2[ksrc].data, dims = 2), dims = 2)
    bd2030s = band!(paneld, TTD_time, ℰ2030smin, ℰ2030smax; color = bdcolor)
    bd2090s = band!(paneld, TTD_time, ℰ2090smin, ℰ2090smax; color = bdcolor2)
    ln2030s = lines!(paneld, TTD_time, ℰ2030s; color = color, linewidth=2, linecap=:round, joinstyle=:round)
    ln2090s = lines!(paneld, TTD_time, ℰ2090s; color = color2, linewidth=2, linecap=:round, joinstyle=:round, linestyle = :dash)
end
ylims!(paneld, ymax, ymin)




# Add zoom lines
rectattrs = (strokecolor = :lightgray, linestyle = :solid)
lineattrs = (color = :lightgray, linestyle = :solid)
zoom_lines!(panelc, paneld; rectattrs, lineattrs)


colsize!(fig.layout, 1, Relative(0.4))
rowgap!(fig.layout, 20)
colgap!(fig.layout, 20)



#############
# Box plots #
#############
# (must come after resizing to get correct width of bar plots)
# Box plot of mean age inside
valuesΓ2030s = reduce(vcat, Γouts)
categories = reduce(vcat, fill(label, 40) for label in texts)
# Place Γ box plot close to top and stagger its display down for each injection locations (max 4)
categorypositions = maxall𝒢s * (1 .- reduce(vcat, fill(ilabel, 40) for ilabel in reverse(eachindex(texts))) / 4)
color = reduce(vcat, fill(color, 40) for color in colors)
# I want a 10 pixel width = width_dataspace / limits_dataspace * limits_figspace
# so width_dataspace = 10 * limits_dataspace / limits_figspace
limits_dataspace = panelb.finallimits[].widths[2]
limits_figspace = fullproject(panelb, panelb.finallimits[]).widths[2]
width = 20 * limits_dataspace / limits_figspace
boxplot!(panelb, categorypositions, valuesΓ2030s; color=color, orientation = :horizontal, width, strokewidth = 1, whiskerwidth = :match)
# Repeat for 2090s
valuesΓ2090s = reduce(vcat, Γouts2)
categories = reduce(vcat, fill(label, 40) for label in texts)
# Place Γ box plot close to top and stagger its display down for each injection locations (max 4)
categorypositions = maxall𝒢s * (1 .- reduce(vcat, fill(ilabel, 40) for ilabel in reverse(eachindex(texts))) / 4)
color = reduce(vcat, fill(color, 40) for color in colors)
color2 = reduce(vcat, fill(color, 40) for color in colors2)
# I want a 10 pixel width = width_dataspace / limits_dataspace * limits_figspace
# so width_dataspace = 10 * limits_dataspace / limits_figspace
limits_dataspace = panelb.finallimits[].widths[2]
limits_figspace = fullproject(panelb, panelb.finallimits[]).widths[2]
width = 20 * limits_dataspace / limits_figspace
boxplot!(panelb, categorypositions, valuesΓ2090s; color=color2, orientation = :horizontal, width, strokewidth = 1, whiskerwidth = :match)
text!(panelb, (minimum(valuesΓ2030s) + maximum(valuesΓ2090s)) / 2, categorypositions[1]; text = rich("mean time, ", Γstr), align = (:center, :bottom), offset = (0, 30))
text!(panelb, (minimum(valuesΓ2030s) + maximum(valuesΓ2030s)) / 2, categorypositions[1]; text = "2030s", align = (:center, :bottom), offset = (0, 15), fontsize = 10)
text!(panelb, (minimum(valuesΓ2090s) + maximum(valuesΓ2090s)) / 2, categorypositions[1]; text = "2090s", align = (:center, :bottom), offset = (0, 15), fontsize = 10)
# text!(panelb, maximum(values), categorypositions[1]; text = rich("mean time ", Γstr), align = (:left, :center), offset = (5, 0))
# given ℰ
for (ℰsdecade, colors) in zip((ℰs, ℰs2), (colors, colors2))
    for panel in (panelc, paneld)
        for ℰval in [50, 90]
            values = [findfirst(100 * ℰi .< ℰval) for ℰ in ℰsdecade for ℰi in eachcol(ℰ.data)]
            categories = reduce(vcat, fill(label, 3) for label in texts)
            categorypositions = fill(ℰval, length(categories))
            color = reduce(vcat, fill(color, 3) for color in colors)
            limits_dataspace = panel.finallimits[].widths[2]
            limits_figspace = fullproject(panel, panel.finallimits[]).widths[2]
            width = 20 * limits_dataspace / limits_figspace
            boxplot!(panel, categorypositions, values; color=color, orientation = :horizontal, width, strokewidth = 1, whiskerwidth = :match)
            (ℰval == 50) && (panel == paneld) && (ℰsdecade == ℰs) && text!(panel, minimum(values), categorypositions[1]; text = "median\ntime", align = (:right, :center), offset = (-10, 0))
            (ℰval == 90) && (panel == panelc) && (ℰsdecade == ℰs) && text!(panel, minimum(values), categorypositions[1]; text = "10th percentile\ntime", align = (:right, :center), offset = (-10, 0))
        end
        # given τ
        for τval = [100, 300, 1000]
            values = [100 * ℰi[τval] for ℰ in ℰsdecade for ℰi in eachcol(ℰ.data)]
            categories = reduce(vcat, fill(label, 3) for label in texts)
            categorypositions = fill(τval, length(categories))
            color = reduce(vcat, fill(color, 3) for color in colors)
            limits_dataspace = panel.finallimits[].widths[1]
            limits_figspace = fullproject(panel, panel.finallimits[]).widths[1]
            width = 20 * limits_dataspace / limits_figspace
            boxplot!(panel, categorypositions, values; color=color, orientation = :vertical, width, strokewidth = 1, whiskerwidth = :match)
        end
    end
end
# band for mean age





labeloptions = (
    font = :bold,
    align = (:left, :top),
    offset = (5, -2),
    space = :relative,
    fontsize = 24
)
for (ax, label) in zip([panela, panelb, panelc, paneld], string.('a':'h'))
    text!(ax, 0, 1; text = label, labeloptions..., strokecolor = :white, strokewidth = 3)
    text!(ax, 0, 1; text = label, labeloptions...)
end


# save plot
outputdir = joinpath(fixedvarsinputdir, "all_members")
mkpath(outputdir)
outputfile = joinpath(outputdir, "TTD_seqeff_timeseries_1loc.png")
@info "Saving injection location as image file:\n  $(outputfile)"
save(outputfile, fig)
outputfile = joinpath(outputdir, "TTD_seqeff_timeseries_1loc.pdf")
@info "Saving injection location as image file:\n  $(outputfile)"
save(outputfile, fig)
