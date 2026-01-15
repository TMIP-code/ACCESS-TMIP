# qsub -I -P y99 -q normal -l mem=47GB -l storage=scratch/xv83+gdata/fs38 -l walltime=02:00:00 -l ncpus=12
# This is like Fig. 3 in Pasquier et al. (GRL, 2025) but for a different location

using Pkg
Pkg.activate(".")
Pkg.instantiate()

ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
using OceanTransportMatrixBuilder
using NetCDF
using YAXArrays
using DataFrames
using DimensionalData
using SparseArrays
using LinearAlgebra
using Unitful
using Unitful: s, yr, d, kyr
using Statistics
using Format
using Dates
using FileIO
using ProgressMeter
try
    using CairoMakie
catch
    using CairoMakie
end
using GeoMakie
using OceanBasins
using NaNStatistics
using NaturalEarth
using GeometryOps
using GeometryBasics
using Format
using LibGEOS
using Shapefile
using StatsBase
using KernelDensity

include("plotting_functions.jl")


# Load matrix and grid metrics
# @show model = "ACCESS-ESM1-5"
# @show experiment = ARGS[1]
# @show member = ARGS[2]
# @show time_window = ARGS[3]
@show model = "ACCESS-ESM1-5"
# @show experiment = "historical"
@show experiment = "ssp370"
@show time_window = "Jan2030-Dec2039"

# Load \Gamma out
κVdeep = 3.0e-5 # m^2/s
κVML = 1.0      # m^2/s
κH = 300.0      # m^2/s
κVdeep_str = "kVdeep" * format(κVdeep, conversion = "e")
κVML_str = "kVML" * format(κVML, conversion = "e")
κH_str = "kH" * format(κH, conversion = "d")
upwind = false
upwind_str = upwind ? "" : "_centered"
upwind_str2 = upwind ? "upwind" : "centered"
yearly = false
yearly_str = yearly ? "_yearly" : ""
yearly_str2 = yearly ? "(yearly)" : ""

# Gadi directory for input files
fixedvarsinputdir = "/scratch/xv83/TMIP/data/$model"
# Load areacello and volcello for grid geometry
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
(; lon_vertices, lat_vertices, v3D) = gridmetrics

# Make indices
indices = makeindices(v3D)
(; N, wet3D) = indices

# model depths
deptho = readcubedata(open_dataset("/g/data/fs38/publications/CMIP6/ScenarioMIP/CSIRO/ACCESS-ESM1-5/ssp370/r1i1p1f1/Ofx/deptho/gn/files/d20191115/deptho_Ofx_ACCESS-ESM1-5_ssp370_r1i1p1f1_gn.nc").deptho)

# Depth Polygons for map of location
depths = [10000, 9000, 8000, 7000, 6000, 5000, 4000, 3000, 2000, 1000, 200, 0]
# maxdepth = 6000
maxdepth = 5000
@time "simplepolygons" simplepolygons = [GeometryOps.simplify(VisvalingamWhyatt(tol = 0.2), NaturalEarth.bathymetry(z).geometry) for z in depths[depths .≤ maxdepth]]
# @time "simplepolygons" simplepolygons = [GeometryOps.simplify(NaturalEarth.bathymetry(z).geometry, tol=0.2) for z in depths[depths .≤ maxdepth]]
# @time "simplepolygons" simplepolygons = [GeometryOps.simplify(NaturalEarth.bathymetry(z).geometry, ratio=0.05) for z in depths[depths .≤ maxdepth]]

@info "Cloud locations"
cloudsourceshapefile = "/scratch/xv83/TMIP/data/ACCESS-ESM1-5/final_positions/opendrift_nearbed_velocity_20yrs.shp"
shp_table = Shapefile.Table(cloudsourceshapefile)
cloudsourcelocations = shp_table.geometry
# Find wet grid cells containing source locations and count number of cloud points in each grid cell
# Trick to mask out land points (make the lat and lon = -1000)
lon2 = copy(lon.data); lon2[.!wet3D[:, :, 1]] .= -1000
lat2 = copy(lat.data); lat2[.!wet3D[:, :, 1]] .= -1000
@time cloudsourceijs = [Tuple(argmin(map(P -> norm(P .- (src_P.x, src_P.y)), zip(lon2, lat2)))) for src_P in cloudsourcelocations]
cloudsourcedict = countmap(cloudsourceijs)
# Print locations, depths, and counts
clouddf = DataFrame(
    lon = Float64[],
    lat = Float64[],
    depth = Float64[],
    deptho = Float64[],
    fraction = Float64[],
    src_i = Int[],
    src_j = Int[],
    src_k = Int[],
)
for (src_ij, count) in pairs(cloudsourcedict)
    src_i, src_j = src_ij
    src_k = findlast(wet3D[src_i, src_j, :])
    bottomcelldepth = (sum(volcello[src_i, src_j, 1:src_k]) - volcello[src_i, src_j, src_k] / 2) / areacello[src_i, src_j]
    bottomcellthickness = volcello[src_i, src_j, src_k] / areacello[src_i, src_j]
    bottomcelltopandbottom = bottomcelldepth .+ [-bottomcellthickness / 2, bottomcellthickness / 2]
    row = (
        lon = lon[src_i, src_j],
        lat = lat[src_i, src_j],
        depth = bottomcelldepth,
        deptho = round(Int, deptho[src_i, src_j]),
        fraction = count / length(cloudsourceijs),
        src_i = src_i,
        src_j = src_j,
        src_k = src_k,
    )
    push!(clouddf, row)
end
@show clouddf




members = map(m -> "r$(m)i1p1f1", 1:40)
# members = map(m -> "r$(m)i1p1f1", 1:3)
Nmembers = length(members)

@info "Grab mean reemergence time at injection locations for all members"
Γouts1cloud = map(members) do member
    @info "loading $member Γ†"
    Γout_file = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)/cyclomonth/reemergence_time$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str).nc"
    Γoutyr4D_ds = open_dataset(Γout_file)
    map(eachrow(clouddf)) do row
        src_i, src_j, src_k = row.src_i, row.src_j, row.src_k
        Γoutyr0D = mean(readcubedata(Γoutyr4D_ds.adjointage[src_i, src_j, src_k, :]))
    end
end


@info "Loading TTDs"
𝒢_files = ["/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)/calgtilde$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str).nc" for member in members]
𝒢_ds = open_mfdataset(DimArray(𝒢_files, Dim{:member}(members)))

@time "Grab 2030s TTDs at injection location" 𝒢s1cloud = map(eachrow(clouddf)) do row
    src_i, src_j = row.src_i, row.src_j
    readcubedata(𝒢_ds.calgtilde[src_i, src_j, :, :]).data
end
𝒢s1cloud = cat(𝒢s1cloud..., dims = 3)


TTD_time = collect(axes(𝒢s1cloud, 1)) # saved every year so index should be same as year
varname = yearly ? "seqeff" : "calE"

# # Concatete into 3D array to facilitate access

@info "Loading sequestration efficiencies"
ℰ_files = ["/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)/$(varname)$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str).nc" for member in members]
ℰ_ds = open_mfdataset(DimArray(ℰ_files, Dim{:member}(members)))
@time "Grab 2030s seq. eff. at injection location" ℰs1cloud = map(eachrow(clouddf)) do row
    src_i, src_j = row.src_i, row.src_j
    readcubedata(ℰ_ds[varname][src_i, src_j, :, :]).data
end
ℰs1cloud = cat(ℰs1cloud..., dims = 3)


# for cool density-colored scatter instead (taken from unpublished website blog, itself from slack convo)
cloudlons = [P.x for P in cloudsourcelocations]
cloudlats = [P.y for P in cloudsourcelocations]
D = kde((cloudlons, cloudlats))
c = [pdf(D, x, y) for (x,y) in zip(cloudlons, cloudlats)]
cidx = sortperm(c)


fig = Figure(size = (800, 600))


#######
# MAP #
#######
# limits = ((107, 158), (-48, -7))
limits = ((170, 200), (-55, -35))
xticks = -180:10:360
yticks = -90:5:90
xticks = (xticks, lonticklabel.(xticks))
yticks = (yticks, latticklabel.(yticks))
panela = Axis(
    fig[1, 1];
    backgroundcolor = :lightgray,
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
colors = cgrad(:jblue, length(depths), categorical = true)[reverse(depths) .≤ maxdepth]
# colors = cgrad(:blues, length(depths), categorical=true)
# colors = cgrad(:grays, length(depths), categorical=true, rev=true)

# isinlimits(P) = !ismissing(P[1]) && !ismissing(P[2]) && (limits[1][1] ≤ P[1] ≤ limits[1][2]) && (limits[2][1] ≤ P[2] ≤ limits[2][2])
# isinlimits(vecP::Vector) = any(isinlimits, vecP)
# for (simplepolygon, color) in zip(reverse(simplepolygons), colors)
#     p = poly!(panela, simplepolygon; color)
#     translate!(p, 0, 0, -100)
# end()
ilon = (limits[1][1] - 10) .< lon.data[:,1] .< (limits[1][2] + 10)
ilat = (limits[2][1] - 10) .< lat.data[1,:] .< (limits[2][2] + 10)
co = contourf!(panela,
    lon.data[ilon, ilat], lat.data[ilon, ilat],
    replace(deptho.data[ilon, ilat], missing => NaN);
    colormap = colors,
    levels = reverse(depths[depths .≤ maxdepth]),
    extendhigh = :auto,
    extendlow = :auto,
)
translate!(co, 0, 0, -90)
# p = poly!(panela2, simplepolygons[2]; color = (:red, 0), strokecolor = :black, strokewidth = 1) # should be 5000m
# translate!(p, 0, 0, -90)
# poly!(panela2, GeoMakie.land(); color = :white, strokecolor = :black, strokewidth = 1)
poly!(panela, GeoMakie.land(); color = :lightgray, strokecolor = :black, strokewidth = 1)

text!(panela, 171.4, -43.5; text = rich("New Zealand", font = :italic), align = (:center, :center), color = :black, fontsize = 14, rotation = pi / 4)

# Add depth colorbar inside
# cbarlabelformat(x) = isinteger(x) ? string(round(Int, x)) : string(x)
# ilow = only(findall(depths .== 0))
# ihigh = only(findall(depths .== maxdepth))
# Colorbar(
#     fig;
#     # limits = (1,length(simplepolygons) - 2),
#     # # ticks = (1:length(simplepolygons) - 2, cbarlabelformat.(1e-3 * reverse(depths[ihigh + 1 : ilow - 1]))),
#     # ticks = (1:length(simplepolygons) - 2, string.(reverse(depths[ihigh + 1 : ilow - 1]))),
#     # colormap = cgrad(reverse(colors)[ihigh + 2 : ilow - 1], categorical = true, rev = true),
#     # highclip = reverse(colors)[ihigh + 1],
#     # lowclip = reverse(colors)[ilow],
#     limits = (1, length(simplepolygons) - 1),
#     # ticks = (1:length(simplepolygons) - 2, cbarlabelformat.(1e-3 * reverse(depths[ihigh + 1 : ilow - 1]))),
#     ticks = (1:(length(simplepolygons) - 1), string.(reverse(depths[ihigh:(ilow - 1)]))),
#     colormap = cgrad(reverse(colors)[(ihigh + 1):(ilow - 1)], categorical = true, rev = true),
#     highclip = reverse(colors)[ihigh],
#     lowclip = reverse(colors)[ilow],
#     label = "seafloor depth (m)",
#     height = 10,
#     width = 160,
#     labelsize = 10,
#     vertical = false,
#     flipaxis = true,
#     ticklabelsize = 10,
#     bbox = panela.scene.viewport,
#     alignmode = Outside(10),
#     valign = :bottom,
#     halign = :left,
#     # topspinecolor = :white,
#     # bottomspinecolor = :white,
#     # rightspinecolor = :white,
#     # leftspinecolor = :white,
#     # ticklabelcolor = :white,
#     # labelcolor = :white,
#     # tickcolor = :white
# )
Colorbar(
    fig, co;
    # limits = (1,length(simplepolygons) - 2),
    # # ticks = (1:length(simplepolygons) - 2, cbarlabelformat.(1e-3 * reverse(depths[ihigh + 1 : ilow - 1]))),
    # ticks = (1:length(simplepolygons) - 2, string.(reverse(depths[ihigh + 1 : ilow - 1]))),
    # colormap = cgrad(reverse(colors)[ihigh + 2 : ilow - 1], categorical = true, rev = true),
    # highclip = reverse(colors)[ihigh + 1],
    # lowclip = reverse(colors)[ilow],
    # limits = (1, length(simplepolygons) - 1),
    # ticks = (1:length(simplepolygons) - 2, cbarlabelformat.(1e-3 * reverse(depths[ihigh + 1 : ilow - 1]))),
    # ticks = (1:(length(simplepolygons) - 1), string.(reverse(depths[ihigh:(ilow - 1)]))),
    # colormap = cgrad(reverse(colors)[(ihigh + 1):(ilow - 1)], categorical = true, rev = true),
    # highclip = reverse(colors)[ihigh],
    # lowclip = reverse(colors)[ilow],
    label = "seafloor depth (m)",
    height = 10,
    width = 160,
    labelsize = 10,
    vertical = false,
    flipaxis = true,
    ticklabelsize = 10,
    bbox = panela.scene.viewport,
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
# colors1 = cgrad(:Egypt, categorical=true)[[3, 4, 1]]
# colors1 = cgrad(:Egypt, categorical=true)[[3]]
colors1 = cgrad(:galah, categorical = true)[[3]]
# color2 = cgrad(:Egypt, categorical=true)[4]
colors2 = cgrad(:galah, categorical = true)[[6]]
# bdcolors1 = [(c + 2 * RGBA(1,1,1,1)) / 3 for c in colors1]
# bdcolors2 = [(c + 2 * RGBA(1,1,1,1)) / 3 for c in colors2]
bdcolors1 = colors1 #[(c + 2 * RGBA(1,1,1,1)) / 3 for c in colors1]
bdcolors2 = colors2 #[(c + 2 * RGBA(1,1,1,1)) / 3 for c in colors2]
Γcolors1 = [(c + 3 * RGBA(1, 1, 1, 1)) / 4 for c in colors1]
Γcolors2 = [(c + 3 * RGBA(1, 1, 1, 1)) / 4 for c in colors2]
# colors = Makie.wong_colors()[[1, 3, 6]]
# offsets = map(x -> x.* 2, [(-2, 1), (-2, -1), (2, -1)])
# offsets = map(x -> x.* 2, [(-2, 1)])
# aligns = [(:right, :bottom), (:right, :top), (:left, :top)]
# aligns = [(:right, :center), (:right, :center), (:left, :center)]
aligns = [(:right, :center)]
# aligns = [(:left, :center)]
# texts = ["A"]
texts = ["injection\nlocation"]

# sc1 = scatter!(panela, cloudlons, cloudlats; marker = :circle, color = colors1[1], markersize = 10, alpha = 0.005)
# Cool density-colored scatter
markersize = 5
cloudsc = scatter!(panela, cloudlons[cidx], cloudlats[cidx]; markersize, color = c[cidx], colormap = cgrad([:black, colors1[1]]))
# Colorbar(fig[2,2], sc; vertical = false, flipaxis = false)
translate!(cloudsc, 0, 0, 80)

# troughcolors = [
#     colorant"#ed683b"
#     colorant"#5465ff"
# ]
# troughlons = [P[1] for P in sourcesBountyTrough]
# troughlats = [P[2] for P in sourcesBountyTrough]
# troughsc = scatter!(panela, troughlons, troughlats; marker = :star5, color = troughcolors, markersize = 20, strokecolor = :black, strokewidth = 1)

# offsets = [
#     (-1, -0.5)
#     (-1, 0.5)
# ]
# troughtxts = [
#     """B ($(round(sourcesBountyTrough[1][1], digits=2))°E, $(-round(sourcesBountyTrough[1][2], digits=2))°S)
#     model depth: $(troughdf.depth[1]) m""",
#     """A ($(round(sourcesBountyTrough[2][1], digits=2))°E, $(-round(sourcesBountyTrough[2][2], digits=2))°S)
#     model depth: $(troughdf.depth[2]) m""",
# ]
# for (itrough, (lon, lat, offset)) in enumerate(zip(troughlons, troughlats, offsets))
#     src_P = (lon, lat)
#     lines!(panela, kinkline(src_P .+ offset, src_P); color = :black)
#     text!(panela, src_P .+ offset; text = troughtxts[itrough], align = (:right, :center), color = :black, strokecolor = :black, offset = (-3, 0))
# end


# for row in eachrow(troughdf)
#     src_i, src_j = row.src_i, row.src_j
#     plotgridcell!(panela, lon_vertices[:, src_i, src_j], lat_vertices[:, src_i, src_j], color = (:white, 0), strokecolor = :white, strokewidth = 0.5)
# end
# for row in eachrow(clouddf)
#     src_i, src_j = row.src_i, row.src_j
#     plotgridcell!(panela, lon_vertices[:, src_i, src_j], lat_vertices[:, src_i, src_j], color = (:white, 0), strokecolor = :white, strokewidth = 0.5)
# end


# for (ksrc, (srcname, offset, align, color, text)) in enumerate(zip(srcnames, offsets, aligns, colors1, texts))
#     src_P = sourcelocation(srcname)
#     # sc1 = scatter!(panela2, src_P; marker=:circle, markersize=10, color=(:black, 0), strokecolor=:black, strokewidth=3)
#     # sc2 = scatter!(panela, src_P; marker=:circle, markersize=10, color=(:black, 0), strokecolor=:black, strokewidth=4)
#     # sc2 = scatter!(panela, src_P; marker=:circle, markersize=10, color=(:black, 0), strokecolor=color, strokewidth=2)
#     # lines!(panela2, [src_P, src_P .+ offset]; color=:white)
#     lines!(panela, kinkline(src_P .+ offset, src_P); color = :black)
#     sc = scatter!(panela, src_P; marker = :star5, markersize = 20, color = colors1[ksrc], strokecolor = :black, strokewidth = 1)
#     # lines!(panela2, kinkline(src_P .+ offset, src_P); color=:black, linewidth=3)
#     # lines!(panela2, kinkline(src_P .+ offset, src_P); color)
#     text!(panela, src_P .+ offset; text, align, color = :black, strokecolor = :black)
#     # text!(panela2, src_P .+ offset; text, align, color=:black, font=:bold, fontsize=18, strokecolor=:black, strokewidth=2)
#     # text!(panela2, src_P .+ offset; text, align, color, font=:bold, fontsize=18)
#     # translate!(sc1, 0, 0, 99)
#     # translate!(sc2, 0, 0, 100)
# end


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
# xmin, xmax = 0, 1500
xmin, xmax = 0, 3000
# ymin, ymax = 0, 1.2
ymin, ymax = 0, 1.0 # It's a pain but this needs to be set for my boxplots to match the size
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

# vspan!(panelb, minimum(minimum.(Γouts1cloud)), maximum(maximum.(Γouts1cloud)), color = Γcolors1[1])
vspan!(panelb, minimum(Γouts1cloud[1, :]), maximum(Γouts1cloud[1, :]), color = (colors[1], 0.2))
vspan!(panelb, minimum(Γouts1cloud[2, :]), maximum(Γouts1cloud[2, :]), color = (colors[2], 0.2))

# 𝒢2030scloudmin = ustrip.(kyr^-1, dropdims(minimum(𝒢s1cloud, dims = (2, 3)), dims = (2, 3)) * s^-1)
# 𝒢2030scloudmax = ustrip.(kyr^-1, dropdims(maximum(𝒢s1cloud, dims = (2, 3)), dims = (2, 3)) * s^-1)
cloudweights = Weights([ustrip(row.fraction) for row in eachrow(clouddf)])
# 𝒢2030scloudmean = ustrip.(kyr^-1, dropdims(mean(mean(𝒢s1cloud, dims = 2), cloudweights, dims = 3), dims = (2,3)) * s^-1)
# bd2030scloud = band!(panelb, TTD_time, 𝒢2030scloudmin, 𝒢2030scloudmax; color = bdcolors1[1], alpha = 0.5)
# ln2030scloud = lines!(panelb, TTD_time, 𝒢2030scloudmean; color = colors1[1], linewidth = 2, linecap = :round, joinstyle = :round)
# text!(panelb, 150, 𝒢2030scloudmax[150]; text = "TTD", align = (:left, :bottom), offset = (3, 1), color = colors1[1])

for (i, row) in enumerate(eachrow(clouddf))
    # 𝒢min = ustrip.(kyr^-1, dropdims(minimum(𝒢s1cloud[:,:,i], dims = (2,)), dims = (2,)) * s^-1)
    # 𝒢max = ustrip.(kyr^-1, dropdims(maximum(𝒢s1cloud[:,:,i], dims = (2,)), dims = (2,)) * s^-1)
    # band!(panelb, TTD_time, 𝒢min, 𝒢max; color = bdcolors1[1], alpha = sqrt(cloudweights[i]) / 2maximum(sqrt.(cloudweights)))
    for m in 1:40
        # lines!(panelb, TTD_time, ustrip.(kyr^-1, 𝒢s1cloud[:, m, i] * s^-1); color = colors1[1], alpha = sqrt(cloudweights[i]) / 2maximum(sqrt.(cloudweights)), linewidth = 1)
        lines!(panelb, TTD_time, ustrip.(kyr^-1, 𝒢s1cloud[:, m, i] * s^-1); color = colors1[1], alpha = cloudweights[i] / 2maximum(cloudweights), linewidth = 1)
    end
end






############################
# sequestration efficiency #
############################
# xmin, xmax = 0, 350
# ymin, ymax = 70, 100
# xmin, xmax = 0, 400
xmin, xmax = 0, 500
# ymin, ymax = 95, 100
ymin, ymax = 80, 100
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
    # xticks = 0:100:xmax,
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

# ℰ2030smin = dropdims(minimum(100 * ℰs1cloud, dims = (2, 3)), dims = (2, 3))
# ℰ2030smax = dropdims(maximum(100 * ℰs1cloud, dims = (2, 3)), dims = (2, 3))
# ℰ2030s = dropdims(mean(mean(100 * ℰs1cloud, dims = 2), cloudweights, dims = 3), dims = (2, 3))
# bd2030s = band!(panelc, TTD_time, ℰ2030smin, ℰ2030smax; color = bdcolors1[1], alpha = 0.5)
# ln2030s = lines!(panelc, TTD_time, ℰ2030s; color = colors1[1], linewidth = 2, linecap = :round, joinstyle = :round)


for (i, row) in enumerate(eachrow(clouddf))
    # 𝒢min = ustrip.(kyr^-1, dropdims(minimum(𝒢s1cloud[:,:,i], dims = (2,)), dims = (2,)) * s^-1)
    # 𝒢max = ustrip.(kyr^-1, dropdims(maximum(𝒢s1cloud[:,:,i], dims = (2,)), dims = (2,)) * s^-1)
    # band!(panelb, TTD_time, 𝒢min, 𝒢max; color = bdcolors1[1], alpha = sqrt(cloudweights[i]) / 2maximum(sqrt.(cloudweights)))
    for m in 1:40
        # lines!(panelb, TTD_time, ustrip.(kyr^-1, 𝒢s1cloud[:, m, i] * s^-1); color = colors1[1], alpha = sqrt(cloudweights[i]) / 2maximum(sqrt.(cloudweights)), linewidth = 1)
        lines!(panelc, TTD_time, 100 * ℰs1cloud[:, m, i]; color = colors1[1], alpha = cloudweights[i] / 2maximum(cloudweights), linewidth = 1)
    end
end
ylims!(panelc, ymax, ymin)






# Repeat for unzoomed time series
xmin, xmax = 0, 3000
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
    yticks = 0:20:100,
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
# vspan!(paneld, minimum(minimum.(Γouts1cloud)), maximum(maximum.(Γouts1cloud)), color = Γcolors1[1])
# vspan!(paneld, minimum(Γouts1trough[1, :]), maximum(Γouts1trough[1, :]), color = (troughcolors[1], 0.2))
# vspan!(paneld, minimum(Γouts1trough[2, :]), maximum(Γouts1trough[2, :]), color = (troughcolors[2], 0.2))

# bd2030s = band!(paneld, TTD_time, ℰ2030smin, ℰ2030smax; color = bdcolors1[1], alpha = 0.5)
# ln2030s = lines!(paneld, TTD_time, ℰ2030s; color = colors1[1], linewidth = 2, linecap = :round, joinstyle = :round)

for (i, row) in enumerate(eachrow(clouddf))
    # 𝒢min = ustrip.(kyr^-1, dropdims(minimum(𝒢s1cloud[:,:,i], dims = (2,)), dims = (2,)) * s^-1)
    # 𝒢max = ustrip.(kyr^-1, dropdims(maximum(𝒢s1cloud[:,:,i], dims = (2,)), dims = (2,)) * s^-1)
    # band!(panelb, TTD_time, 𝒢min, 𝒢max; color = bdcolors1[1], alpha = sqrt(cloudweights[i]) / 2maximum(sqrt.(cloudweights)))
    for m in 1:40
        # lines!(panelb, TTD_time, ustrip.(kyr^-1, 𝒢s1cloud[:, m, i] * s^-1); color = colors1[1], alpha = sqrt(cloudweights[i]) / 2maximum(sqrt.(cloudweights)), linewidth = 1)
        lines!(paneld, TTD_time, 100 * ℰs1cloud[:, m, i]; color = colors1[1], alpha = cloudweights[i] / 2maximum(cloudweights), linewidth = 1)
    end
end
ylims!(paneld, ymax, ymin)





# Add zoom lines
rectattrs = (strokecolor = :lightgray, linestyle = :solid)
lineattrs = (color = :lightgray, linestyle = :solid)
zoom_lines!(panelc, paneld; rectattrs, lineattrs)


colsize!(fig.layout, 1, Relative(0.4))
rowgap!(fig.layout, 20)
colgap!(fig.layout, 20)


# save temp plot
outputdir = joinpath(fixedvarsinputdir, "all_members")
mkpath(outputdir)
outputfile = joinpath(outputdir, "TTD_seqeff_timeseries_Gigablue_final_positions_map.png")
@info "Saving injection location as image file:\n  $(outputfile)"
save(outputfile, fig)
foo





#############
# Box plots #
#############
# (must come after resizing to get correct width of bar plots)
# Box plot of mean age inside
valuesΓ2030strough1 = reduce(vcat, Γouts1trough[1, :])
valuesΓ2030strough2 = reduce(vcat, Γouts1trough[2, :])
categories = reduce(vcat, fill(label, Nmembers) for label in texts)
# Place Γ box plot close to top and stagger its display down for each injection locations (max 4)
categorypositions1 = max𝒢s1trough[1] * (1 .- reduce(vcat, fill(ilabel, Nmembers) for ilabel in reverse(eachindex(texts))) / 4)
categorypositions2 = max𝒢s1trough[2] * (1 .- reduce(vcat, fill(ilabel, Nmembers) for ilabel in reverse(eachindex(texts))) / 4)
color1 = reduce(vcat, fill(color, Nmembers) for color in troughcolors[[1]])
color2 = reduce(vcat, fill(color, Nmembers) for color in troughcolors[[2]])
# I want a 10 pixel width = width_dataspace / limits_dataspace * limits_figspace
# so width_dataspace = 10 * limits_dataspace / limits_figspace
limits_dataspace = panelb.finallimits[].widths[2]
limits_figspace = fullproject(panelb, panelb.finallimits[]).widths[2]
width = 20 * limits_dataspace / limits_figspace
boxplot!(panelb, categorypositions1, valuesΓ2030strough1; color = color1, orientation = :horizontal, width, strokewidth = 1, whiskerwidth = :match)
boxplot!(panelb, categorypositions2, valuesΓ2030strough2; color = color2, orientation = :horizontal, width, strokewidth = 1, whiskerwidth = :match)


meantimestr1 = "$(round(Int, mean(valuesΓ2030strough1))) ± $(round(Int, std(valuesΓ2030strough1))) years"
text!(panelb, mean(valuesΓ2030strough1), categorypositions1[1]; text = rich("location B mean time:\n", meantimestr1), align = (:center, :bottom), offset = (0, +15))
meantimestr2 = "$(round(Int, mean(valuesΓ2030strough2))) ± $(round(Int, std(valuesΓ2030strough2))) years"
text!(panelb, mean(valuesΓ2030strough2), categorypositions2[1]; text = rich("location A mean time:\n", meantimestr2), align = (:center, :top), offset = (0, -15))
# text!(panelb, (minimum(valuesΓ2030strough) + maximum(valuesΓ2030strough)) / 2, categorypositions[1]; text = "2030s", align = (:center, :bottom), offset = (0, 15), fontsize = 10)
# text!(panelb, (minimum(valuesΓ2090strough) + maximum(valuesΓ2090strough)) / 2, categorypositions[1]; text = "2090s", align = (:center, :bottom), offset = (0, 15), fontsize = 10)
# text!(panelb, maximum(values), categorypositions[1]; text = rich("mean time ", Γstr), align = (:left, :center), offset = (5, 0))


mediantimemeantrough1 = mean([findfirst(100 * ℰi .< 50) for ℰi in eachcol(ℰs1trough[:, :, 1])])
mediantimestdtrough1 = std([findfirst(100 * ℰi .< 50) for ℰi in eachcol(ℰs1trough[:, :, 1])])
mediantimemaxtrough1 = maximum([findfirst(100 * ℰi .< 50) for ℰi in eachcol(ℰs1trough[:, :, 1])])
mediantimetroughtxtpos1 = mediantimemaxtrough1
mediantimetroughtxt1 = "$(round(Int, mediantimemeantrough1)) ± $(round(Int, mediantimestdtrough1)) years"
mediantimemeantrough2 = mean([findfirst(100 * ℰi .< 50) for ℰi in eachcol(ℰs1trough[:, :, 2])])
mediantimestdtrough2 = std([findfirst(100 * ℰi .< 50) for ℰi in eachcol(ℰs1trough[:, :, 2])])
mediantimemaxtrough2 = maximum([findfirst(100 * ℰi .< 50) for ℰi in eachcol(ℰs1trough[:, :, 2])])
mediantimetroughtxtpos2 = mediantimemaxtrough2
mediantimetroughtxt2 = "$(round(Int, mediantimemeantrough2)) ± $(round(Int, mediantimestdtrough2)) years"
mediantimetroughtxtpos = [mediantimetroughtxtpos1, mediantimetroughtxtpos2]
mediantimetroughtxts = [mediantimetroughtxt1, mediantimetroughtxt2]

tenthpercentilemeantrough1 = mean([findfirst(100 * ℰi .< 90) for ℰi in eachcol(ℰs1trough[:, :, 1])])
tenthpercentilemaxtrough1 = maximum([findfirst(100 * ℰi .< 90) for ℰi in eachcol(ℰs1trough[:, :, 1])])
tenthpercentilestdtrough1 = std([findfirst(100 * ℰi .< 90) for ℰi in eachcol(ℰs1trough[:, :, 1])])
tenthpercentiletroughtxtpos1 = tenthpercentilemeantrough1
tenthpercentiletroughtxt1 = "$(round(Int, tenthpercentilemeantrough1)) ± $(round(Int, tenthpercentilestdtrough1)) years"
tenthpercentilemeantrough2 = mean([findfirst(100 * ℰi .< 90) for ℰi in eachcol(ℰs1trough[:, :, 2])])
tenthpercentilemaxtrough2 = maximum([findfirst(100 * ℰi .< 90) for ℰi in eachcol(ℰs1trough[:, :, 2])])
tenthpercentilestdtrough2 = std([findfirst(100 * ℰi .< 90) for ℰi in eachcol(ℰs1trough[:, :, 2])])
tenthpercentiletroughtxtpos2 = tenthpercentilemeantrough2
tenthpercentiletroughtxt2 = "$(round(Int, tenthpercentilemeantrough2)) ± $(round(Int, tenthpercentilestdtrough2)) years"
tenthpercentiletxtpos = [tenthpercentilemaxtrough1, tenthpercentilemaxtrough2]
tenthpercentiletxts = [tenthpercentiletroughtxt1, tenthpercentiletroughtxt2]

ℰ100txtpostrough1 = 100 * mean(ℰs1trough[100, :, 1])
ℰ300txtpostrough1 = 100 * mean(ℰs1trough[300, :, 1])
ℰ1000txtpostrough1 = 100 * mean(ℰs1trough[1000, :, 1])
ℰ100txtpostrough2 = 100 * mean(ℰs1trough[100, :, 2])
ℰ300txtpostrough2 = 100 * mean(ℰs1trough[300, :, 2])
ℰ1000txtpostrough2 = 100 * mean(ℰs1trough[1000, :, 2])

ℰ100meantrough1 = mean(100 * ℰs1trough[100, :, 1])
ℰ100stdtrough1 = std(100 * ℰs1trough[100, :, 1])
ℰ100txttrough1 = "$(round(Int, ℰ100meantrough1)) ± $(round(Int, ℰ100stdtrough1)) %"
println("Trough 1 ℰ(100yr) = ", ℰ100txttrough1)
ℰ100meantrough2 = mean(100 * ℰs1trough[100, :, 2])
ℰ100stdtrough2 = std(100 * ℰs1trough[100, :, 2])
ℰ100txttrough2 = "$(round(Int, ℰ100meantrough2)) ± $(round(Int, ℰ100stdtrough2)) %"
println("Trough 2 ℰ(100yr) = ", ℰ100txttrough2)
ℰ300meantrough1 = mean(100 * ℰs1trough[300, :, 1])
ℰ300stdtrough1 = std(100 * ℰs1trough[300, :, 1])
ℰ300txttrough1 = "$(round(Int, ℰ300meantrough1)) ± $(round(Int, ℰ300stdtrough1)) %"
println("Trough 1 ℰ(300yr) = ", ℰ300txttrough1)

ℰ300meantrough2 = mean(100 * ℰs1trough[300, :, 2])
ℰ300stdtrough2 = std(100 * ℰs1trough[300, :, 2])
ℰ300txttrough2 = "$(round(Int, ℰ300meantrough2)) ± $(round(Int, ℰ300stdtrough2)) %"
println("Trough 2 ℰ(300yr) = ", ℰ300txttrough2)
ℰ1000meantrough1 = mean(100 * ℰs1trough[1000, :, 1])
ℰ1000stdtrough1 = std(100 * ℰs1trough[1000, :, 1])
ℰ1000txttrough1 = "$(round(Int, ℰ1000meantrough1)) ± $(round(Int, ℰ1000stdtrough1)) %"
println("Trough 1 ℰ(1000yr) = ", ℰ1000txttrough1)
ℰ1000meantrough2 = mean(100 * ℰs1trough[1000, :, 2])
ℰ1000stdtrough2 = std(100 * ℰs1trough[1000, :, 2])
ℰ1000txttrough2 = "$(round(Int, ℰ1000meantrough2)) ± $(round(Int, ℰ1000stdtrough2)) %"
println("Trough 2 ℰ(1000yr) = ", ℰ1000txttrough2)

ℰ100txttrough = [ℰ100txttrough1, ℰ100txttrough2]
ℰ300txttrough = [ℰ300txttrough1, ℰ300txttrough2]
ℰ1000txttrough = [ℰ1000txttrough1, ℰ1000txttrough2]
ℰ100txtpostrough = [ℰ100txtpostrough1, ℰ100txtpostrough2]
ℰ300txtpostrough = [ℰ300txtpostrough1, ℰ300txtpostrough2]
ℰ1000txtpostrough = [ℰ1000txtpostrough1, ℰ1000txtpostrough2]


troughlabels = ["B", "A"]

for itrough in 1:2
    ℰstrough = ℰs1trough[:, :, itrough]
    troughcolor = troughcolors[itrough]
    troughlabel = troughlabels[itrough]

    for panel in (panelc, paneld)
        # given ℰ
        for ℰval in [50, 90]
            values = [findfirst(100 * ℰi .< ℰval) for ℰi in eachcol(ℰstrough)]
            local categorypositions = fill(ℰval, length(values))
            local color = fill(troughcolor, length(values))
            local limits_dataspace = panel.finallimits[].widths[2]
            local limits_figspace = abs(fullproject(panel, panel.finallimits[]).widths[2])
            local width = 20 * limits_dataspace / limits_figspace
            boxplot!(panel, categorypositions, values; color, orientation = :horizontal, width, strokewidth = 1, whiskerwidth = :match)
            # (ℰval == 50) && (panel == paneld) && text!(panel, mediantimetroughtxtpos[itrough], 50; text = "$(mediantimetroughtxts[itrough])", align = (:left, :center), offset = (10, 0), fontsize=12, color = troughcolor)
            # (ℰval == 90) && (panel == panelc) && text!(panel, tenthpercentiletxtpos[itrough], 90; text = "$(tenthpercentiletxts[itrough])", align = (:left, :center), offset = (15, 0), fontsize=12, color = troughcolor)
            txtpos = mediantimetroughtxtpos[itrough]
            (ℰval == 50) && (panel == panelc) && (txtpos < 200) && text!(panel, txtpos, 50; text = "$(troughlabel) median time:\n$(mediantimetroughtxts[itrough])", align = (:left, :center), offset = (10, 0), fontsize=12)
            (ℰval == 50) && (panel == paneld) && (txtpos > 200) && text!(panel, txtpos, 50; text = "$(troughlabel) median time:\n$(mediantimetroughtxts[itrough])", align = (:left, :center), offset = (10, 0), fontsize=12)
            txtpos = tenthpercentiletxtpos[itrough]
            (ℰval == 90) && (panel == panelc) && (txtpos < 200) && text!(panel, txtpos, 90; text = "$(troughlabel) 10%ile time:\n$(tenthpercentiletxts[itrough])", align = (:center, :bottom), offset = (0, 15), fontsize=12)
            (ℰval == 90) && (panel == paneld) && (txtpos > 200) && text!(panel, txtpos, 90; text = "$(troughlabel) 10%ile time:\n$(tenthpercentiletxts[itrough])", align = (:center, :bottom), offset = (0, 15), fontsize=12)
        end
        # given τ
        for τval in [100, 300, 1000]
            values = [100 * ℰi[τval] for ℰi in eachcol(ℰstrough)]
            local categorypositions = fill(τval, length(values))
            local color = fill(troughcolor, length(values))
            local limits_dataspace = panel.finallimits[].widths[1]
            local limits_figspace = fullproject(panel, panel.finallimits[]).widths[1]
            local width = 20 * limits_dataspace / limits_figspace
            # boxplot!(panel, categorypositions, values; color, orientation = :vertical, width, strokewidth = 1, whiskerwidth = :match)
            # (τval == 300) && (panel == panelc) && (ℰsdecade == ℰs1cloud) && text!(panel, 300, ℰ300txtpos; text = rich(ℰstr, "(300 yr)"), align = (:center, :bottom), offset = (-15, 0), rotation = π/2)
            # (τval == 1000) && (panel == paneld) && (ℰsdecade == ℰs1cloud) && text!(panel, 1000, ℰ1000txtpos; text = rich(ℰstr, "(1000 yr)"), align = (:center, :bottom), offset = (-15, 0), rotation = π/2)
            # scatter!(panel, fill(τval, length(values)), values)
        end
    end
end
# band for mean age


labeloptions = (
    font = :bold,
    align = (:left, :top),
    offset = (5, -2),
    space = :relative,
    fontsize = 24,
)
for (ax, label) in zip([panela, panelb, panelc, paneld], string.('a':'h'))
    text!(ax, 0, 1; text = label, labeloptions..., strokecolor = :white, strokewidth = 3)
    text!(ax, 0, 1; text = label, labeloptions...)
end


# save plot
outputdir = joinpath(fixedvarsinputdir, "all_members")
mkpath(outputdir)
outputfile = joinpath(outputdir, "TTD_seqeff_timeseries_1loc_GigablueCloud.png")
@info "Saving injection location as image file:\n  $(outputfile)"
save(outputfile, fig)
outputfile = joinpath(outputdir, "TTD_seqeff_timeseries_1loc_GigablueCloud.pdf")
@info "Saving injection location as image file:\n  $(outputfile)"
save(outputfile, fig)

# Save the data for Gigablue
metadata = Dict(
    "description" => """
        Sequestration efficiency and TTD as plotted in the figure requested by Gigablue.
        Note that the `location_index` corresponds to the index of the particle bins.
        That is, the particle injection locations have been binned onto the model grid,
        and each `location_index` corresponds to one of these grid cells where at least
        one particle was binned.
        In the figure, the fraction of particles binned in each grid cell was used
        as an opacity weight for each curve.
        Made by Benoît Pasquier Thursday 15 Jan 2026.
    """,
    "model" => model,
    "experiment" => experiment,
    "Ti unit" => "yr",
    "Sequestration efficiency unit" => "%",
    "Transit time distribution unit" => "kyr^-1",
    "Dimensions" => "(time, member, location_index)",
)
axlist = (
    dims(DimArray(ones(size(TTD_time)), Ti(TTD_time)))[1],
    dims(DimArray(ones(size(members)), Dim{:member}(members)))[1],
    dims(DimArray(ones(size(clouddf, 1)), Dim{:location_index}(1:size(clouddf, 1))))[1],
)
TTDproperties = Dict(
    "description" => "Transit time distribution",
    "unit" => "kyr^-1",
)
seqeffproperties = Dict(
    "description" => "Sequestration efficiency",
    "unit" => "%",
)
TTDcube = YAXArray(axlist, ustrip.(kyr^-1, 𝒢s1cloud * s^-1), TTDproperties)
seqeffcube = YAXArray(axlist, 100 * ℰs1cloud, seqeffproperties)

lonproperties = Dict(
    "description" => "Longitude of injection location grid cell center",
    "unit" => "degrees_east",
)
loncube = YAXArray(axlist[[3]], clouddf.lon, lonproperties)
latproperties = Dict(
    "description" => "Latitude of injection location grid cell center",
    "unit" => "degrees_north",
)
latcube = YAXArray(axlist[[3]], clouddf.lat, latproperties)
depthproperties = Dict(
    "description" => "Depth of injection location grid cell center",
    "unit" => "m",
)
depthcube = YAXArray(axlist[[3]], clouddf.depth, depthproperties)
fractionproperties = Dict(
    "description" => "Fraction of particles binned in grid cell",
    "unit" => "%",
)
fractioncube = YAXArray(axlist[[3]], 100 * clouddf.fraction, fractionproperties)

arrays = Dict(
    :TTD => TTDcube,
    :seqeff => seqeffcube,
    :lon => loncube,
    :lat => latcube,
    :depth => depthcube,
    :fraction => fractioncube,
)
ds = Dataset(; properties = metadata, arrays...)

# Save to netCDF file
outputfile = joinpath(outputdir, "Pasquier_data_for_GigaBlueCloud.nc")
@info "Saving TTD and sequestration efficiency time series as netCDF file:\n  $(outputfile)"
# ds_chunked = setchunks(ds, (x = 60, y = 60, Ti = length(ds.Ti)))
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)
