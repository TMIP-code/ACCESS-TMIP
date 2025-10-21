# qsub -I -P y99 -q express -l mem=47GB -l storage=scratch/gh0+scratch/xv83+scratch/p66 -l walltime=01:00:00 -l ncpus=12

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
using Format
using KernelDensity

include("plotting_functions.jl")

#############
# Load data #
#############


# script options
@show model = "ACCESS-OM2-025"
if isempty(ARGS)
    experiment = "omip2"
    member = "r1i1p1f1"
    time_window = "Jan0200-Dec0209"
else
    experiment, member, time_window = ARGS
end

@show experiment
@show member
@show time_window

# preferred diffusivities
κVdeep = 3.0e-5 # m^2/s
κVML = 1.0      # m^2/s
κH = 300.0 / 4  # m^2/s (grid-scaling by sqrt(area))
@show κVdeep
@show κVML
@show κH
κVdeep_str = "kVdeep" * format(κVdeep, conversion="e")
κVML_str = "kVML" * format(κVML, conversion="e")
κH_str = "kH" * format(κH, conversion="d")

upwind = false
@show upwind
upwind_str = upwind ? "" : "_centered"
upwind_str2 = upwind ? "upwind" : "centered"

# Load areacello and volcello for grid geometry
inputdir = "/scratch/y99/TMIP/data/ACCESS-OM2-025/omip2/r1i1p1f1/Jan0200-Dec0209/"
areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))
volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))

# Load fixed variables in memory
areacello = readcubedata(areacello_ds.areacello)
volcello = readcubedata(volcello_ds.volcello)
lon = readcubedata(volcello_ds.lon)
lat = readcubedata(volcello_ds.lat)
lev = volcello_ds.lev
# Identify the vertices keys (vary across CMIPs / models)
# FIXME using test CMORized data for vertices. Hopefully they match!
CMORtestfile = "/scratch/p66/yz9299/OM2_CMORised/umo_Omon_ACCESS-OM2-025_omip1_r1i1p1f1_gn_190001-190112.nc"
CMORtest_ds = open_dataset(CMORtestfile)
volcello_keys = propertynames(CMORtest_ds)
lon_vertices_key = volcello_keys[findfirst(x -> occursin("lon", x) & occursin("vert", x), string.(volcello_keys))]
lat_vertices_key = volcello_keys[findfirst(x -> occursin("lat", x) & occursin("vert", x), string.(volcello_keys))]
lon_vertices = readcubedata(getproperty(CMORtest_ds, lon_vertices_key))
lat_vertices = readcubedata(getproperty(CMORtest_ds, lat_vertices_key))


# Make makegridmetrics
gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
(; lon_vertices, lat_vertices, v3D) = gridmetrics

# Make indices
indices = makeindices(v3D)
(; N, wet3D) = indices

agefile1 = joinpath(inputdir, "steady_age_$(κVdeep_str)_$(κH_str)_$(κVML_str).nc")
age1 = readcubedata(open_dataset(joinpath(inputdir, agefile1)).age)
agefile2 = joinpath(inputdir, "steady_age_$(κVdeep_str)_$(κH_str)_$(κVML_str)_OM2-025-2x2.nc")
age2 = readcubedata(open_dataset(joinpath(inputdir, agefile2)).age)

OCEANS = OceanBasins.oceanpolygons()

zt = lev |> Array


##################
# Zonal averages #
##################



basin_keys = (:ATL, :PAC, :IND)
basin_strs = ("Atlantic", "Pacific", "Indian")
basin_functions = (isatlantic, ispacific, isindian)
basin_values = (reshape(f(lat[:], lon[:], OCEANS), size(lat)) for f in basin_functions)
basins = (; (basin_keys .=> basin_values)...)
basin_latlims_values = [clamp.((-5, +5) .+ extrema(lat[.!isnan.(v3D[:,:,1]) .& basin[:,:,1]]), -80, 80) for basin in basins]
basin_latlims = (; (basin_keys .=> basin_latlims_values)...)
basin_sumlatranges = sum(x[2] - x[1] for x in basin_latlims_values)

contouroptions1 = let
    levels = 0:100:1600
    colormap = cgrad(:viridis, length(levels); categorical=true)
    extendlow = nothing
    extendhigh = colormap[end]
    colormap = cgrad(colormap[1:end-1]; categorical=true)
    nan_color = :lightgray
    (; levels, colormap, extendlow, extendhigh, nan_color)
end
contouroptionsdiff = let
    levels = -100:10:100
    colormap = cgrad(:balance, length(levels); categorical=true)[[1:end÷2+1; end÷2+1:end]]
    extendlow = colormap[1]
    extendhigh = colormap[end]
    colormap = cgrad(colormap[2:end-1]; categorical=true)
    nan_color = :lightgray
    (; levels, colormap, extendlow, extendhigh, nan_color)
end
yticks = 0:1000:6000
yticks = (yticks, [t < 6000 ? string(t) : "" for t in yticks])
axisoptions = (;
    backgroundcolor = :lightgray,
    xgridvisible = true,
    ygridvisible = true,
    ylabel = "depth (m)",
    yticks,
)

data = (age1, age2, age2 - age1)
strs = [model, "$model\ncoarsened", "Difference"]
Nrows = length(strs)
Ncols = length(basins)
fig = Figure(size = (3 * basin_sumlatranges, 250 * Nrows), fontsize = 18)
axs = Array{Any,2}(undef, (Nrows, Ncols))
contours = Array{Any,2}(undef, (Nrows, Ncols))

lat2 = dropdims(maximum(lat, dims=1), dims=1) |> Array # <- for plotting ZAVG (inexact)

for (irow, (x3D, str)) in enumerate(zip(data, strs))

    for (icol, (basin_key, basin)) in enumerate(pairs(basins))

        x2D = zonalaverage(x3D |> Array, gridmetrics; mask = basin)

        local ax = Axis(fig[irow, icol]; axisoptions...)

        contouroptions = irow < 3 ? contouroptions1 : contouroptionsdiff
        local co = contourf!(ax, lat2, zt, x2D; contouroptions...)

        translate!(co, 0, 0, -100)
        contours[irow, icol] = co

        xlim = basin_latlims[basin_key]
        # basin2 = LONGTEXT[basin]

        ax.yticks = (ztick, zticklabel)
        xticks = -90:30:90
        ax.xticks = (xticks, latticklabel.(xticks))
        ylims!(ax, zlim)
        # xlims!(ax, (-90, 90))
        xlims!(ax, xlim)

        myhidexdecorations!(ax, irow < Nrows)
        myhideydecorations!(ax, icol > 1)

        axs[irow, icol] = ax
    end

    Label(fig[irow, 0], text = str, fontsize=20, tellheight=false, rotation=π/2)

end

cb1 = Colorbar(fig[1:Nrows - 1, Ncols + 1], contours[1, 1];
    vertical = true, flipaxis = true,
    # ticks = (, cbarticklabelformat.(levels)),
    label = rich("age (years)"),
)
cb1.height = Relative(0.6)

cbdiff = Colorbar(fig[Nrows, Ncols + 1], contours[Nrows, 1];
    vertical = true, flipaxis = true,
    # ticks = (, cbarticklabelformat.(levels)),
    label = rich("Δage (years)"),
    tickformat = x -> divergingcbarticklabel.(x),
)
cbdiff.height = Relative(0.8)

for (icol, (basin_str, xlims)) in enumerate(zip(basin_strs, basin_latlims))
    Label(fig[0, icol], basin_str, fontsize=20, tellwidth=false)
    colsize!(fig.layout, icol, Auto(xlims[2] - xlims[1]))
end

title = "Effect of coarsening on ideal age from $model (0.25°) to (0.5°)"

Label(fig[-1, 1:3], text = title, fontsize=20, tellwidth=false)
labels = permutedims(reshape(string.('a':'a' + length(axs) - 1), size(axs')))
labeloptions = (
    font = :bold,
    align = (:left, :bottom),
    offset = (5, 2),
    space = :relative,
    fontsize = 24
)
for (ax, label) in zip(axs, labels)
    text!(ax, 0, 0; text = label, labeloptions..., strokecolor = :white, strokewidth = 3)
    text!(ax, 0, 0; text = label, labeloptions...)
end




rowgap!(fig.layout, 20)
colgap!(fig.layout, 20)
rowgap!(fig.layout, 1, 10)
colgap!(fig.layout, 1, 10)

# save plot
outputfile = joinpath(inputdir, "coarsening_effect_ZAVGs.png")
@info "Saving as image file:\n  $(outputfile)"
save(outputfile, fig)




#####################
# Meridional slices #
#####################

# Redo the lon bands from Chamberlain et al. 2019
basin_keys = (:ATL, :PAC)
basin_strs = ("Atlantic 30–40°W", "Pacific 170–180°W")
isatlanticband(lat, lon, OCEANS) = isatlantic(lat, lon, OCEANS) .& (320 .≤ mod.(lon, 360) .≤ 330)
ispacificband(lat, lon, OCEANS) = ispacific(lat, lon, OCEANS) .& (180 .≤ mod.(lon, 360) .≤ 190)
basin_functions = (isatlanticband, ispacificband)
basin_values = (reshape(f(lat[:], lon[:], OCEANS), size(lat)) for f in basin_functions)
basins = (; (basin_keys .=> basin_values)...)
basin_latlims_values = [clamp.((-5, +5) .+ extrema(lat[.!isnan.(v3D[:,:,1]) .& basin[:,:,1]]), -80, 80) for basin in basins]
basin_latlims = (; (basin_keys .=> basin_latlims_values)...)
basin_sumlatranges = sum(x[2] - x[1] for x in basin_latlims_values)


axisoptions = (
    backgroundcolor = :lightgray,
    xgridvisible = true,
    ygridvisible = true,
    ylabel = "depth (m)",
    yticks,
)


Nrows = length(strs)
Ncols = length(basins)
fig = Figure(size = (3 * basin_sumlatranges, 250 * Nrows), fontsize = 18)
axs = Array{Any,2}(undef, (Nrows, Ncols))
contours = Array{Any,2}(undef, (Nrows, Ncols))

lat2 = dropdims(maximum(lat, dims=1), dims=1) |> Array # <- for plotting ZAVG (inexact)

for (irow, (x3D, str)) in enumerate(zip(data, strs))

    for (icol, (basin_key, basin)) in enumerate(pairs(basins))

        x2D = zonalaverage(x3D |> Array, gridmetrics; mask = basin)

        local ax = Axis(fig[irow, icol]; axisoptions...)

        contouroptions = irow < 3 ? contouroptions1 : contouroptionsdiff
        local co = contourf!(ax, lat2, zt, x2D; contouroptions...)

        translate!(co, 0, 0, -100)
        contours[irow, icol] = co

        xlim = basin_latlims[basin_key]
        # basin2 = LONGTEXT[basin]

        ax.yticks = (ztick, zticklabel)
        xticks = -90:30:90
        ax.xticks = (xticks, latticklabel.(xticks))
        ylims!(ax, zlim)
        # xlims!(ax, (-90, 90))
        xlims!(ax, xlim)

        myhidexdecorations!(ax, irow < Nrows)
        myhideydecorations!(ax, icol > 1)

        axs[irow, icol] = ax
    end

    Label(fig[irow, 0], text = str, fontsize=20, tellheight=false, rotation=π/2)

end

cb1 = Colorbar(fig[1:Nrows - 1, Ncols + 1], contours[1, 1];
    vertical = true, flipaxis = true,
    # ticks = (, cbarticklabelformat.(levels)),
    label = rich("age (years)"),
)
cb1.height = Relative(0.6)

cbdiff = Colorbar(fig[Nrows, Ncols + 1], contours[Nrows, 1];
    vertical = true, flipaxis = true,
    # ticks = (, cbarticklabelformat.(levels)),
    label = rich("Δage (years)"),
    tickformat = x -> divergingcbarticklabel.(x),
)
cbdiff.height = Relative(0.8)

for (icol, (basin_str, xlims)) in enumerate(zip(basin_strs, basin_latlims))
    Label(fig[0, icol], basin_str, fontsize=20, tellwidth=false)
    colsize!(fig.layout, icol, Auto(xlims[2] - xlims[1]))
end

title = "Effect of coarsening on ideal age from $model (0.25°) to (0.5°)"
Label(fig[-1, 1:3], text = title, fontsize=20, tellwidth=false)

labels = permutedims(reshape(string.('a':'a' + length(axs) - 1), size(axs')))
labeloptions = (
    font = :bold,
    align = (:left, :bottom),
    offset = (5, 2),
    space = :relative,
    fontsize = 24
)
for (ax, label) in zip(axs, labels)
    text!(ax, 0, 0; text = label, labeloptions..., strokecolor = :white, strokewidth = 3)
    text!(ax, 0, 0; text = label, labeloptions...)
end


rowgap!(fig.layout, 20)
colgap!(fig.layout, 20)
rowgap!(fig.layout, 1, 10)
colgap!(fig.layout, 1, 10)

# save plot
outputfile = joinpath(inputdir, "coarsening_effect_slices.png")
@info "Saving ideal age ZAVGs as image file:\n  $(outputfile)"
save(outputfile, fig)