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

include("plotting_functions.jl")

model = "ACCESS-ESM1-5"
# model = "ACCESS-CM2"
# model = "ACCESS1-3"

# CMIP_version = "CMIP5"
CMIP_version = "CMIP6"

experiment = "historical"
# experiment = "piControl"

time_window = "Jan1990-Dec1999"
# time_window = "Jan1071-Dec1100" # <- last 30 years of ACCESS-ESM1-5 piControl
# time_window = "Jan1420-Dec1449" # <- last 30 years of ACCESS-CM2 piControl


# Gadi directory for input files
# inputdirfun(member) = "/scratch/xv83/TMIP/data/$model/$experiment/all members/$(time_window)"
inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/all_members/$(time_window)/cyclomonth"
outputdir = inputdir
mkpath(inputdir)


gridinputdir = "/scratch/xv83/TMIP/data/$model/$experiment/r1i1p1f1/$(time_window)"
areacello_ds = open_dataset(joinpath(gridinputdir, "areacello.nc"))
volcello_ds = open_dataset(joinpath(gridinputdir, "volcello.nc"))
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
(; lon_vertices, lat_vertices, lon, lat, zt, v3D,) = gridmetrics
lev = zt
# Make indices
indices = makeindices(gridmetrics.v3D)
(; wet3D, N) = indices



Γdown = rich("Γ", superscript("↓"))
Γup = rich("Γ", superscript("↑"))






# Plot zonal averages

basin_keys = (:ATL, :PAC, :IND)
basin_strs = ("Atlantic", "Pacific", "Indian")
basin_functions = (isatlantic, ispacific, isindian)
basin_values = (reshape(f(lat[:], lon[:], OCEANS), size(lat)) for f in basin_functions)
basins = (; (basin_keys .=> basin_values)...)
basin_latlims_values = [clamp.((-5, +5) .+ extrema(lat[.!isnan.(v3D[:,:,1]) .& basin[:,:,1]]), -80, 80) for basin in basins]
basin_latlims = (; (basin_keys .=> basin_latlims_values)...)


Γinyr3D_mean = readcubedata(open_dataset(joinpath(inputdir, "age_ensemblemean.nc")).age_ensemblemean)
Γinyr3D_std = readcubedata(open_dataset(joinpath(inputdir, "age_ensemblestd.nc")).age_ensemblestd)
Γinyr3D_max = readcubedata(open_dataset(joinpath(inputdir, "age_ensemblemax.nc")).age_ensemblemax)
Γinyr3D_min = readcubedata(open_dataset(joinpath(inputdir, "age_ensemblemin.nc")).age_ensemblemin)
Γinyr3D_maxdiff = Γinyr3D_max - Γinyr3D_min
Γoutyr3D_mean = readcubedata(open_dataset(joinpath(inputdir, "adjointage_ensemblemean.nc")).adjointage_ensemblemean)
Γoutyr3D_std = readcubedata(open_dataset(joinpath(inputdir, "adjointage_ensemblestd.nc")).adjointage_ensemblestd)
Γoutyr3D_max = readcubedata(open_dataset(joinpath(inputdir, "adjointage_ensemblemax.nc")).adjointage_ensemblemax)
Γoutyr3D_min = readcubedata(open_dataset(joinpath(inputdir, "adjointage_ensemblemin.nc")).adjointage_ensemblemin)
Γoutyr3D_maxdiff = Γoutyr3D_max - Γoutyr3D_min

# Plot Γ↓ zonal averages

fig = Figure(size = (1200, 600), fontsize = 18)
axs = Array{Any,2}(undef, (2, 3))
contours = Array{Any,2}(undef, (2, 3))

for (irow, (x3D, str)) in enumerate(zip((Γinyr3D_mean, Γinyr3D_std), ("mean", "std")))

    if str == "mean" # mean
        levels = 0:100:1500
        colormap = cgrad(:viridis, length(levels); categorical=true)
        extendlow = nothing
        extendhigh = colormap[end]
        colormap = cgrad(colormap[1:end-1]; categorical=true)
    elseif str == "std"
        levels = 0:50:400
        colormap = cgrad(:magma, length(levels); categorical=true)
        extendlow = nothing
        extendhigh = colormap[end]
        colormap = cgrad(colormap[1:end-1]; categorical=true)
    end

    for (icol, (basin_key, basin)) in enumerate(pairs(basins))

        x2D = zonalaverage(x3D, gridmetrics; mask = basin)

        local ax = Axis(fig[irow, icol],
            backgroundcolor=:lightgray,
            xgridvisible=false, ygridvisible=false,
            ylabel = "depth (m)")

        X = dropdims(maximum(lat, dims=1), dims=1)
        Y = zt
        Z = x2D

        co = contourf!(ax, X, Y, Z;
            levels,
            colormap,
            nan_color = :lightgray,
            extendlow,
            extendhigh,
        )
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

        hidexdecorations!(ax,
            label = irow < 2, ticklabels = irow < 2,
            ticks = irow < 2, grid = false)
        hideydecorations!(ax,
            label = icol > 1, ticklabels = icol > 1,
            ticks = icol > 1, grid = false)


        axs[irow, icol] = ax

    end

    cb = Colorbar(fig[irow, 4], contours[irow, 1];
        vertical = true, flipaxis = true,
        # ticks = (, cbarticklabelformat.(levels)),
        label = rich(str, " ", Γdown, " (yr)"),
        )
    cb.height = Relative(1)
end


for (icol, (basin_str, xlims)) in enumerate(zip(basin_strs, basin_latlims))
    Label(fig[0, icol], basin_str, fontsize=20, tellwidth=false)
    colsize!(fig.layout, icol, Auto(xlims[2] - xlims[1]))
end

title = "$model $experiment $(time_window) ideal age"
Label(fig[-1, 1:3], text = title, fontsize=20, tellwidth=false)

rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)

# save plot
outputfile = joinpath(outputdir, "ideal_age_ZAVGs_v4.png")
@info "Saving ideal age ZAVGs as image file:\n  $(outputfile)"
save(outputfile, fig)






# Plot Γ↑ zonal averages


fig = Figure(size = (1200, 600), fontsize = 18)
axs = Array{Any,2}(undef, (2, 3))
contours = Array{Any,2}(undef, (2, 3))

for (irow, (x3D, str)) in enumerate(zip((Γoutyr3D_mean, Γoutyr3D_std), ("mean", "std")))

    if str == "mean" # mean
        levels = 0:100:1500
        colormap = cgrad(:viridis, length(levels); categorical=true)
        extendlow = nothing
        extendhigh = colormap[end]
        colormap = cgrad(colormap[1:end-1]; categorical=true)
    elseif str == "std"
        levels = 0:50:400
        colormap = cgrad(:magma, length(levels); categorical=true)
        extendlow = nothing
        extendhigh = colormap[end]
        colormap = cgrad(colormap[1:end-1]; categorical=true)
    end

    for (icol, (basin_key, basin)) in enumerate(pairs(basins))

        x2D = zonalaverage(x3D, gridmetrics; mask = basin)

        local ax = Axis(fig[irow, icol],
            backgroundcolor=:lightgray,
            xgridvisible=false, ygridvisible=false,
            ylabel = "depth (m)")

        X = dropdims(maximum(lat, dims=1), dims=1)
        Y = zt
        Z = x2D
        co = contourf!(ax, X, Y, Z;
            levels,
            colormap,
            nan_color = :lightgray,
            extendlow,
            extendhigh,
        )
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

        hidexdecorations!(ax,
            label = irow < 2, ticklabels = irow < 2,
            ticks = irow < 2, grid = false)
        hideydecorations!(ax,
            label = icol > 1, ticklabels = icol > 1,
            ticks = icol > 1, grid = false)

        axs[irow, icol] = ax
    end
    cb = Colorbar(fig[irow, 4], contours[irow, 1];
        vertical = true, flipaxis = true,
        # ticks = (, cbarticklabelformat.(levels)),
        label = rich(str, " ", Γup, " (yr)"),
        )
    cb.height = Relative(1)

end



for (icol, (basin_str, xlims)) in enumerate(zip(basin_strs, basin_latlims))
    Label(fig[0, icol], basin_str, fontsize=20, tellwidth=false)
    colsize!(fig.layout, icol, Auto(xlims[2] - xlims[1]))
end

title = "$model $experiment $(time_window) reemergence time"
Label(fig[-1, 1:3], text = title, fontsize=20, tellwidth=false)

rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)

# save plot
outputfile = joinpath(outputdir, "reemergence_time_ZAVGs_v4.png")
@info "Saving reemergence time ZAVGs as image file:\n  $(outputfile)"
save(outputfile, fig)









fig = Figure(size = (1200, 1200), fontsize = 18)
axs = Array{Any,2}(undef, (2, 1))
contours = Array{Any,2}(undef, (2, 1))

for (irow, (x3D, str)) in enumerate(zip((Γinyr3D_mean, Γinyr3D_std), ("mean", "std")))

    if str == "mean" # mean
        levels = 0:100:1500
        colormap = :viridis
        colorrange = extrema(levels)
    elseif str == "std"
        levels = 0:50:400
        colorrange = extrema(levels)
        colormap = :magma
    end

    # Plot mean age at the seafloor level
    local title = "$model $experiment $(time_window) mean age at seafloor"

    ax = Axis(fig[irow,1]; title, xtickformat, ytickformat)


    # plot
    x2D = seafloorvalue(x3D, wet3D)
    plt = plotmap!(ax, x2D, gridmetrics; colorrange, colormap)

    Colorbar(fig[irow,2], plt, label=rich(str, " ", Γdown, " at seafloor (yr)"))
end

rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)

# save plot
outputfile = joinpath(outputdir, "mean_age_at_seafloor_v4.png")
@info "Saving ideal mean age at sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)















fig = Figure(size = (1200, 1200), fontsize = 18)
axs = Array{Any,2}(undef, (2, 1))
contours = Array{Any,2}(undef, (2, 1))

for (irow, (x3D, str)) in enumerate(zip((Γoutyr3D_mean, Γoutyr3D_std), ("mean", "std")))

    if str == "mean" # mean
        levels = 0:100:1500
        colormap = :viridis
        colorrange = extrema(levels)
    elseif str == "std"
        levels = 0:50:400
        colorrange = extrema(levels)
        colormap = :magma
    end

    # Plot mean age at the seafloor level
    local title = "$model $experiment $(time_window) reemergence time at seafloor"

    ax = Axis(fig[irow,1]; title, xtickformat, ytickformat)


    # plot
    x2D = seafloorvalue(x3D, wet3D)
    plt = plotmap!(ax, x2D, gridmetrics; colorrange, colormap)

    # poly!(ax, reverse.(OCEANS[OceanBasins.atlantic()].polygon))
    # poly!(ax, reverse.(OCEANS[OceanBasins.indian()].polygon))
    # poly!(ax, reverse.(OCEANS[OceanBasins.east_pacific()].polygon))
    # poly!(ax, reverse.(OCEANS[OceanBasins.west_pacific()].polygon))

    Colorbar(fig[irow,2], plt, label=rich(str, " ", Γup, " at seafloor (yr)"))
end

rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)

# save plot
outputfile = joinpath(outputdir, "reemergence_time_at_seafloor_v4.png")
@info "Saving mean reemergence time at sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)











# Redo the same for max minus min


# Plot Γ↓ zonal averages

fig = Figure(size = (1200, 600), fontsize = 18)
axs = Array{Any,2}(undef, (2, 3))
contours = Array{Any,2}(undef, (2, 3))

for (irow, (x3D, str)) in enumerate(zip((Γinyr3D_mean, Γinyr3D_maxdiff), ("mean", "maxΔ")))

    if str == "mean" # mean
        levels = 0:100:1500
        colormap = cgrad(:viridis, length(levels); categorical=true)
        extendlow = nothing
        extendhigh = colormap[end]
        colormap = cgrad(colormap[1:end-1]; categorical=true)
    elseif str == "maxΔ"
        levels = 0:50:400
        colormap = cgrad(:magma, length(levels); categorical=true)
        extendlow = nothing
        extendhigh = colormap[end]
        colormap = cgrad(colormap[1:end-1]; categorical=true)
    end

    for (icol, (basin_key, basin)) in enumerate(pairs(basins))

        x2D = zonalaverage(x3D, gridmetrics; mask = basin)

        local ax = Axis(fig[irow, icol],
            backgroundcolor=:lightgray,
            xgridvisible=false, ygridvisible=false,
            ylabel = "depth (m)")

        X = dropdims(maximum(lat, dims=1), dims=1)
        Y = zt
        Z = x2D

        co = contourf!(ax, X, Y, Z;
            levels,
            colormap,
            nan_color = :lightgray,
            extendlow,
            extendhigh,
        )
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

        hidexdecorations!(ax,
            label = irow < 2, ticklabels = irow < 2,
            ticks = irow < 2, grid = false)
        hideydecorations!(ax,
            label = icol > 1, ticklabels = icol > 1,
            ticks = icol > 1, grid = false)


        axs[irow, icol] = ax

    end

    cb = Colorbar(fig[irow, 4], contours[irow, 1];
        vertical = true, flipaxis = true,
        # ticks = (, cbarticklabelformat.(levels)),
        label = rich(str, " ", Γdown, " (yr)"),
        )
    cb.height = Relative(1)
end


for (icol, (basin_str, xlims)) in enumerate(zip(basin_strs, basin_latlims))
    Label(fig[0, icol], basin_str, fontsize=20, tellwidth=false)
    colsize!(fig.layout, icol, Auto(xlims[2] - xlims[1]))
end

title = "$model $experiment $(time_window) ideal age"
Label(fig[-1, 1:3], text = title, fontsize=20, tellwidth=false)

rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)

# save plot
outputfile = joinpath(outputdir, "ideal_age_ZAVGs_vsmaxdiff.png")
@info "Saving ideal age ZAVGs as image file:\n  $(outputfile)"
save(outputfile, fig)






# Plot Γ↑ zonal averages


fig = Figure(size = (1200, 600), fontsize = 18)
axs = Array{Any,2}(undef, (2, 3))
contours = Array{Any,2}(undef, (2, 3))

for (irow, (x3D, str)) in enumerate(zip((Γoutyr3D_mean, Γoutyr3D_maxdiff), ("mean", "maxΔ")))

    if str == "mean" # mean
        levels = 0:100:1500
        colormap = cgrad(:viridis, length(levels); categorical=true)
        extendlow = nothing
        extendhigh = colormap[end]
        colormap = cgrad(colormap[1:end-1]; categorical=true)
    elseif str == "maxΔ"
        levels = 0:50:400
        colormap = cgrad(:magma, length(levels); categorical=true)
        extendlow = nothing
        extendhigh = colormap[end]
        colormap = cgrad(colormap[1:end-1]; categorical=true)
    end

    for (icol, (basin_key, basin)) in enumerate(pairs(basins))

        x2D = zonalaverage(x3D, gridmetrics; mask = basin)

        local ax = Axis(fig[irow, icol],
            backgroundcolor=:lightgray,
            xgridvisible=false, ygridvisible=false,
            ylabel = "depth (m)")

        X = dropdims(maximum(lat, dims=1), dims=1)
        Y = zt
        Z = x2D
        co = contourf!(ax, X, Y, Z;
            levels,
            colormap,
            nan_color = :lightgray,
            extendlow,
            extendhigh,
        )
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

        hidexdecorations!(ax,
            label = irow < 2, ticklabels = irow < 2,
            ticks = irow < 2, grid = false)
        hideydecorations!(ax,
            label = icol > 1, ticklabels = icol > 1,
            ticks = icol > 1, grid = false)

        axs[irow, icol] = ax
    end
    cb = Colorbar(fig[irow, 4], contours[irow, 1];
        vertical = true, flipaxis = true,
        # ticks = (, cbarticklabelformat.(levels)),
        label = rich(str, " ", Γup, " (yr)"),
        )
    cb.height = Relative(1)

end



for (icol, (basin_str, xlims)) in enumerate(zip(basin_strs, basin_latlims))
    Label(fig[0, icol], basin_str, fontsize=20, tellwidth=false)
    colsize!(fig.layout, icol, Auto(xlims[2] - xlims[1]))
end

title = "$model $experiment $(time_window) reemergence time"
Label(fig[-1, 1:3], text = title, fontsize=20, tellwidth=false)

rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)

# save plot
outputfile = joinpath(outputdir, "reemergence_time_ZAVGs_vsmaxdiff.png")
@info "Saving reemergence time ZAVGs as image file:\n  $(outputfile)"
save(outputfile, fig)









fig = Figure(size = (1200, 1200), fontsize = 18)
axs = Array{Any,2}(undef, (2, 1))
contours = Array{Any,2}(undef, (2, 1))

for (irow, (x3D, str)) in enumerate(zip((Γinyr3D_mean, Γinyr3D_maxdiff), ("mean", "maxΔ")))

    if str == "mean" # mean
        levels = 0:100:1500
        colormap = :viridis
        colorrange = extrema(levels)
    elseif str == "maxΔ"
        levels = 0:50:400
        colorrange = extrema(levels)
        colormap = :magma
    end

    # Plot mean age at the seafloor level
    local title = "$model $experiment $(time_window) mean age at seafloor"

    ax = Axis(fig[irow,1]; title, xtickformat, ytickformat)


    # plot
    x2D = seafloorvalue(x3D, wet3D)
    plt = plotmap!(ax, x2D, gridmetrics; colorrange, colormap)

    Colorbar(fig[irow,2], plt, label=rich(str, " ", Γdown, " at seafloor (yr)"))
end

rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)

# save plot
outputfile = joinpath(outputdir, "mean_age_at_seafloor_vsmaxdiff.png")
@info "Saving ideal mean age at sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)















fig = Figure(size = (1200, 1200), fontsize = 18)
axs = Array{Any,2}(undef, (2, 1))
contours = Array{Any,2}(undef, (2, 1))

for (irow, (x3D, str)) in enumerate(zip((Γoutyr3D_mean, Γoutyr3D_maxdiff), ("mean", "maxΔ")))

    if str == "mean" # mean
        levels = 0:100:1500
        colormap = :viridis
        colorrange = extrema(levels)
    elseif str == "maxΔ"
        levels = 0:50:400
        colorrange = extrema(levels)
        colormap = :magma
    end

    # Plot mean age at the seafloor level
    local title = "$model $experiment $(time_window) reemergence time at seafloor"

    ax = Axis(fig[irow,1]; title, xtickformat, ytickformat)


    # plot
    x2D = seafloorvalue(x3D, wet3D)
    plt = plotmap!(ax, x2D, gridmetrics; colorrange, colormap)

    # poly!(ax, reverse.(OCEANS[OceanBasins.atlantic()].polygon))
    # poly!(ax, reverse.(OCEANS[OceanBasins.indian()].polygon))
    # poly!(ax, reverse.(OCEANS[OceanBasins.east_pacific()].polygon))
    # poly!(ax, reverse.(OCEANS[OceanBasins.west_pacific()].polygon))

    Colorbar(fig[irow,2], plt, label=rich(str, " ", Γup, " at seafloor (yr)"))
end

rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)

# save plot
outputfile = joinpath(outputdir, "reemergence_time_at_seafloor_vsmaxdiff.png")
@info "Saving mean reemergence time at sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)
