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
using CairoMakie
using GeoMakie
using Interpolations
using OceanBasins
using Statistics
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
inputdirfun(member) = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)"
outputdir = "/scratch/xv83/TMIP/data/$model/$experiment/all/$(time_window)"
mkpath(outputdir)

# find all members for which the inputdir contains umo.nc, vmo.nc, mlotst.nc, volcello.nc, and areacello.nc
transportmatrix_filename = "transportmatrix_resolved_GM_submeso.jld2"
circulation_timescales_filename = "circulation_timescales_resolved_GM_submeso.nc"
requiredfiles = [circulation_timescales_filename, transportmatrix_filename]
hasrequireddata(member, file_name) = isfile(joinpath(inputdirfun(member), file_name))
hasrequireddata(member) = all(file_name -> hasrequireddata(member, file_name), requiredfiles)


# sort members by r, i, p[, f]
member_regex = CMIP_version == "CMIP6" ? r"r(\d+)i(\d+)p(\d+)f(\d+)" : r"r(\d+)i(\d+)p(\d+)"
members = [m for m in readdir("/scratch/xv83/TMIP/data/$model/$experiment") if !isnothing(match(member_regex, m))]
parse_member(member) = parse.(Int, match(member_regex, member).captures)
members = sort(members, by = x -> parse_member(x))
dataavailability = DataFrame(
    :member => members,
    :has_it_all => hasrequireddata.(members),
    [Symbol(f) => [hasrequireddata(m, f) for m in members] for f in requiredfiles]...,
)
show(dataavailability; allrows = true)
println()

@info "Merging valid members into single data cube"
@show valid_members = members[dataavailability.has_it_all]
cubes = [open_dataset(joinpath(inputdirfun(m), circulation_timescales_filename)).ideal_age for m in valid_members]
Γinyr3Ds = readcubedata(concatenatecubes(cubes, Dim{:member}(valid_members)))
cubes = [open_dataset(joinpath(inputdirfun(m), circulation_timescales_filename)).reemergence_time for m in valid_members]
Γoutyr3Ds = readcubedata(concatenatecubes(cubes, Dim{:member}(valid_members)))

@info "Computing mean and std"
Γinyr3D_mean = dropdims(mean(Γinyr3Ds; dims = :member), dims = :member)
Γinyr3D_std = dropdims(stdm(Γinyr3Ds, Γinyr3D_mean; dims = :member), dims = :member)
Γoutyr3D_mean = dropdims(mean(Γoutyr3Ds; dims = :member), dims = :member)
Γoutyr3D_std = dropdims(stdm(Γoutyr3Ds, Γoutyr3D_mean; dims = :member), dims = :member)


Γdown = rich("Γ", superscript("↓"))
Γup = rich("Γ", superscript("↑"))

member = first(valid_members)

inputdir = inputdirfun(member)

transportmatrix_filepath = joinpath(inputdir, transportmatrix_filename)

gridmetrics, indices = load(transportmatrix_filepath, "gridmetrics", "indices")

# unpack model grid
(; lon, lat, zt, v3D) = gridmetrics
lev = zt
# unpack indices
(; wet3D, N) = indices


# Plot zonal averages

basin_keys = (:ATL, :PAC, :IND)
basin_strs = ("Atlantic", "Pacific", "Indian")
basin_functions = (isatlantic, ispacific, isindian)
basin_values = (reshape(f(lat[:], lon[:], OCEANS), size(lat)) for f in basin_functions)
basins = (; (basin_keys .=> basin_values)...)
basin_latlims_values = [clamp.((-5, +5) .+ extrema(lat[.!isnan.(v3D[:, :, 1]) .& basin[:, :, 1]]), -80, 80) for basin in basins]
basin_latlims = (; (basin_keys .=> basin_latlims_values)...)


# Plot Γ↓ zonal averages

fig = Figure(size = (1200, 600), fontsize = 18)
axs = Array{Any, 2}(undef, (2, 3))
contours = Array{Any, 2}(undef, (2, 3))

for (irow, (x3D, str)) in enumerate(zip((Γinyr3D_mean, Γinyr3D_std), ("mean", "std")))

    if str == "mean" # mean
        levels = 0:100:2500
        colormap = cgrad(:viridis, length(levels); categorical = true)
        extendlow = nothing
        extendhigh = colormap[end]
        colormap = cgrad(colormap[1:(end - 1)]; categorical = true)
    elseif str == "std"
        levels = 0:10:200
        colormap = cgrad(:viridis, length(levels); categorical = true)
        extendlow = nothing
        extendhigh = colormap[end]
        colormap = cgrad(colormap[1:(end - 1)]; categorical = true)
    end

    for (icol, (basin_key, basin)) in enumerate(pairs(basins))

        x2D = zonalaverage(x3D, gridmetrics; mask = basin)

        local ax = Axis(
            fig[irow, icol],
            backgroundcolor = :lightgray,
            xgridvisible = false, ygridvisible = false,
            ylabel = "depth (m)"
        )

        X = dropdims(maximum(lat, dims = 1), dims = 1)
        Y = zt
        Z = x2D

        co = contourf!(
            ax, X, Y, Z;
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

        hidexdecorations!(
            ax,
            label = irow < 2, ticklabels = irow < 2,
            ticks = irow < 2, grid = false
        )
        hideydecorations!(
            ax,
            label = icol > 1, ticklabels = icol > 1,
            ticks = icol > 1, grid = false
        )


        axs[irow, icol] = ax

    end

    cb = Colorbar(
        fig[irow, 4], contours[irow, 1];
        vertical = true, flipaxis = true,
        # ticks = (, cbarticklabelformat.(levels)),
        label = rich(str, " ", Γdown, " (yr)"),
    )
    cb.height = Relative(1)
end


for (icol, (basin_str, xlims)) in enumerate(zip(basin_strs, basin_latlims))
    Label(fig[0, icol], basin_str, fontsize = 20, tellwidth = false)
    colsize!(fig.layout, icol, Auto(xlims[2] - xlims[1]))
end

title = "$model $experiment all members $(time_window) ideal age"
Label(fig[-1, 1:3], text = title, fontsize = 20, tellwidth = false)

rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)

# save plot
outputfile = joinpath(outputdir, "ideal_age_ZAVGs_v3.png")
@info "Saving ideal age ZAVGs as image file:\n  $(outputfile)"
save(outputfile, fig)


# Plot Γ↑ zonal averages


fig = Figure(size = (1200, 600), fontsize = 18)
axs = Array{Any, 2}(undef, (2, 3))
contours = Array{Any, 2}(undef, (2, 3))

for (irow, (x3D, str)) in enumerate(zip((Γoutyr3D_mean, Γoutyr3D_std), ("mean", "std")))

    if str == "mean" # mean
        levels = 0:100:2500
        colormap = cgrad(:viridis, length(levels); categorical = true)
        extendlow = nothing
        extendhigh = colormap[end]
        colormap = cgrad(colormap[1:(end - 1)]; categorical = true)
    elseif str == "std"
        levels = 0:10:200
        colormap = cgrad(:viridis, length(levels); categorical = true)
        extendlow = nothing
        extendhigh = colormap[end]
        colormap = cgrad(colormap[1:(end - 1)]; categorical = true)
    end

    for (icol, (basin_key, basin)) in enumerate(pairs(basins))

        x2D = zonalaverage(x3D, gridmetrics; mask = basin)

        local ax = Axis(
            fig[irow, icol],
            backgroundcolor = :lightgray,
            xgridvisible = false, ygridvisible = false,
            ylabel = "depth (m)"
        )

        X = dropdims(maximum(lat, dims = 1), dims = 1)
        Y = zt
        Z = x2D
        co = contourf!(
            ax, X, Y, Z;
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

        hidexdecorations!(
            ax,
            label = irow < 2, ticklabels = irow < 2,
            ticks = irow < 2, grid = false
        )
        hideydecorations!(
            ax,
            label = icol > 1, ticklabels = icol > 1,
            ticks = icol > 1, grid = false
        )

        axs[irow, icol] = ax
    end
    cb = Colorbar(
        fig[irow, 4], contours[irow, 1];
        vertical = true, flipaxis = true,
        # ticks = (, cbarticklabelformat.(levels)),
        label = rich(str, " ", Γup, " (yr)"),
    )
    cb.height = Relative(1)

end


for (icol, (basin_str, xlims)) in enumerate(zip(basin_strs, basin_latlims))
    Label(fig[0, icol], basin_str, fontsize = 20, tellwidth = false)
    colsize!(fig.layout, icol, Auto(xlims[2] - xlims[1]))
end

title = "$model $experiment all members $(time_window) reemergence time"
Label(fig[-1, 1:3], text = title, fontsize = 20, tellwidth = false)

rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)

# save plot
outputfile = joinpath(outputdir, "reemergence_time_ZAVGs_v3.png")
@info "Saving reemergence time ZAVGs as image file:\n  $(outputfile)"
save(outputfile, fig)


fig = Figure(size = (1200, 1200), fontsize = 18)
axs = Array{Any, 2}(undef, (2, 1))
contours = Array{Any, 2}(undef, (2, 1))

for (irow, (x3D, str)) in enumerate(zip((Γinyr3D_mean, Γinyr3D_std), ("mean", "std")))

    if str == "mean" # mean
        levels = 0:100:2500
        colormap = :viridis
        colorrange = extrema(levels)
    elseif str == "std"
        levels = 0:10:200
        colorrange = extrema(levels)
        colormap = :viridis
    end

    # Plot mean age at the seafloor level
    local title = "$model $experiment $member $(time_window) mean age at seafloor"

    ax = Axis(fig[irow, 1]; title, xtickformat, ytickformat)


    # plot
    x2D = seafloorvalue(x3D, wet3D)
    plt = plotmap!(ax, x2D, gridmetrics; colorrange, colormap)

    Colorbar(fig[irow, 2], plt, label = rich(str, " ", Γdown, " at seafloor (yr)"))
end

rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)

# save plot
outputfile = joinpath(outputdir, "mean_age_at_seafloor_v3.png")
@info "Saving ideal mean age at sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)


fig = Figure(size = (1200, 1200), fontsize = 18)
axs = Array{Any, 2}(undef, (2, 1))
contours = Array{Any, 2}(undef, (2, 1))

for (irow, (x3D, str)) in enumerate(zip((Γoutyr3D_mean, Γoutyr3D_std), ("mean", "std")))

    if str == "mean" # mean
        levels = 0:100:2500
        colormap = :viridis
        colorrange = extrema(levels)
    elseif str == "std"
        levels = 0:10:200
        colorrange = extrema(levels)
        colormap = :viridis
    end

    # Plot mean age at the seafloor level
    local title = "$model $experiment $member $(time_window) reemergence time at seafloor"

    ax = Axis(fig[irow, 1]; title, xtickformat, ytickformat)


    # plot
    x2D = seafloorvalue(x3D, wet3D)
    plt = plotmap!(ax, x2D, gridmetrics; colorrange, colormap)

    # poly!(ax, reverse.(OCEANS[OceanBasins.atlantic()].polygon))
    # poly!(ax, reverse.(OCEANS[OceanBasins.indian()].polygon))
    # poly!(ax, reverse.(OCEANS[OceanBasins.east_pacific()].polygon))
    # poly!(ax, reverse.(OCEANS[OceanBasins.west_pacific()].polygon))

    Colorbar(fig[irow, 2], plt, label = rich(str, " ", Γup, " at seafloor (yr)"))
end

rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)

# save plot
outputfile = joinpath(outputdir, "reemergence_time_at_seafloor_v3.png")
@info "Saving mean reemergence time at sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)
