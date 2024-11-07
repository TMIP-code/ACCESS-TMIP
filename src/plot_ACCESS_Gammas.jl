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
using NaNStatistics
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

# find all members for which the inputdir contains umo.nc, vmo.nc, mlotst.nc, volcello.nc, and areacello.nc
transportmatrix_filename = "transportmatrix_resolved_GM_submeso.jld2"
circulation_timescales_filename = "circulation_timescales_resolved_GM_submeso.nc"
requiredfiles = [circulation_timescales_filename, transportmatrix_filename]
hasrequireddata(member, file_name) = isfile(joinpath(inputdirfun(member), file_name))
hasrequireddata(member) = all(file_name -> hasrequireddata(member, file_name), requiredfiles)
members = readdir("/scratch/xv83/TMIP/data/$model/$experiment")

# sort members by r, i, p[, f]

member_regex = CMIP_version == "CMIP6" ? r"r(\d+)i(\d+)p(\d+)f(\d+)" : r"r(\d+)i(\d+)p(\d+)"
parse_member(member) = parse.(Int, match(member_regex, member).captures)
members = sort(members, by = x -> parse_member(x))
dataavailability = DataFrame(
    :member => members,
    :has_it_all => hasrequireddata.(members),
    [Symbol(f) => [hasrequireddata(m, f) for m in members] for f in requiredfiles]...,
)
show(dataavailability; allrows = true)
println()

Γdown = rich("Γ", superscript("↓"))
Γup = rich("Γ", superscript("↑"))

for member in members[dataavailability.has_it_all][4:end]
# for member in [last(members)]

    inputdir = inputdirfun(member)

    transportmatrix_filepath = joinpath(inputdir, transportmatrix_filename)

    gridmetrics, indices = load(transportmatrix_filepath, "gridmetrics", "indices")

    # unpack model grid
    (; lon, lat, zt, v3D,) = gridmetrics
    lev = zt
    # unpack indices
    (; wet3D, N) = indices

    # Load ideal mean age and reemergence time
    circulation_timescales_filepath = joinpath(inputdir, circulation_timescales_filename)
    circulation_timescales_ds = open_dataset(circulation_timescales_filepath)
    Γinyr3D = circulation_timescales_ds.ideal_age |> Array
    Γoutyr3D = circulation_timescales_ds.reemergence_time |> Array


    # Plot zonal averages

    basin_keys = (:ATL, :PAC, :IND)
    basin_strs = ("Atlantic", "Pacific", "Indian")
    basin_functions = (isatlantic, ispacific, isindian)
    basin_values = (reshape(f(lat[:], lon[:], OCEANS), size(lat)) for f in basin_functions)
    basins = (; (basin_keys .=> basin_values)...)
    basin_latlims_values = [clamp.((-5, +5) .+ extrema(lat[.!isnan.(v3D[:,:,1]) .& basin[:,:,1]]), -80, 80) for basin in basins]
    basin_latlims = (; (basin_keys .=> basin_latlims_values)...)

    levels = 0:100:1500
    colormap = cgrad(:viridis, length(levels); categorical=true)
    extendlow = nothing
    extendhigh = colormap[end]
    colormap = cgrad(colormap[1:end-1]; categorical=true)

    # Plot Γ↓ zonal averages

    fig = Figure(size = (1200, 300), fontsize = 18)
    axs = Array{Any,2}(undef, (1, 3))
    contours = Array{Any,2}(undef, (1, 3))
    for (icol, (basin_key, basin)) in enumerate(pairs(basins))

        irow = 1

        x2D = zonalaverage(Γinyr3D, gridmetrics; mask = basin)

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
            label = irow < 1, ticklabels = irow < 1,
            ticks = irow < 1, grid = false)
        hideydecorations!(ax,
            label = icol > 1, ticklabels = icol > 1,
            ticks = icol > 1, grid = false)


        axs[irow, icol] = ax
    end

    cb = Colorbar(fig[1, 4], contours[1, 1];
        vertical = true, flipaxis = true,
        # ticks = (, cbarticklabelformat.(levels)),
        label = rich(Γdown, " (yr)"),
        )
    cb.height = Relative(1)


    for (icol, (basin_str, xlims)) in enumerate(zip(basin_strs, basin_latlims))
        Label(fig[0, icol], basin_str, fontsize=20, tellwidth=false)
        colsize!(fig.layout, icol, Auto(xlims[2] - xlims[1]))
    end

    title = "$model $experiment $member $(time_window) ideal age"
    Label(fig[-1, 1:3], text = title, fontsize=20, tellwidth=false)

    rowgap!(fig.layout, 10)
    colgap!(fig.layout, 10)

    # save plot
    outputfile = joinpath(inputdir, "ideal_age_ZAVGs_v3.png")
    @info "Saving ideal age ZAVGs as image file:\n  $(outputfile)"
    save(outputfile, fig)







    # Plot Γ↑ zonal averages

    fig = Figure(size = (1200, 300), fontsize = 18)
    axs = Array{Any,2}(undef, (1, 3))
    contours = Array{Any,2}(undef, (1, 3))
    for (icol, (basin_key, basin)) in enumerate(pairs(basins))

        x2D = zonalaverage(Γoutyr3D, gridmetrics; mask = basin)
        irow = 1

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
            label = irow < 1, ticklabels = irow < 1,
            ticks = irow < 1, grid = false)
        hideydecorations!(ax,
            label = icol > 1, ticklabels = icol > 1,
            ticks = icol > 1, grid = false)

        axs[irow, icol] = ax
    end

    cb = Colorbar(fig[1, 4], contours[1, 1];
        vertical = true, flipaxis = true,
        # ticks = (, cbarticklabelformat.(levels)),
        label = rich(Γup, " (yr)"),
        )
    cb.height = Relative(1)


    for (icol, (basin_str, xlims)) in enumerate(zip(basin_strs, basin_latlims))
        Label(fig[0, icol], basin_str, fontsize=20, tellwidth=false)
        colsize!(fig.layout, icol, Auto(xlims[2] - xlims[1]))
    end

    title = "$model $experiment $member $(time_window) reemergence time"
    Label(fig[-1, 1:3], text = title, fontsize=20, tellwidth=false)

    rowgap!(fig.layout, 10)
    colgap!(fig.layout, 10)

    # save plot
    outputfile = joinpath(inputdir, "reemergence_time_ZAVGs_v3.png")
    @info "Saving reemergence time ZAVGs as image file:\n  $(outputfile)"
    save(outputfile, fig)


    # Plot mean age at the seafloor level
    Γinyrseafloor = seafloorvalue(Γinyr3D, wet3D)
    title = "$model $experiment $member $(time_window) mean age at seafloor"
    # plot options
    colorrange = (0, 1500)
    colormap = :viridis
    # plot
    fig = Figure(size = (1200, 600), fontsize = 18)
    ax = Axis(fig[1,1]; title, xtickformat, ytickformat)
    plt = plotmap!(ax, Γinyrseafloor, gridmetrics; colorrange, colormap)
    Colorbar(fig[1,2], plt, label=rich(Γup, " at seafloor (yr)"))
    # save plot
    outputfile = joinpath(inputdir, "mean_age_at_seafloor_v3.png")
    @info "Saving ideal mean age at sea floor as image file:\n  $(outputfile)"
    save(outputfile, fig)



    # Plot reemergence time at the seafloor level
    Γoutyrseafloor = seafloorvalue(Γoutyr3D, wet3D)
    title = "$model $experiment $member $(time_window) reemergence time at seafloor"
    # plot options
    colorrange = (0, 1500)
    colormap = :viridis
    # plot
    fig = Figure(size = (1200, 600), fontsize = 18)
    ax = Axis(fig[1,1]; title, xtickformat, ytickformat)
    plt = plotmap!(ax, Γoutyrseafloor, gridmetrics; colorrange, colormap)
    Colorbar(fig[1,2], plt, label=rich(Γup, " at seafloor (yr)"))
    # save plot
    outputfile = joinpath(inputdir, "reemergence_time_at_seafloor_v3.png")
    @info "Saving mean reemergence time at seafloor as image file:\n  $(outputfile)"
    save(outputfile, fig)


end


