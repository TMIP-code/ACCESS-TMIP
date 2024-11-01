using Pkg
Pkg.activate(".")
Pkg.instantiate()

using OceanTransportMatrixBuilder
using NetCDF
using YAXArrays
using DataFrames
using DimensionalData
using SparseArrays
using LinearAlgebra
using Unitful
using Unitful: s, yr
using NaNStatistics
using Format
using CairoMakie
using GeoMakie
using Interpolations
using OceanBasins
using NaNStatistics

# Load functions for GM terms
include("GentMcWilliams.jl")
include("plotting_functions.jl")

model = "ACCESS-ESM1-5"
member = "r1i1p1f1"
CMIP_version = "CMIP5"
experiment = "historical"
time_window = "Jan1990-Dec1999"

# Gadi directory for input files
inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)"

# Load umo, vmo, mlotst, volcello, and areacello
mlotst_ds = open_dataset(joinpath(inputdir, "mlotst.nc"))
volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))

# Load variables in memory
mlotst = readcubedata(mlotst_ds.mlotst)
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

# Make indices
indices = makeindices(gridmetrics.v3D)

# unpack model grid
(; lon, lat, zt, v3D,) = gridmetrics
lev = zt
# unpack indices
(; wet3D, N) = indices

strs = ("resolved", "resolved_GM", "resolved_GM_submeso")

# Load ideal mean age and reemergence time
Γinyr3D = [open_dataset(joinpath(inputdir, "ideal_mean_age_$k.nc"))["age"] |> Array for k in strs]
Γoutyr3D = [open_dataset(joinpath(inputdir, "mean_reemergence_time_$k.nc"))["age"] |> Array for k in strs]
# TODO: rename "age" to something else for reexposure time (Γout)?

# Plot zonal averages

basin_keys = (:ATL, :PAC, :IND)
basin_strs = ("Atlantic", "Pacific", "Indian")
basin_functions = (isatlantic, ispacific, isindian)
basin_values = (reshape(f(lat[:], lon[:], OCEANS), size(lat)) for f in basin_functions)
basins = (; (basin_keys .=> basin_values)...)
basin_latlims_values = [clamp.((-5, +5) .+ extrema(lat[.!isnan.(v3D[:,:,1]) .& basin[:,:,1]]), -80, 80) for basin in basins]
basin_latlims = (; (basin_keys .=> basin_latlims_values)...)

levels = 0:250:2500
colormap = cgrad(:viridis, length(levels); categorical=true)
extendlow = nothing
extendhigh = colormap[end]
colormap = cgrad(colormap[1:end-1]; categorical=true)

Δlevels = -100:20:100
Δcolormap = cgrad(:RdBu, length(Δlevels) + 1; categorical=true, rev=true)
Δextendlow = Δcolormap[1]
Δextendhigh = Δcolormap[end]
Δcolormap = cgrad(Δcolormap[2:end-1]; categorical=true)

# Plot Γ↓ zonal averages

fig = Figure(size = (1200, 1200), fontsize = 18)
axs = Array{Any,2}(undef, (5, 3))
contours = Array{Any,2}(undef, (5, 3))
for (icol, (basin_key, basin)) in enumerate(pairs(basins))

    for (irow, (x3D, str)) in enumerate(zip(Γinyr3D, strs))

        x2D = zonalaverage(x3D, gridmetrics; mask = basin)

        local ax = Axis(fig[irow, icol],
            backgroundcolor=:lightgray,
            xgridvisible=true, ygridvisible=true,
            xgridcolor=(:black, 0.05), ygridcolor=(:black, 0.05),
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
            label = irow < 5, ticklabels = irow < 5,
            ticks = irow < 5, grid = false)
        hideydecorations!(ax,
            label = icol > 1, ticklabels = icol > 1,
            ticks = icol > 1, grid = false)


        axs[irow, icol] = ax
    end
end

Γdown = rich("Γ", superscript("↓"))
cb = Colorbar(fig[1:3, 4], contours[1, 1];
    vertical = true, flipaxis = true,
    # ticks = (, cbarticklabelformat.(levels)),
    label = rich(Γdown, " (yr)"),
    )
cb.height = Relative(0.666)

for (icol, (basin_key, basin)) in enumerate(pairs(basins))

    for (irow, (x3D, str)) in enumerate(zip(Γinyr3D[2:end], strs[2:end]))

        irow += 3

        x2D = zonalaverage(x3D - Γinyr3D[1], gridmetrics; mask = basin)

        local ax = Axis(fig[irow, icol],
            backgroundcolor=:lightgray,
            xgridvisible=true, ygridvisible=true,
            xgridcolor=(:black, 0.05), ygridcolor=(:black, 0.05),
            ylabel = "depth (m)")

        X = dropdims(maximum(lat, dims=1), dims=1)
        Y = zt
        Z = x2D
        co = contourf!(ax, X, Y, Z;
            levels = Δlevels,
            colormap = Δcolormap,
            nan_color = :lightgray,
            extendlow = Δextendlow,
            extendhigh = Δextendhigh,
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
            label = irow < 5, ticklabels = irow < 5,
            ticks = irow < 5, grid = false)
        hideydecorations!(ax,
            label = icol > 1, ticklabels = icol > 1,
            ticks = icol > 1, grid = false)


        axs[irow, icol] = ax
    end
end

cb = Colorbar(fig[4:5, 4], contours[4, 1];
    vertical = true, flipaxis = true,
    # ticks = (, cbarticklabelformat.(levels)),
    label = rich("Δ", Γdown, " (yr)"),
    )
cb.height = Relative(4/5)

for (icol, (basin_str, xlims)) in enumerate(zip(basin_strs, basin_latlims))
    Label(fig[0, icol], basin_str, fontsize=20, tellwidth=false)
    colsize!(fig.layout, icol, Auto(xlims[2] - xlims[1]))
end
Label(fig[1, 0], text = "resolved", fontsize=20, tellheight=false, rotation=π/2)
Label(fig[2, 0], text = "+GM", fontsize=20, tellheight=false, rotation=π/2)
Label(fig[3, 0], text = "+GM+submeso", fontsize=20, tellheight=false, rotation=π/2)
Label(fig[4, 0], text = "diff GM", fontsize=20, tellheight=false, rotation=π/2)
Label(fig[5, 0], text = "diff GM+submeso", fontsize=20, tellheight=false, rotation=π/2)

title = "$model $experiment $member $(time_window) ideal age"
Label(fig[-1, 1:3], text = title, fontsize=20, tellwidth=false)

# text = rich("Upstream sweeping time, ", ΓupΩ, ", for Ω = $(LONGTEXT[Ωz]) $(LONGTEXT[Ωbasin])")
# Label(f[-1, :]; text, fontsize=20)

# for (irun, run) in enumerate(runs)
#     prefix = irun > 1 ? "future " : ""
#     # Label(f[irun + 1, 0], text="future $(LONGTEXT[run])", fontsize=20, tellheight=false, rotation=π/2)
#     # Label(f[irun + 2, 0], text=LONGTEXT[run], fontsize=20, tellheight=false, rotation=π/2)
#     (irun == 1) && continue
#     Label(f[irun + 3, 0], text=LONGTEXT[run], fontsize=20, tellheight=false, rotation=π/2)
# end
# Label(f[5:6, -1], text="preindustrial-to-future change", fontsize=20, tellheight=false, rotation=π/2)


# # rowgap!(fig.layout, 1, 10)
# # rowgap!(fig.layout, 2, 5)
# # rowgap!(fig.layout, 5)
rowgap!(fig.layout, 10)
# rowgap!(fig.layout, 4, 15)
# # rowgap!(fig.layout, 5, 10)

colgap!(fig.layout, 10)
# save plot
outputfile = joinpath(inputdir, "ideal_age_ZAVGs_GMcomparison.png")
@info "Saving ideal age ZAVGs as image file:\n  $(outputfile)"
save(outputfile, fig)





Δlevels = -300:50:300
Δcolormap = cgrad(:RdBu, length(Δlevels) + 1; categorical=true, rev=true)
Δextendlow = Δcolormap[1]
Δextendhigh = Δcolormap[end]
Δcolormap = cgrad(Δcolormap[2:end-1]; categorical=true)


# Plot Γ↑ zonal averages

fig = Figure(size = (1200, 1200), fontsize = 18)
axs = Array{Any,2}(undef, (5, 3))
contours = Array{Any,2}(undef, (5, 3))

for (icol, (basin_key, basin)) in enumerate(pairs(basins))

    for (irow, (x3D, str)) in enumerate(zip(Γoutyr3D, strs))

        x2D = zonalaverage(x3D, gridmetrics; mask = basin)

        local ax = Axis(fig[irow, icol],
            backgroundcolor=:lightgray,
            xgridvisible=true, ygridvisible=true,
            xgridcolor=(:black, 0.05), ygridcolor=(:black, 0.05),
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
            label = irow < 5, ticklabels = irow < 5,
            ticks = irow < 5, grid = false)
        hideydecorations!(ax,
            label = icol > 1, ticklabels = icol > 1,
            ticks = icol > 1, grid = false)

        axs[irow, icol] = ax
    end
end
Γup = rich("Γ", superscript("↑"))
cb = Colorbar(fig[1:3, 4], contours[1, 1];
    vertical = true, flipaxis = true,
    # ticks = (, cbarticklabelformat.(levels)),
    label = rich(Γup, " (yr)"),
    )
cb.height = Relative(0.666)

for (icol, (basin_key, basin)) in enumerate(pairs(basins))

    for (irow, (x3D, str)) in enumerate(zip(Γoutyr3D[2:end], strs[2:end]))

        irow += 3

        x2D = zonalaverage(x3D - Γoutyr3D[1], gridmetrics; mask = basin)

        local ax = Axis(fig[irow, icol],
            backgroundcolor=:lightgray,
            xgridvisible=true, ygridvisible=true,
            xgridcolor=(:black, 0.05), ygridcolor=(:black, 0.05),
            ylabel = "depth (m)")

        X = dropdims(maximum(lat, dims=1), dims=1)
        Y = zt
        Z = x2D
        co = contourf!(ax, X, Y, Z;
            levels = Δlevels,
            colormap = Δcolormap,
            nan_color = :lightgray,
            extendlow = Δextendlow,
            extendhigh = Δextendhigh,
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
            label = irow < 5, ticklabels = irow < 5,
            ticks = irow < 5, grid = false)
        hideydecorations!(ax,
            label = icol > 1, ticklabels = icol > 1,
            ticks = icol > 1, grid = false)


        axs[irow, icol] = ax
    end
end

cb = Colorbar(fig[4:5, 4], contours[4, 1];
    vertical = true, flipaxis = true,
    # ticks = (, cbarticklabelformat.(levels)),
    label = rich("Δ", Γdown, " (yr)"),
    )
cb.height = Relative(4/5)

for (icol, (basin_str, xlims)) in enumerate(zip(basin_strs, basin_latlims))
    Label(fig[0, icol], basin_str, fontsize=20, tellwidth=false)
    colsize!(fig.layout, icol, Auto(xlims[2] - xlims[1]))
end

Label(fig[1, 0], text = "resolved", fontsize=20, tellheight=false, rotation=π/2)
Label(fig[2, 0], text = "+GM", fontsize=20, tellheight=false, rotation=π/2)
Label(fig[3, 0], text = "+GM+submeso", fontsize=20, tellheight=false, rotation=π/2)
Label(fig[4, 0], text = "diff GM", fontsize=20, tellheight=false, rotation=π/2)
Label(fig[5, 0], text = "diff GM+submeso", fontsize=20, tellheight=false, rotation=π/2)


title = "$model $experiment $member $(time_window) reemergence time"
Label(fig[-1, 1:3], text = title, fontsize=20, tellwidth=false)

rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)

# save plot
outputfile = joinpath(inputdir, "reemergence_time_ZAVGs_GMcomparison.png")
@info "Saving reemergence time ZAVGs as image file:\n  $(outputfile)"
save(outputfile, fig)





# Plot mean age at the seafloor level
title = "$model $experiment $member $(time_window) mean age at seafloor"
# plot options
colorrange = extrema(levels)
Δcolorrange = extrema(Δlevels)
# plot
fig = Figure(size = (800, 1400), fontsize = 18)
axs = Array{Any,2}(undef, (5, 1))
hms = Array{Any,2}(undef, (5, 1))
for (irow, (x3D, str)) in enumerate(zip(Γinyr3D, strs))

    icol = 1

    ax = Axis(fig[irow, icol]; xtickformat, ytickformat)
    Γinyrseafloor = seafloorvalue(x3D, wet3D)
    hms[irow, icol] = plotmap!(ax, Γinyrseafloor, gridmetrics; colorrange, colormap)


    hidexdecorations!(ax,
        label = irow < 5, ticklabels = irow < 5,
        ticks = irow < 5, grid = false)
    hideydecorations!(ax,
        label = icol > 1, ticklabels = icol > 1,
        ticks = icol > 1, grid = false)


end
cb = Colorbar(fig[1:3,2], hms[1,1], label=rich(Γdown, " at seafloor (yr)"))
cb.height = Relative(3/4)

for (irow, (x3D, str)) in enumerate(zip(Γinyr3D[2:end], strs[2:end]))

    icol = 1
    irow2 = irow + 3
    ax = Axis(fig[irow2, icol]; xtickformat, ytickformat)
    Γinyrseafloor = seafloorvalue(x3D - Γinyr3D[1], wet3D)
    hms[irow2, icol] = plotmap!(ax, Γinyrseafloor, gridmetrics;
        colorrange = Δcolorrange,
        colormap = Δcolormap)

    hidexdecorations!(ax,
        label = irow2 < 5, ticklabels = irow2 < 5,
        ticks = irow2 < 5, grid = false)
    hideydecorations!(ax,
        label = icol > 1, ticklabels = icol > 1,
        ticks = icol > 1, grid = false)

end
cb = Colorbar(fig[4:5,2], hms[4,1], label = rich("Δ", Γdown, " at seafloor (yr)"))
cb.height = Relative(5/6)

Label(fig[1, 0], text = "resolved", fontsize=20, tellheight=false, rotation=π/2)
Label(fig[2, 0], text = "+GM", fontsize=20, tellheight=false, rotation=π/2)
Label(fig[3, 0], text = "+GM+submeso", fontsize=20, tellheight=false, rotation=π/2)
Label(fig[4, 0], text = "diff GM", fontsize=20, tellheight=false, rotation=π/2)
Label(fig[5, 0], text = "diff GM+submeso", fontsize=20, tellheight=false, rotation=π/2)


title = "$model $experiment $member $(time_window) mean age at seafloor"
Label(fig[0, 1], text = title, fontsize=20, tellwidth=false)

rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)


# save plot
outputfile = joinpath(inputdir, "mean_age_at_seafloor_GMcomparison.png")
@info "Saving ideal mean age at sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)










# Plot reemergence time at the seafloor level
fig = Figure(size = (800, 1400), fontsize = 18)
axs = Array{Any,2}(undef, (5, 1))
hms = Array{Any,2}(undef, (5, 1))
for (irow, (x3D, str)) in enumerate(zip(Γoutyr3D, strs))

    icol = 1

    ax = Axis(fig[irow, icol]; xtickformat, ytickformat)
    Γinyrseafloor = seafloorvalue(x3D, wet3D)
    hms[irow, icol] = plotmap!(ax, Γinyrseafloor, gridmetrics; colorrange, colormap)


    hidexdecorations!(ax,
        label = irow < 5, ticklabels = irow < 5,
        ticks = irow < 5, grid = false)
    hideydecorations!(ax,
        label = icol > 1, ticklabels = icol > 1,
        ticks = icol > 1, grid = false)


end
cb = Colorbar(fig[1:3,2], hms[1,1], label=rich(Γup, " at seafloor (yr)"))
cb.height = Relative(3/4)

for (irow, (x3D, str)) in enumerate(zip(Γoutyr3D[2:end], strs[2:end]))

    icol = 1
    irow2 = irow + 3
    ax = Axis(fig[irow2, icol]; xtickformat, ytickformat)
    Γinyrseafloor = seafloorvalue(x3D - Γoutyr3D[1], wet3D)
    hms[irow2, icol] = plotmap!(ax, Γinyrseafloor, gridmetrics;
        colorrange = Δcolorrange,
        colormap = Δcolormap)

    hidexdecorations!(ax,
        label = irow2 < 5, ticklabels = irow2 < 5,
        ticks = irow2 < 5, grid = false)
    hideydecorations!(ax,
        label = icol > 1, ticklabels = icol > 1,
        ticks = icol > 1, grid = false)

end
cb = Colorbar(fig[4:5,2], hms[4,1], label = rich("Δ", Γup, " at seafloor (yr)"))
cb.height = Relative(5/6)

Label(fig[1, 0], text = "resolved", fontsize=20, tellheight=false, rotation=π/2)
Label(fig[2, 0], text = "+GM", fontsize=20, tellheight=false, rotation=π/2)
Label(fig[3, 0], text = "+GM+submeso", fontsize=20, tellheight=false, rotation=π/2)
Label(fig[4, 0], text = "diff GM", fontsize=20, tellheight=false, rotation=π/2)
Label(fig[5, 0], text = "diff GM+submeso", fontsize=20, tellheight=false, rotation=π/2)


title = "$model $experiment $member $(time_window) reemergence time at seafloor"
Label(fig[0, 1], text = title, fontsize=20, tellwidth=false)

rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)


# save plot
outputfile = joinpath(inputdir, "reemergence_time_at_seafloor_GMcomparison.png")
@info "Saving mean reemergence time at seafloor as image file:\n  $(outputfile)"
save(outputfile, fig)




