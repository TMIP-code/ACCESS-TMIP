# using Pkg
# Pkg.activate(".")
# Pkg.instantiate()

# using OceanTransportMatrixBuilder
# using NetCDF
# using YAXArrays
# using DataFrames
# using DimensionalData
# # using SparseArrays
# # using LinearAlgebra
# using Unitful
# using Unitful: s, yr
# # Must run twice v1.11.1 https://github.com/JuliaLang/julia/issues/56216
# try
#     using CairoMakie
# catch
#     using CairoMakie
# end
# using GeoMakie
# using Interpolations
# using OceanBasins
# using Statistics
# using StatsBase
# using FileIO
# using Dates
# using NaNStatistics

# include("plotting_functions.jl")

# model = "ACCESS-ESM1-5"
# # model = "ACCESS-CM2"
# # model = "ACCESS1-3"

# # CMIP_version = "CMIP5"
# CMIP_version = "CMIP6"

# experiment = "historical"
# # experiment = "piControl"

# time_window = "Jan1990-Dec1999"
# # time_window = "Jan1071-Dec1100" # <- last 30 years of ACCESS-ESM1-5 piControl
# # time_window = "Jan1420-Dec1449" # <- last 30 years of ACCESS-CM2 piControl

# lumpby = "month"
# months = 1:12
# Nmonths = length(months)

# # Gadi directory for input files
# # Gadi directory for input files
# inputdirfun(member) = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)"
# cycloinputdirfun(member) = joinpath(inputdirfun(member), "cyclo$lumpby")

# # find all members for which the inputdir contains umo.nc, vmo.nc, mlotst.nc, volcello.nc, and areacello.nc
# circulation_timescales_filename = "ideal_mean_age.nc"
# requiredfiles = [circulation_timescales_filename]
# hasrequireddata(member, file_name) = isfile(joinpath(cycloinputdirfun(member), file_name))
# hasrequireddata(member) = all(file_name -> hasrequireddata(member, file_name), requiredfiles)


# # sort members by r, i, p[, f]
# member_regex = CMIP_version == "CMIP6" ? r"r(\d+)i(\d+)p(\d+)f(\d+)" : r"r(\d+)i(\d+)p(\d+)"
# members = [m for m in readdir("/scratch/xv83/TMIP/data/$model/$experiment") if !isnothing(match(member_regex, m))]
# parse_member(member) = parse.(Int, match(member_regex, member).captures)
# members = sort(members, by = x -> parse_member(x))
# dataavailability = DataFrame(
#     :member => members,
#     :has_it_all => hasrequireddata.(members),
#     [Symbol(f) => [hasrequireddata(m, f) for m in members] for f in requiredfiles]...,
# )
# show(dataavailability; allrows = true)
# println()

# @show valid_members = members[dataavailability.has_it_all]

# contouroptions = let
#     levels = 0:100:2500
#     colormap = cgrad(:viridis, length(levels); categorical=true)
#     extendlow = nothing
#     extendhigh = colormap[end]
#     colormap = cgrad(colormap[1:end-1]; categorical=true)
#     nan_color = :lightgray
#     (; levels, colormap, extendlow, extendhigh, nan_color)
# end

# @show axisoptions = (
#     backgroundcolor = :lightgray,
#     xgridvisible = true,
#     ygridvisible = true,
#     ylabel = "depth (m)"
# )


# member = first(valid_members)

# inputdir = inputdirfun(member)
# cycloinputdir = cycloinputdirfun(member)
# outputdir = cycloinputdir

# # Load age into memory
# ds = open_dataset(joinpath(cycloinputdir, circulation_timescales_filename))
# Γinyr4D = readcubedata(ds.age)

# # Load areacello and volcello for grid geometry
# volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
# areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))

# # Load fixed variables in memory
# areacello = readcubedata(areacello_ds.areacello)
# volcello = readcubedata(volcello_ds.volcello)
# lon = readcubedata(volcello_ds.lon)
# lat = readcubedata(volcello_ds.lat)
# lat2 = dropdims(maximum(lat, dims=1), dims=1) |> Array # <- for plotting ZAVG (inexact)
# lev = volcello_ds.lev
# # Identify the vertices keys (vary across CMIPs / models)
# volcello_keys = propertynames(volcello_ds)
# lon_vertices_key = volcello_keys[findfirst(x -> occursin("lon", x) & occursin("vert", x), string.(volcello_keys))]
# lat_vertices_key = volcello_keys[findfirst(x -> occursin("lat", x) & occursin("vert", x), string.(volcello_keys))]
# lon_vertices = readcubedata(getproperty(volcello_ds, lon_vertices_key))
# lat_vertices = readcubedata(getproperty(volcello_ds, lat_vertices_key))

# # Make makegridmetrics
# gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
# (; lon, lat, zt, v3D,) = gridmetrics
# lev = zt
# # indices
# indices = makeindices(v3D)
# (; wet3D, N) = indices

# Γdown = rich("Γ", superscript("↓"))
# Γup = rich("Γ", superscript("↑"))


# # Plot zonal averages

# basin_keys = (:ATL, :PAC, :IND)
# basin_strs = ("Atlantic", "Pacific", "Indian")
# basin_functions = (isatlantic, ispacific, isindian)
# basin_values = (reshape(f(lat[:], lon[:], OCEANS), size(lat)) for f in basin_functions)
# basins = (; (basin_keys .=> basin_values)...)
# basin_latlims_values = [clamp.((-5, +5) .+ extrema(lat[.!isnan.(v3D[:,:,1]) .& basin[:,:,1]]), -80, 80) for basin in basins]
# basin_latlims = (; (basin_keys .=> basin_latlims_values)...)

# # Plot Γ↓ zonal averages
# fig = Figure(size = (3*400, 13*280), fontsize = 18)
# axs = Array{Any,2}(undef, (Nmonths + 1, 3))
# contours = Array{Any,2}(undef, (Nmonths + 1, 3))


# for (irow, month) in enumerate(months)

#     Γinyr3D = Γinyr4D[Ti=At(month)]

#     for (icol, (basin_key, basin)) in enumerate(pairs(basins))

#         Γinyr2D = zonalaverage(Γinyr3D, gridmetrics; mask = basin)

#         local ax = Axis(fig[irow, icol]; axisoptions...)

#         co = contourf!(ax, lat2, zt, Γinyr2D; contouroptions...)

#         translate!(co, 0, 0, -100)
#         contours[irow, icol] = co

#         xlim = basin_latlims[basin_key]
#         # basin2 = LONGTEXT[basin]

#         ax.yticks = (ztick, zticklabel)
#         xticks = -90:30:90
#         ax.xticks = (xticks, latticklabel.(xticks))
#         ylims!(ax, zlim)
#         # xlims!(ax, (-90, 90))
#         xlims!(ax, xlim)

#         myhidexdecorations!(ax, irow < Nmonths + 1)
#         myhideydecorations!(ax, icol > 1)

#         axs[irow, icol] = ax

#     end

#     Label(fig[irow, 0], monthabbr(month), fontsize=20, tellheight=false, rotation=π/2)


# end
# # The plot the annual mean
# irow = Nmonths + 1

# for (icol, (basin_key, basin)) in enumerate(pairs(basins))

#     Γinyr3D = average(Γinyr4D, gridmetrics; dims=4)

#     Γinyr2D = zonalaverage(Γinyr3D, gridmetrics; mask = basin)

#     local ax = Axis(fig[irow, icol]; axisoptions...)


#     co = contourf!(ax, lat2, zt, Γinyr2D; contouroptions...)
#     translate!(co, 0, 0, -100)
#     contours[irow, icol] = co

#     xlim = basin_latlims[basin_key]
#     # basin2 = LONGTEXT[basin]

#     ax.yticks = (ztick, zticklabel)
#     xticks = -90:30:90
#     ax.xticks = (xticks, latticklabel.(xticks))
#     ylims!(ax, zlim)
#     # xlims!(ax, (-90, 90))
#     xlims!(ax, xlim)

#     myhidexdecorations!(ax, irow < Nmonths + 1)
#     myhideydecorations!(ax, icol > 1)

#     Label(fig[irow, 0], "mean", fontsize=20, tellheight=false, rotation=π/2)


#     axs[irow, icol] = ax

# end

# cb = Colorbar(fig[1:13, 4], contours[irow, 1];
#     vertical = true, flipaxis = true,
#     # ticks = (, cbarticklabelformat.(levels)),
#     label = rich(Γdown, " (yr)"),
#     )
# cb.height = Relative(1/3)


# for (icol, (basin_str, xlims)) in enumerate(zip(basin_strs, basin_latlims))
#     Label(fig[0, icol], basin_str, fontsize=20, tellwidth=false)
#     colsize!(fig.layout, icol, Auto(xlims[2] - xlims[1]))
# end

# title = "$model $experiment all members $(time_window) ideal age"
# Label(fig[-1, 1:3], text = title, fontsize=20, tellwidth=false)

# rowgap!(fig.layout, 10)
# colgap!(fig.layout, 10)

# # save plot
# outputfile = joinpath(outputdir, "cyclo_stationary_age_ZAVGs.png")
# @info "Saving ideal age ZAVGs as image file:\n  $(outputfile)"
# save(outputfile, fig)

# # Check seasonality
# foo = average(Γinyr4D, gridmetrics; dims = (1, 2, 3))
# fig = Figure()
# ax = Axis(fig[1,1];
#     xlabel = "Month",
#     ylabel = rich("Global mean ", Γdown),
#     xticks = (months, monthabbr.(months))
# )
# plt = lines!(ax, foo)
# outputfile = joinpath(outputdir, "cyclo_stationary_age_global_mean.png")
# @info "Saving ideal age ZAVGs as image file:\n  $(outputfile)"
# save(outputfile, fig)


# # Plot seasonal variation
# contouroptionsvariability = let
#     levels = [1, 2, 5, 10, 20, 50]
#     colormap = cgrad([:white, Makie.wong_colors()[6]], length(levels) + 1; categorical=true)
#     extendhigh = :auto
#     extendlow = :auto
#     # colormap = cgrad(colormap[2:end-1]; categorical=true)
#     nan_color = :lightgray
#     (; levels, colormap, extendlow, extendhigh, nan_color, colorscale=log10)
# end
# Γinyr3Dcyclomean = average(Γinyr4D, gridmetrics; dims = 4)
# Γinyr3Dmin = dropdims(nanminimum(Γinyr4D; dims = 4), dims = 4)
# Γinyr3Dmax = dropdims(nanmaximum(Γinyr4D; dims = 4), dims = 4)
# Γinyr3Ddiff = Γinyr3Dmax - Γinyr3Dmin
# fig = Figure(size = (3*400, 1*280), fontsize = 18)
# for (icol, (basin_key, basin)) in enumerate(pairs(basins))
#     Γinyr2Ddiff = zonalaverage(Γinyr3Ddiff, gridmetrics; mask = basin)
#     Γinyr2Dmean = zonalaverage(Γinyr3Dcyclomean, gridmetrics; mask = basin)
#     Γinyr2Dpercentdiff = 100Γinyr2Ddiff ./ Γinyr2Dmean
#     local irow = 1
#     local ax = Axis(fig[irow, icol]; axisoptions...)
#     co = contourf!(ax, lat2, zt, Γinyr2Dpercentdiff; contouroptionsvariability...)
#     translate!(co, 0, 0, -100)
#     contours[irow, icol] = co

#     xlim = basin_latlims[basin_key]
#     # basin2 = LONGTEXT[basin]

#     ax.yticks = (ztick, zticklabel)
#     xticks = -90:30:90
#     ax.xticks = (xticks, latticklabel.(xticks))
#     ylims!(ax, zlim)
#     # xlims!(ax, (-90, 90))
#     xlims!(ax, xlim)

#     myhidexdecorations!(ax, irow < 1)
#     myhideydecorations!(ax, icol > 1)

#     Label(fig[irow, 0], "seasonal variability", fontsize=20, tellheight=false, rotation=π/2)
# end
# for (icol, (basin_str, xlims)) in enumerate(zip(basin_strs, basin_latlims))
#     Label(fig[0, icol], basin_str, fontsize=20, tellwidth=false)
#     colsize!(fig.layout, icol, Auto(xlims[2] - xlims[1]))
# end
# cb = Colorbar(fig[1, 4], contours[1,2];
#     vertical = true, flipaxis = true,
#     ticks = contouroptionsvariability.levels,
#     label = rich("variability (%)"),
#     # scale = Ident
# )
# cb.height = Relative(5/6)
# outputfile = joinpath(outputdir, "cyclo_stationary_age_global_variability.png")
# @info "Saving ideal age ZAVGs as image file:\n  $(outputfile)"
# save(outputfile, fig)


# # Compare age built from cyclo(T) vs T(mean) vs mean(T)
# strs = ("cyclo", "mean(T)", "T(mean)")
# # Load ideal mean age and reemergence time
# Γinyr3DTmean = open_dataset(joinpath(inputdir, "ideal_mean_age_resolved_GM_submeso.nc"))["age"] |> Array
# Γinyr3DmeanT = open_dataset(joinpath(inputdir, "ideal_mean_age_monthlymatrices.nc"))["age"] |> Array
# Γ3Ds = (Γinyr3Dcyclomean, Γinyr3DmeanT, Γinyr3DTmean)

# contouroptionscomparison = let
#     levels = -1000:100:1000
#     colormap = cgrad(:balance)
#     extendhigh = :auto
#     extendlow = :auto
#     # colormap = cgrad(colormap[2:end-1]; categorical=true)
#     nan_color = :lightgray
#     (; levels, colormap, extendlow, extendhigh, nan_color)
# end


Nrows = length(strs)
fig = Figure(size = (1200, 250 * (2Nrows - 1)), fontsize = 18)
axs = Array{Any, 2}(undef, (2Nrows - 1, 3))
contours = Array{Any, 2}(undef, (2Nrows - 1, 3))
for (icol, (basin_key, basin)) in enumerate(pairs(basins))

    for (irow, (x3D, str)) in enumerate(zip(Γ3Ds, strs))

        x2D = zonalaverage(x3D, gridmetrics; mask = basin)

        local ax = Axis(fig[irow, icol]; axisoptions...)

        co = contourf!(ax, lat2, zt, x2D; contouroptions...)

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

        myhidexdecorations!(ax, irow < 2Nrows - 1)
        myhideydecorations!(ax, icol > 1)

        axs[irow, icol] = ax
    end
end

Γdown = rich("Γ", superscript("↓"))
cb = Colorbar(
    fig[1:Nrows, 4], contours[1, 1];
    vertical = true, flipaxis = true,
    # ticks = (, cbarticklabelformat.(levels)),
    label = rich(Γdown, " (yr)"),
)
cb.height = Relative(Nrows / (Nrows + 1))


for (icol, (basin_key, basin)) in enumerate(pairs(basins))

    for (irow, (x3D, str)) in enumerate(zip(Γ3Ds[2:end], strs[2:end]))

        irow += Nrows

        x2D = zonalaverage(x3D - Γinyr3Dcyclomean, gridmetrics; mask = basin)

        local ax = Axis(fig[irow, icol]; axisoptions...)

        co = contourf!(ax, lat2, zt, x2D; contouroptionscomparison...)
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

        myhidexdecorations!(ax, irow < 2Nrows - 1)
        myhideydecorations!(ax, icol > 1)

        axs[irow, icol] = ax
    end
end

cb = Colorbar(
    fig[(Nrows + 1):(2Nrows - 1), 4], contours[Nrows + 1, 1];
    vertical = true, flipaxis = true,
    # ticks = (, cbarticklabelformat.(levels)),
    label = rich("Δ", Γdown, " (yr)"),
)
cb.height = Relative((2Nrows - 2) / (2Nrows - 1))

for (icol, (basin_str, xlims)) in enumerate(zip(basin_strs, basin_latlims))
    Label(fig[0, icol], basin_str, fontsize = 20, tellwidth = false)
    colsize!(fig.layout, icol, Auto(xlims[2] - xlims[1]))
end
Label(fig[1, 0], text = "cyclo", fontsize = 20, tellheight = false, rotation = π / 2)
Label(fig[2, 0], text = "mean(T)", fontsize = 20, tellheight = false, rotation = π / 2)
Label(fig[3, 0], text = "T(mean)", fontsize = 20, tellheight = false, rotation = π / 2)
Label(fig[4, 0], text = "mean(T) − cyclo", fontsize = 20, tellheight = false, rotation = π / 2)
Label(fig[5, 0], text = "T(mean) − cyclo", fontsize = 20, tellheight = false, rotation = π / 2)

title = "$model $experiment $member $(time_window) ideal age"
Label(fig[-1, 1:3], text = title, fontsize = 20, tellwidth = false)


rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)
rowgap!(fig.layout, 5, 20)

# save plot
outputfile = joinpath(inputdir, "ideal_age_ZAVGs_cyclo_vs_meanT_vs_Tmean.png")
@info "Saving ideal age ZAVGs as image file:\n  $(outputfile)"
save(outputfile, fig)


contouroptionsrelativecomparison = let
    levels = -100:10:100
    colormap = cgrad(:balance)
    # colormap = cgrad(colormap[2:end-1]; categorical=true)
    nan_color = :lightgray
    (; levels, colormap, nan_color)
end


Nrows = length(strs)
fig = Figure(size = (1200, 250 * (2Nrows - 1)), fontsize = 18)
axs = Array{Any, 2}(undef, (2Nrows - 1, 3))
contours = Array{Any, 2}(undef, (2Nrows - 1, 3))
for (icol, (basin_key, basin)) in enumerate(pairs(basins))

    for (irow, (x3D, str)) in enumerate(zip(Γ3Ds, strs))

        x2D = zonalaverage(x3D, gridmetrics; mask = basin)

        local ax = Axis(fig[irow, icol]; axisoptions...)

        co = contourf!(ax, lat2, zt, x2D; contouroptions...)

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

        myhidexdecorations!(ax, irow < 2Nrows - 1)
        myhideydecorations!(ax, icol > 1)

        axs[irow, icol] = ax
    end
end

Γdown = rich("Γ", superscript("↓"))
cb = Colorbar(
    fig[1:Nrows, 4], contours[1, 1];
    vertical = true, flipaxis = true,
    # ticks = (, cbarticklabelformat.(levels)),
    label = rich(Γdown, " (yr)"),
)
cb.height = Relative(Nrows / (Nrows + 1))


for (icol, (basin_key, basin)) in enumerate(pairs(basins))

    for (irow, (x3D, str)) in enumerate(zip(Γ3Ds[2:end], strs[2:end]))

        irow += Nrows

        x2D = zonalaverage(x3D - Γinyr3Dcyclomean, gridmetrics; mask = basin)
        x2Dbase = zonalaverage(Γinyr3Dcyclomean, gridmetrics; mask = basin)

        local ax = Axis(fig[irow, icol]; axisoptions...)

        co = contourf!(ax, lat2, zt, 100x2D ./ x2Dbase; contouroptionsrelativecomparison...)
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

        myhidexdecorations!(ax, irow < 2Nrows - 1)
        myhideydecorations!(ax, icol > 1)

        axs[irow, icol] = ax
    end
end

cb = Colorbar(
    fig[(Nrows + 1):(2Nrows - 1), 4], contours[Nrows + 1, 1];
    vertical = true, flipaxis = true,
    # ticks = (, cbarticklabelformat.(levels)),
    label = rich("Δ", Γdown, " / ", Γdown, " (%)"),
)
cb.height = Relative((2Nrows - 2) / (2Nrows - 1))

for (icol, (basin_str, xlims)) in enumerate(zip(basin_strs, basin_latlims))
    Label(fig[0, icol], basin_str, fontsize = 20, tellwidth = false)
    colsize!(fig.layout, icol, Auto(xlims[2] - xlims[1]))
end
Label(fig[1, 0], text = "cyclo", fontsize = 20, tellheight = false, rotation = π / 2)
Label(fig[2, 0], text = "mean(T)", fontsize = 20, tellheight = false, rotation = π / 2)
Label(fig[3, 0], text = "T(mean)", fontsize = 20, tellheight = false, rotation = π / 2)
Label(fig[4, 0], text = "mean(T) − cyclo", fontsize = 20, tellheight = false, rotation = π / 2)
Label(fig[5, 0], text = "T(mean) − cyclo", fontsize = 20, tellheight = false, rotation = π / 2)

title = "$model $experiment $member $(time_window) ideal age"
Label(fig[-1, 1:3], text = title, fontsize = 20, tellwidth = false)


rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)
rowgap!(fig.layout, 5, 20)

# save plot
outputfile = joinpath(inputdir, "ideal_age_ZAVGs_cyclo_vs_meanT_vs_Tmean_reldiff.png")
@info "Saving ideal age ZAVGs as image file:\n  $(outputfile)"
save(outputfile, fig)


# end


# # Plot Γ↑ zonal averages

# fig = Figure(size = (1200, 600), fontsize = 18)
# axs = Array{Any,2}(undef, (2, 3))
# contours = Array{Any,2}(undef, (2, 3))

# for (irow, (x3D, str)) in enumerate(zip((Γoutyr3D_mean, Γoutyr3D_std), ("mean", "std")))

#     if str == "mean" # mean
#         levels = 0:100:2500
#         colormap = cgrad(:viridis, length(levels); categorical=true)
#         extendlow = nothing
#         extendhigh = colormap[end]
#         colormap = cgrad(colormap[1:end-1]; categorical=true)
#     elseif str == "std"
#         levels = 0:10:200
#         colormap = cgrad(:viridis, length(levels); categorical=true)
#         extendlow = nothing
#         extendhigh = colormap[end]
#         colormap = cgrad(colormap[1:end-1]; categorical=true)
#     end

#     for (icol, (basin_key, basin)) in enumerate(pairs(basins))

#         x2D = zonalaverage(x3D, gridmetrics; mask = basin)

#         local ax = Axis(fig[irow, icol],
#             backgroundcolor=:lightgray,
#             xgridvisible=false, ygridvisible=false,
#             ylabel = "depth (m)")

#         X = dropdims(maximum(lat, dims=1), dims=1)
#         Y = zt
#         Z = x2D
#         co = contourf!(ax, X, Y, Z;
#             levels,
#             colormap,
#             nan_color = :lightgray,
#             extendlow,
#             extendhigh,
#         )
#         translate!(co, 0, 0, -100)
#         contours[irow, icol] = co

#         xlim = basin_latlims[basin_key]
#         # basin2 = LONGTEXT[basin]

#         ax.yticks = (ztick, zticklabel)
#         xticks = -90:30:90
#         ax.xticks = (xticks, latticklabel.(xticks))
#         ylims!(ax, zlim)
#         # xlims!(ax, (-90, 90))
#         xlims!(ax, xlim)

#         hidexdecorations!(ax,
#             label = irow < Nmonths + 1, ticklabels = irow < Nmonths + 1,
#             ticks = irow < Nmonths + 1, grid = false)
#         hideydecorations!(ax,
#             label = icol > 1, ticklabels = icol > 1,
#             ticks = icol > 1, grid = false)

#         axs[irow, icol] = ax
#     end
#     cb = Colorbar(fig[irow, 4], contours[irow, 1];
#         vertical = true, flipaxis = true,
#         # ticks = (, cbarticklabelformat.(levels)),
#         label = rich(str, " ", Γup, " (yr)"),
#         )
#     cb.height = Relative(1)

# end


# for (icol, (basin_str, xlims)) in enumerate(zip(basin_strs, basin_latlims))
#     Label(fig[0, icol], basin_str, fontsize=20, tellwidth=false)
#     colsize!(fig.layout, icol, Auto(xlims[2] - xlims[1]))
# end

# title = "$model $experiment all members $(time_window) reemergence time"
# Label(fig[-1, 1:3], text = title, fontsize=20, tellwidth=false)

# rowgap!(fig.layout, 10)
# colgap!(fig.layout, 10)

# # save plot
# outputfile = joinpath(outputdir, "reemergence_time_ZAVGs_v3.png")
# @info "Saving reemergence time ZAVGs as image file:\n  $(outputfile)"
# save(outputfile, fig)


# fig = Figure(size = (1200, 1200), fontsize = 18)
# axs = Array{Any,2}(undef, (2, 1))
# contours = Array{Any,2}(undef, (2, 1))

# for (irow, (x3D, str)) in enumerate(zip((Γinyr3D_mean, Γinyr3D_std), ("mean", "std")))

#     if str == "mean" # mean
#         levels = 0:100:2500
#         colormap = :viridis
#         colorrange = extrema(levels)
#     elseif str == "std"
#         levels = 0:10:200
#         colorrange = extrema(levels)
#         colormap = :viridis
#     end

#     # Plot mean age at the seafloor level
#     local title = "$model $experiment $member $(time_window) mean age at seafloor"

#     ax = Axis(fig[irow,1]; title, xtickformat, ytickformat)


#     # plot
#     x2D = seafloorvalue(x3D, wet3D)
#     plt = plotmap!(ax, x2D, gridmetrics; colorrange, colormap)

#     Colorbar(fig[irow,2], plt, label=rich(str, " ", Γdown, " at seafloor (yr)"))
# end

# rowgap!(fig.layout, 10)
# colgap!(fig.layout, 10)

# # save plot
# outputfile = joinpath(outputdir, "mean_age_at_seafloor_v3.png")
# @info "Saving ideal mean age at sea floor as image file:\n  $(outputfile)"
# save(outputfile, fig)


# fig = Figure(size = (1200, 1200), fontsize = 18)
# axs = Array{Any,2}(undef, (2, 1))
# contours = Array{Any,2}(undef, (2, 1))

# for (irow, (x3D, str)) in enumerate(zip((Γoutyr3D_mean, Γoutyr3D_std), ("mean", "std")))

#     if str == "mean" # mean
#         levels = 0:100:2500
#         colormap = :viridis
#         colorrange = extrema(levels)
#     elseif str == "std"
#         levels = 0:10:200
#         colorrange = extrema(levels)
#         colormap = :viridis
#     end

#     # Plot mean age at the seafloor level
#     local title = "$model $experiment $member $(time_window) reemergence time at seafloor"

#     ax = Axis(fig[irow,1]; title, xtickformat, ytickformat)


#     # plot
#     x2D = seafloorvalue(x3D, wet3D)
#     plt = plotmap!(ax, x2D, gridmetrics; colorrange, colormap)

#     # poly!(ax, reverse.(OCEANS[OceanBasins.atlantic()].polygon))
#     # poly!(ax, reverse.(OCEANS[OceanBasins.indian()].polygon))
#     # poly!(ax, reverse.(OCEANS[OceanBasins.east_pacific()].polygon))
#     # poly!(ax, reverse.(OCEANS[OceanBasins.west_pacific()].polygon))

#     Colorbar(fig[irow,2], plt, label=rich(str, " ", Γup, " at seafloor (yr)"))
# end

# rowgap!(fig.layout, 10)
# colgap!(fig.layout, 10)

# # save plot
# outputfile = joinpath(outputdir, "reemergence_time_at_seafloor_v3.png")
# @info "Saving mean reemergence time at sea floor as image file:\n  $(outputfile)"
# save(outputfile, fig)
