# qsub -I -P xv83 -l mem=64GB -l storage=scratch/gh0+scratch/xv83 -l walltime=02:00:00 -l ncpus=12

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

include("plotting_functions.jl")


model = "ACCESS-ESM1-5"

# for time_window in ["Jan1850-Dec1859", "Jan1990-Dec1999", "Jan2030-Dec2039", "Jan2090-Dec2099"]
time_window = "Jan2030-Dec2039"
    experiment = parse(Int, time_window[4:7]) ≤ 2010 ? "historical" : "ssp370"
    # experiment = "historical"
    # time_window = "Jan1850-Dec1859"
    # time_window = "Jan1990-Dec1999"
    # experiment = "ssp370"
    # time_window = "Jan2030-Dec2039"
    # time_window = "Jan2090-Dec2099"


    # Gadi directory for input files
    # inputdirfun(member) = "/scratch/xv83/TMIP/data/$model/$experiment/all members/$(time_window)"
    inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/all_members/$(time_window)/cyclomonth"
    outputdir = inputdir
    mkpath(inputdir)


    # Load areacello and volcello for grid geometry
    fixedvarsinputdir = "/scratch/xv83/TMIP/data/$model"
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
    (; lon_vertices, lat_vertices, lon, lat, zt, v3D, thkcello) = gridmetrics
    lev = zt
    # Make indices
    indices = makeindices(gridmetrics.v3D)
    (; wet3D, N) = indices



    Γdown = rich("Γ", superscript("↓"))
    Γup = rich("Γ", superscript("↑"))

    basin_keys = (:ATL, :PAC, :IND)
    basin_strs = ("Atlantic", "Pacific", "Indian")
    basin_functions = (isatlantic, ispacific, isindian)
    basin_values = (reshape(f(lat[:], lon[:], OCEANS), size(lat)) for f in basin_functions)
    basins = (; (basin_keys .=> basin_values)...)
    basin_latlims_values = [clamp.((-5, +5) .+ extrema(lat[.!isnan.(v3D[:,:,1]) .& basin[:,:,1]]), -80, 80) for basin in basins]
    basin_latlims = (; (basin_keys .=> basin_latlims_values)...)


    # Γinyr3D_mean = readcubedata(open_dataset(joinpath(inputdir, "age_ensemblemean.nc")).age_ensemblemean)
    # Γinyr3D_std = readcubedata(open_dataset(joinpath(inputdir, "age_ensemblestd.nc")).age_ensemblestd)
    # Γinyr3D_max = readcubedata(open_dataset(joinpath(inputdir, "age_ensemblemax.nc")).age_ensemblemax)
    # Γinyr3D_min = readcubedata(open_dataset(joinpath(inputdir, "age_ensemblemin.nc")).age_ensemblemin)
    # Γinyr3D_maxdiff = Γinyr3D_max - Γinyr3D_min
    Γoutyr3D_timemean = readcubedata(open_dataset(joinpath(inputdir, "adjointage_timemean.nc")).adjointage_timemean)
    Γoutyr3D_mean = readcubedata(open_dataset(joinpath(inputdir, "adjointage_ensemblemean.nc")).adjointage_ensemblemean)
    Γoutyr3D_std = readcubedata(open_dataset(joinpath(inputdir, "adjointage_ensemblestd.nc")).adjointage_ensemblestd)
    Γoutyr3D_max = readcubedata(open_dataset(joinpath(inputdir, "adjointage_ensemblemax.nc")).adjointage_ensemblemax)
    Γoutyr3D_min = readcubedata(open_dataset(joinpath(inputdir, "adjointage_ensemblemin.nc")).adjointage_ensemblemin)
    Γoutyr3D_maxdiff = Γoutyr3D_max - Γoutyr3D_min

    Γoutyr3D_argmin = dropdims(map(x -> Float64(x[4]), argmin(Γoutyr3D_timemean, dims = 4)), dims = 4)
    Γoutyr3D_argmax = dropdims(map(x -> Float64(x[4]), argmax(Γoutyr3D_timemean, dims = 4)), dims = 4)

    Γoutyr3D_argmin[.!wet3D] .= NaN
    Γoutyr3D_argmax[.!wet3D] .= NaN



    # # Plot zonal averages

    # # # Plot Γ↓ zonal averages

    # # fig = Figure(size = (1200, 600), fontsize = 18)
    # # axs = Array{Any,2}(undef, (2, 3))
    # # contours = Array{Any,2}(undef, (2, 3))

    # # for (irow, (x3D, str)) in enumerate(zip((Γinyr3D_mean, Γinyr3D_std), ("mean", "std")))

    # #     if str == "mean" # mean
    # #         levels = 0:100:1500
    # #         colormap = cgrad(:viridis, length(levels); categorical=true)
    # #         extendlow = nothing
    # #         extendhigh = colormap[end]
    # #         colormap = cgrad(colormap[1:end-1]; categorical=true)
    # #     elseif str == "std"
    # #         levels = 0:50:400
    # #         colormap = cgrad(:magma, length(levels); categorical=true)
    # #         extendlow = nothing
    # #         extendhigh = colormap[end]
    # #         colormap = cgrad(colormap[1:end-1]; categorical=true)
    # #     end

    # #     for (icol, (basin_key, basin)) in enumerate(pairs(basins))

    # #         x2D = zonalaverage(x3D, gridmetrics; mask = basin)

    # #         local ax = Axis(fig[irow, icol],
    # #             backgroundcolor=:lightgray,
    # #             xgridvisible=false, ygridvisible=false,
    # #             ylabel = "depth (m)")

    # #         X = dropdims(maximum(lat, dims=1), dims=1)
    # #         Y = zt
    # #         Z = x2D

    # #         co = contourf!(ax, X, Y, Z;
    # #             levels,
    # #             colormap,
    # #             nan_color = :lightgray,
    # #             extendlow,
    # #             extendhigh,
    # #         )
    # #         translate!(co, 0, 0, -100)
    # #         contours[irow, icol] = co

    # #         xlim = basin_latlims[basin_key]
    # #         # basin2 = LONGTEXT[basin]

    # #         ax.yticks = (ztick, zticklabel)
    # #         xticks = -90:30:90
    # #         ax.xticks = (xticks, latticklabel.(xticks))
    # #         ylims!(ax, zlim)
    # #         # xlims!(ax, (-90, 90))
    # #         xlims!(ax, xlim)

    # #         hidexdecorations!(ax,
    # #             label = irow < 2, ticklabels = irow < 2,
    # #             ticks = irow < 2, grid = false)
    # #         hideydecorations!(ax,
    # #             label = icol > 1, ticklabels = icol > 1,
    # #             ticks = icol > 1, grid = false)


    # #         axs[irow, icol] = ax

    # #     end

    # #     cb = Colorbar(fig[irow, 4], contours[irow, 1];
    # #         vertical = true, flipaxis = true,
    # #         # ticks = (, cbarticklabelformat.(levels)),
    # #         label = rich(str, " ", Γdown, " (yr)"),
    # #         )
    # #     cb.height = Relative(1)
    # # end


    # # for (icol, (basin_str, xlims)) in enumerate(zip(basin_strs, basin_latlims))
    # #     Label(fig[0, icol], basin_str, fontsize=20, tellwidth=false)
    # #     colsize!(fig.layout, icol, Auto(xlims[2] - xlims[1]))
    # # end

    # # title = "$model $experiment $(time_window) ideal age"
    # # Label(fig[-1, 1:3], text = title, fontsize=20, tellwidth=false)

    # # rowgap!(fig.layout, 10)
    # # colgap!(fig.layout, 10)

    # # # save plot
    # # outputfile = joinpath(outputdir, "ideal_age_ZAVGs_v4_$(time_window).png")
    # # @info "Saving ideal age ZAVGs as image file:\n  $(outputfile)"
    # # save(outputfile, fig)






    # # Plot Γ↑ zonal averages


    # fig = Figure(size = (1200, 600), fontsize = 18)
    # axs = Array{Any,2}(undef, (2, 3))
    # contours = Array{Any,2}(undef, (2, 3))

    # for (irow, (x3D, str)) in enumerate(zip((Γoutyr3D_mean, Γoutyr3D_std), ("mean", "std")))

    #     if str == "mean" # mean
    #         levels = 0:100:1500
    #         colormap = cgrad(:viridis, length(levels); categorical=true)
    #         extendlow = nothing
    #         extendhigh = colormap[end]
    #         colormap = cgrad(colormap[1:end-1]; categorical=true)
    #     elseif str == "std"
    #         levels = 0:50:400
    #         colormap = cgrad(:magma, length(levels); categorical=true)
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
    #             label = irow < 2, ticklabels = irow < 2,
    #             ticks = irow < 2, grid = false)
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

    # title = "$model $experiment $(time_window) reemergence time"
    # Label(fig[-1, 1:3], text = title, fontsize=20, tellwidth=false)

    # rowgap!(fig.layout, 10)
    # colgap!(fig.layout, 10)

    # # save plot
    # outputfile = joinpath(outputdir, "reemergence_time_ZAVGs_v4_$(time_window).png")
    # @info "Saving reemergence time ZAVGs as image file:\n  $(outputfile)"
    # save(outputfile, fig)









    # # fig = Figure(size = (1200, 1200), fontsize = 18)
    # # axs = Array{Any,2}(undef, (2, 1))
    # # contours = Array{Any,2}(undef, (2, 1))

    # # for (irow, (x3D, str)) in enumerate(zip((Γinyr3D_mean, Γinyr3D_std), ("mean", "std")))

    # #     if str == "mean" # mean
    # #         levels = 0:100:1500
    # #         colormap = :viridis
    # #         colorrange = extrema(levels)
    # #     elseif str == "std"
    # #         levels = 0:50:400
    # #         colorrange = extrema(levels)
    # #         colormap = :magma
    # #     end

    # #     # Plot mean age at the seafloor level
    # #     local title = "$model $experiment $(time_window) mean age at seafloor"

    # #     ax = Axis(fig[irow,1]; title, xtickformat, ytickformat)


    # #     # plot
    # #     x2D = seafloorvalue(x3D, wet3D, gridmetrics)
    # #     plt = plotmap!(ax, x2D, gridmetrics; colorrange, colormap)

    # #     Colorbar(fig[irow,2], plt, label=rich(str, " ", Γdown, " at seafloor (yr)"))
    # # end

    # # rowgap!(fig.layout, 10)
    # # colgap!(fig.layout, 10)

    # # # save plot
    # # outputfile = joinpath(outputdir, "mean_age_at_seafloor_v4_$(time_window).png")
    # # @info "Saving ideal mean age at sea floor as image file:\n  $(outputfile)"
    # # save(outputfile, fig)















    # fig = Figure(size = (1200, 1200), fontsize = 18)
    # axs = Array{Any,2}(undef, (2, 1))
    # contours = Array{Any,2}(undef, (2, 1))

    # for (irow, (x3D, str)) in enumerate(zip((Γoutyr3D_mean, Γoutyr3D_std), ("mean", "std")))

    #     if str == "mean" # mean
    #         levels = 0:100:1500
    #         colormap = :viridis
    #         colorrange = extrema(levels)
    #     elseif str == "std"
    #         levels = 0:50:400
    #         colorrange = extrema(levels)
    #         colormap = :magma
    #     end

    #     # Plot mean age at the seafloor level
    #     local title = "$model $experiment $(time_window) reemergence time at seafloor"

    #     ax = Axis(fig[irow,1]; title, xtickformat, ytickformat)


    #     # plot
    #     x2D = seafloorvalue(x3D, wet3D, gridmetrics)
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
    # outputfile = joinpath(outputdir, "reemergence_time_at_seafloor_v4_$(time_window).png")
    # @info "Saving mean reemergence time at sea floor as image file:\n  $(outputfile)"
    # save(outputfile, fig)











    # # Redo the same for max minus min


    # # # Plot Γ↓ zonal averages

    # # fig = Figure(size = (1200, 600), fontsize = 18)
    # # axs = Array{Any,2}(undef, (2, 3))
    # # contours = Array{Any,2}(undef, (2, 3))

    # # for (irow, (x3D, str)) in enumerate(zip((Γinyr3D_mean, Γinyr3D_maxdiff), ("mean", "maxΔ")))

    # #     if str == "mean" # mean
    # #         levels = 0:100:1500
    # #         colormap = cgrad(:viridis, length(levels); categorical=true)
    # #         extendlow = nothing
    # #         extendhigh = colormap[end]
    # #         colormap = cgrad(colormap[1:end-1]; categorical=true)
    # #     elseif str == "maxΔ"
    # #         levels = 0:50:400
    # #         colormap = cgrad(:magma, length(levels); categorical=true)
    # #         extendlow = nothing
    # #         extendhigh = colormap[end]
    # #         colormap = cgrad(colormap[1:end-1]; categorical=true)
    # #     end

    # #     for (icol, (basin_key, basin)) in enumerate(pairs(basins))

    # #         x2D = zonalaverage(x3D, gridmetrics; mask = basin)

    # #         local ax = Axis(fig[irow, icol],
    # #             backgroundcolor=:lightgray,
    # #             xgridvisible=false, ygridvisible=false,
    # #             ylabel = "depth (m)")

    # #         X = dropdims(maximum(lat, dims=1), dims=1)
    # #         Y = zt
    # #         Z = x2D

    # #         co = contourf!(ax, X, Y, Z;
    # #             levels,
    # #             colormap,
    # #             nan_color = :lightgray,
    # #             extendlow,
    # #             extendhigh,
    # #         )
    # #         translate!(co, 0, 0, -100)
    # #         contours[irow, icol] = co

    # #         xlim = basin_latlims[basin_key]
    # #         # basin2 = LONGTEXT[basin]

    # #         ax.yticks = (ztick, zticklabel)
    # #         xticks = -90:30:90
    # #         ax.xticks = (xticks, latticklabel.(xticks))
    # #         ylims!(ax, zlim)
    # #         # xlims!(ax, (-90, 90))
    # #         xlims!(ax, xlim)

    # #         hidexdecorations!(ax,
    # #             label = irow < 2, ticklabels = irow < 2,
    # #             ticks = irow < 2, grid = false)
    # #         hideydecorations!(ax,
    # #             label = icol > 1, ticklabels = icol > 1,
    # #             ticks = icol > 1, grid = false)


    # #         axs[irow, icol] = ax

    # #     end

    # #     cb = Colorbar(fig[irow, 4], contours[irow, 1];
    # #         vertical = true, flipaxis = true,
    # #         # ticks = (, cbarticklabelformat.(levels)),
    # #         label = rich(str, " ", Γdown, " (yr)"),
    # #         )
    # #     cb.height = Relative(1)
    # # end


    # # for (icol, (basin_str, xlims)) in enumerate(zip(basin_strs, basin_latlims))
    # #     Label(fig[0, icol], basin_str, fontsize=20, tellwidth=false)
    # #     colsize!(fig.layout, icol, Auto(xlims[2] - xlims[1]))
    # # end

    # # title = "$model $experiment $(time_window) ideal age"
    # # Label(fig[-1, 1:3], text = title, fontsize=20, tellwidth=false)

    # # rowgap!(fig.layout, 10)
    # # colgap!(fig.layout, 10)

    # # # save plot
    # # outputfile = joinpath(outputdir, "ideal_age_ZAVGs_vsmaxdiff_$(time_window).png")
    # # @info "Saving ideal age ZAVGs as image file:\n  $(outputfile)"
    # # save(outputfile, fig)






    # Plot Γ↑ zonal averages


    # fig = Figure(size = (1200, 600), fontsize = 18)
    # axs = Array{Any,2}(undef, (2, 3))
    # contours = Array{Any,2}(undef, (2, 3))

    # for (irow, (x3D, str)) in enumerate(zip((Γoutyr3D_mean, Γoutyr3D_maxdiff), ("mean", "maxΔ")))

    #     if str == "mean" # mean
    #         levels = 0:100:1500
    #         colormap = cgrad(:viridis, length(levels); categorical=true)
    #         extendlow = nothing
    #         extendhigh = colormap[end]
    #         colormap = cgrad(colormap[1:end-1]; categorical=true)
    #     elseif str == "maxΔ"
    #         levels = 0:50:400
    #         colormap = cgrad(:magma, length(levels); categorical=true)
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
    #             label = irow < 2, ticklabels = irow < 2,
    #             ticks = irow < 2, grid = false)
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

    # title = "$model $experiment $(time_window) reemergence time"
    # Label(fig[-1, 1:3], text = title, fontsize=20, tellwidth=false)

    # rowgap!(fig.layout, 10)
    # colgap!(fig.layout, 10)

    # # save plot
    # outputfile = joinpath(outputdir, "reemergence_time_ZAVGs_vsmaxdiff_$(time_window).png")
    # @info "Saving reemergence time ZAVGs as image file:\n  $(outputfile)"
    # save(outputfile, fig)









    # # fig = Figure(size = (1200, 1200), fontsize = 18)
    # # axs = Array{Any,2}(undef, (2, 1))
    # # contours = Array{Any,2}(undef, (2, 1))

    # # for (irow, (x3D, str)) in enumerate(zip((Γinyr3D_mean, Γinyr3D_maxdiff), ("mean", "maxΔ")))

    # #     if str == "mean" # mean
    # #         levels = 0:100:1500
    # #         colormap = :viridis
    # #         colorrange = extrema(levels)
    # #     elseif str == "maxΔ"
    # #         levels = 0:50:400
    # #         colorrange = extrema(levels)
    # #         colormap = :magma
    # #     end

    # #     # Plot mean age at the seafloor level
    # #     local title = "$model $experiment $(time_window) mean age at seafloor"

    # #     ax = Axis(fig[irow,1]; title, xtickformat, ytickformat)


    # #     # plot
    # #     x2D = seafloorvalue(x3D, wet3D, gridmetrics)
    # #     plt = plotmap!(ax, x2D, gridmetrics; colorrange, colormap)

    # #     Colorbar(fig[irow,2], plt, label=rich(str, " ", Γdown, " at seafloor (yr)"))
    # # end

    # # rowgap!(fig.layout, 10)
    # # colgap!(fig.layout, 10)

    # # # save plot
    # # outputfile = joinpath(outputdir, "mean_age_at_seafloor_vsmaxdiff_$(time_window).png")
    # # @info "Saving ideal mean age at sea floor as image file:\n  $(outputfile)"
    # # save(outputfile, fig)















    # fig = Figure(size = (800, 800), fontsize = 18)
    # axs = Array{Any,2}(undef, (2, 1))
    # contours = Array{Any,2}(undef, (2, 1))
    # yticks = -60:30:60

    # strs = [
    #     "$(time_window[4:7])s Mean Reemergence Time"
    #     "Internal Variability"
    # ]

    # for (irow, x3D) in enumerate((Γoutyr3D_mean, Γoutyr3D_maxdiff))

    #     if irow == 1 # mean
    #         levels = 0:100:1500
    #         colormap = :viridis
    #         colorrange = extrema(levels)
    #         label = rich("seafloor ensemble mean ", Γup, " (yr)")
    #     else
    #         levels = 0:50:250
    #         colorrange = extrema(levels)
    #         colormap = :magma
    #         label = rich("seafloor ensemble max − min ", Γup, " (yr)")
    #     end

    #     # Plot mean age at the seafloor level
    #     local title = "$model $experiment $(time_window) reemergence time at seafloor"

    #     axs[irow, 1] = ax = Axis(fig[irow,1]; yticks, xtickformat, ytickformat)

    #     # plot
    #     x2D = seafloorvalue(x3D, wet3D, gridmetrics)
    #     plt = plotmap!(ax, x2D, gridmetrics; colorrange, colormap) # <- need to fix wrapping longitude for contour levels

    #     myhidexdecorations!(ax, irow < 2)

    #     # poly!(ax, reverse.(OCEANS[OceanBasins.atlantic()].polygon))
    #     # poly!(ax, reverse.(OCEANS[OceanBasins.indian()].polygon))
    #     # poly!(ax, reverse.(OCEANS[OceanBasins.east_pacific()].polygon))
    #     # poly!(ax, reverse.(OCEANS[OceanBasins.west_pacific()].polygon))

    #     cb = Colorbar(fig[irow,2], plt; label)
    #     cb.height = Relative(2/3)

    #     text = strs[irow]
    #     Label(fig[irow, 0]; text, rotation = π/2, tellheight = false, fontsize = 24)
    # end

    # for (ax, label) in zip(axs, ["a", "b"])
    #     text!(ax, 0, 1; text = label, labeloptions..., strokecolor = :white, strokewidth = 3)
    #     text!(ax, 0, 1; text = label, labeloptions...)
    # end

    # # rowgap!(fig.layout, 10)
    # # colgap!(fig.layout, 10)

    # # save plot
    # outputfile = joinpath(outputdir, "reemergence_time_at_seafloor_vsmaxdiff_$(time_window)_v2.png")
    # @info "Saving mean reemergence time at sea floor as image file:\n  $(outputfile)"
    # save(outputfile, fig)





    # Same but with elevation plot as well
    axs = Array{Any,2}(undef, (3, 1))
    contours = Array{Any,2}(undef, (3, 1))
    fig = Figure(size = (800, size(axs, 1) * 300), fontsize = 18)
    yticks = -60:30:60
    xticks = -180:60:180

    strs = [
        # "$(time_window[4:7])s Mean reemergence time"
        "Ensemble Mean"
        "Internal Variability"
        "Fractional Internal Variability"
    ]

    # for (irow, x3D) in enumerate((cumsum(thkcello, dims=3), Γoutyr3D_mean, Γoutyr3D_maxdiff, Γoutyr3D_maxdiff ./ Γoutyr3D_mean))
    for (irow, x3D) in enumerate((Γoutyr3D_mean, Γoutyr3D_maxdiff, 100 * Γoutyr3D_maxdiff ./ Γoutyr3D_mean))

        # if irow == 1 # bathy
        #     levels = [0, 6000] # Note levels seems unused here (except for colorrange)
        #     colormap = cgrad(:linear_blue_5_95_c73_n256, 6, categorical = true)
        #     colorrange = extrema(levels)
        #     label = rich("seafloor depth (m)")
        # elseif irow == 2 # mean
        if irow == 1 # mean
            levels = 0:100:1500
            colormap = cgrad(:linear_bgy_10_95_c74_n256, 15, categorical = true)
            colorrange = extrema(levels)
            label = rich("seafloor ensemble mean ", Γup, " (yr)")
        elseif irow == 2 # mean
            levels = 0:50:250
            colorrange = extrema(levels)
            colormap = cgrad(:linear_bmy_10_95_c71_n256, 5, categorical = true)
            label = rich("ensemble max − min (yr)")
        else
            levels = 100 .* (0:0.05:0.5)
            colormap = cgrad(:linear_bmw_5_95_c86_n256, 5, categorical = true)
            colorrange = extrema(levels)
            label = rich("(max − min) / mean (%)")
        end

        axs[irow, 1] = ax = Axis(fig[irow,1]; yticks, xticks, xtickformat, ytickformat)

        # plot
        x2D = seafloorvalue(x3D, wet3D, gridmetrics)
        plt = plotmap!(ax, x2D, gridmetrics; colorrange, colormap) # <- need to fix wrapping longitude for contour levels

        myhidexdecorations!(ax, irow < size(axs, 1))

        # poly!(ax, reverse.(OCEANS[OceanBasins.atlantic()].polygon))
        # poly!(ax, reverse.(OCEANS[OceanBasins.indian()].polygon))
        # poly!(ax, reverse.(OCEANS[OceanBasins.east_pacific()].polygon))
        # poly!(ax, reverse.(OCEANS[OceanBasins.west_pacific()].polygon))

        cb = Colorbar(fig[irow,2], plt; label)
        cb.height = Relative(2/3)

        text = strs[irow]
        Label(fig[irow, 0]; text, rotation = π/2, tellheight = false, fontsize = 24)
    end

    for (ax, label) in zip(axs, ["a", "b", "c"])
        txt = text!(ax, 0, 1; text = "$label", labeloptions..., strokecolor = :white, strokewidth = 3)
        translate!(txt, 0, 0, 100)
        txt = text!(ax, 0, 1; text = "$label", labeloptions...)
        translate!(txt, 0, 0, 100)
    end

    Label(fig[0, 1]; text = "$(time_window[4:7])s Mean reemergence time", fontsize = 24, tellwidth = false)
    # rowgap!(fig.layout, 10)
    # colgap!(fig.layout, 10)

    # save plot
    outputfile = joinpath(outputdir, "reemergence_time_at_seafloor_vsmaxdiff_$(time_window)_v4.png")
    @info "Saving mean reemergence time at sea floor as image file:\n  $(outputfile)"
    save(outputfile, fig)




    # # Quick plot of joint PDF of max−min vs mean
    # fig = Figure()
    # # x = seafloorvalue(Γoutyr3D_mean, wet3D, gridmetrics)
    # # y = seafloorvalue(Γoutyr3D_maxdiff, wet3D, gridmetrics)
    # x = seafloorvalue(Γoutyr3D_mean .* (Z3D .> 3000), wet3D, gridmetrics)
    # y = seafloorvalue(Γoutyr3D_maxdiff .* (Z3D .> 3000), wet3D, gridmetrics)
    # # x = Γoutyr3D_mean[:,:,21:end]
    # # y = Γoutyr3D_maxdiff[:,:,21:end]
    # ikeep = findall(x -> !isnan(x) & !iszero(x), x)
    # x = x[ikeep]
    # y = y[ikeep]
    # ax = Axis(fig[1,1], xlabel = "ensemble mean", ylabel = "count", limits = (0, nothing, 0, nothing))
    # hist!(ax, x)
    # ax = Axis(fig[2,2], xlabel = "count", ylabel = "max minus min", limits = (0, nothing, 0, nothing))
    # hist!(ax, y; direction = :x)
    # ax = Axis(fig[2,1], xlabel = "ensemble mean", ylabel = "max minus min", limits = (0, nothing, 0, nothing))
    # points = Makie.StructArray{Point2f}((x, y))
    # ds = datashader!(ax, points; colormap = :binary, async = false)
    # for slope in [0.1, 0.2, 0.5, 1]
    #     ablines!(ax, 0, slope)
    # end
    # slopes = y ./ x
    # ablines!(ax, 0, mean(slopes), linestyle = :dash)
    # ablines!(ax, 0, median(slopes), linestyle = :dash)
    # ax = Axis(fig[1,2], xlabel = "slope", ylabel = "count", limits = (-1, 2, 0, nothing))
    # hist!(ax, slopes, bins = -1:0.05:2, color=:red)
    # for slope in [0.1, 0.2, 0.5, 1]
    #     vlines!(ax, slope)
    # end
    # vlines!(ax, mean(slopes), label = "mean slope", linestyle = :dash)
    # vspan!(ax, (mean(slopes) .+ std(slopes) * [-1, 1])..., color=(:blue, 0.1))
    # vlines!(ax, median(slopes), label = "median slope", linestyle = :dash)
    # # save plot
    # outputfile = joinpath(outputdir, "jointPDF_mean_vs_maxminusmin_adjointage_$(time_window).png")
    # @info "Saving mean reemergence time at sea floor as image file:\n  $(outputfile)"
    # save(outputfile, fig)
    # @show mean(slopes)
    # @show median(slopes)






    # # Plot argmin and argmax to see which members

    # fig = Figure(size = (800, 800), fontsize = 18)
    # axs = Array{Any,2}(undef, (2, 1))
    # contours = Array{Any,2}(undef, (2, 1))
    # yticks = -60:30:60

    # strs = [
    #     rich("Member with Minimal ", Γup)
    #     rich("Member with Maximal ", Γup)
    # ]

    # for (irow, x3D) in enumerate((Γoutyr3D_argmin, Γoutyr3D_argmax))

    #     levels = 0.5:40.5
    #     colormap = cgrad([collect(cgrad(:tab20b, categorical = true)); collect(cgrad(:tab20c, categorical = true))]; categorical = true)
    #     colorrange = extrema(levels)
    #     label = rich("member")

    #     # Plot mean age at the seafloor level
    #     axs[irow, 1] = ax = Axis(fig[irow,1]; yticks, xtickformat, ytickformat)

    #     # plot
    #     x2D = seafloorvalue(x3D, wet3D)
    #     plt = plotmap!(ax, x2D, gridmetrics; colorrange, colormap) # <- need to fix wrapping longitude for contour levels

    #     myhidexdecorations!(ax, irow < 2)

    #     # poly!(ax, reverse.(OCEANS[OceanBasins.atlantic()].polygon))
    #     # poly!(ax, reverse.(OCEANS[OceanBasins.indian()].polygon))
    #     # poly!(ax, reverse.(OCEANS[OceanBasins.east_pacific()].polygon))
    #     # poly!(ax, reverse.(OCEANS[OceanBasins.west_pacific()].polygon))


    #     text = strs[irow]
    #     Label(fig[irow, 0]; text, rotation = π/2, tellheight = false, fontsize = 24)

    #     if irow == 2
    #         cb = Colorbar(fig[:,2], plt; label, ticks = 1:4:40)
    #         cb.height = Relative(2/3)
    #     end
    # end


    # for (ax, label) in zip(axs, ["a", "b"])
    #     text!(ax, 0, 1; text = label, labeloptions..., strokecolor = :white, strokewidth = 3)
    #     text!(ax, 0, 1; text = label, labeloptions...)
    # end

    # # rowgap!(fig.layout, 10)
    # # colgap!(fig.layout, 10)

    # # save plot
    # outputfile = joinpath(outputdir, "reemergence_time_argminmax_$(time_window).png")
    # @info "Saving mean reemergence time at sea floor as image file:\n  $(outputfile)"
    # save(outputfile, fig)








    # # # Same but show min and max separately

    # # axs = Array{Any,2}(undef, (3, 1))
    # # contours = Array{Any,2}(undef, (3, 1))
    # # fig = Figure(size = (1200, size(axs, 1) * 600), fontsize = 18)

    # # for (irow, (x3D, str)) in enumerate(zip((Γoutyr3D_mean, Γoutyr3D_max, Γoutyr3D_min), ("mean", "max", "min")))

    # #     if str ∈ ("mean", "max", "min") # mean
    # #         levels = 0:100:1500
    # #         colormap = cgrad(:viridis, length(levels) - 1; categorical = true)
    # #         colorrange = extrema(levels)
    # #     elseif str == "maxΔ"
    # #         levels = 0:50:400
    # #         colorrange = extrema(levels)
    # #         colormap = :magma
    # #     end

    # #     # Plot mean age at the seafloor level
    # #     local title = "$model $experiment $(time_window) reemergence time at seafloor"

    # #     ax = Axis(fig[irow,1]; title, xtickformat, ytickformat)


    # #     # plot
    # #     x2D = seafloorvalue(x3D, wet3D, gridmetrics)
    # #     plt = plotmap!(ax, x2D, gridmetrics; colorrange, colormap)

    # #     # poly!(ax, reverse.(OCEANS[OceanBasins.atlantic()].polygon))
    # #     # poly!(ax, reverse.(OCEANS[OceanBasins.indian()].polygon))
    # #     # poly!(ax, reverse.(OCEANS[OceanBasins.east_pacific()].polygon))
    # #     # poly!(ax, reverse.(OCEANS[OceanBasins.west_pacific()].polygon))

    # #     contours[irow] = plt

    # # end
    # # Colorbar(fig[:,2], contours[1], label=rich(Γup, " at seafloor (yr)"))

    # # rowgap!(fig.layout, 10)
    # # colgap!(fig.layout, 10)

    # # # save plot
    # # outputfile = joinpath(outputdir, "reemergence_time_at_seafloor_meanmaxmin_$(time_window).png")
    # # @info "Saving mean reemergence time at sea floor as image file:\n  $(outputfile)"
    # # save(outputfile, fig)







    # # Same but show min and max separately
    # # FOcus on AU only
    # # and ploit injection locations

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

    # axs = Array{Any,2}(undef, (3, 1))
    # contours = Array{Any,2}(undef, (3, 1))
    # fig = Figure(size = (600, size(axs, 1) * 400), fontsize = 18)
    # limits = ((107, 158), (-48, -7))

    # for (irow, (x3D, str)) in enumerate(zip((Γoutyr3D_mean, Γoutyr3D_max, Γoutyr3D_min), ("mean", "max", "min")))

    #     if str ∈ ("mean", "max", "min") # mean
    #         # levels = 0:100:1500
    #         levels = 0:100:1200
    #         colormap = cgrad(:viridis, length(levels) - 1; categorical = true)
    #         colorrange = extrema(levels)
    #     elseif str == "maxΔ"
    #         levels = 0:50:400
    #         colorrange = extrema(levels)
    #         colormap = :magma
    #     end

    #     # Plot mean age at the seafloor level
    #     local title = "$model $experiment $(time_window) reemergence time at seafloor"

    #     ax = Axis(fig[irow,1]; xtickformat, ytickformat)


    #     # plot
    #     x2D = seafloorvalue(x3D, wet3D, gridmetrics)
    #     # plt = plotmap!(ax, x2D, gridmetrics; colorrange, colormap)
    #     plt = plotmap!(ax, x2D, gridmetrics; colorrange, colormap, levels)

    #     # poly!(ax, reverse.(OCEANS[OceanBasins.atlantic()].polygon))
    #     # poly!(ax, reverse.(OCEANS[OceanBasins.indian()].polygon))
    #     # poly!(ax, reverse.(OCEANS[OceanBasins.east_pacific()].polygon))
    #     # poly!(ax, reverse.(OCEANS[OceanBasins.west_pacific()].polygon))

    #     # Add injection locations
    #     colors = cgrad(:Egypt, categorical=true)[[3, 4, 1]]
    #     # colors = Makie.wong_colors()[[1, 3, 6]]
    #     offsets = map(x -> x.* 2, [(-2, 1), (-2, -1), (2, -1)])
    #     # aligns = [(:right, :bottom), (:right, :top), (:left, :top)]
    #     aligns = [(:right, :center), (:right, :center), (:left, :center)]
    #     texts = ["A", "B", "C"]
    #     srcnames = ["Karratha", "Portland", "Marlo"]

    #     for (ksrc, (srcname, offset, align, color, text)) in enumerate(zip(srcnames, offsets, aligns, colors, texts))
    #         src_P = sourcelocation(srcname)
    #         # sc = scatter!(ax, src_P; marker=:star5, markersize=20, color=colors[ksrc], strokecolor=:black, strokewidth=1)
    #         # sc1 = scatter!(ax, src_P; marker=:circle, markersize=10, color=(:black, 0), strokecolor=:black, strokewidth=3)
    #         sc2 = scatter!(ax, src_P; marker=:circle, markersize=10, color=(:black, 0), strokecolor=:black, strokewidth=4)
    #         sc2 = scatter!(ax, src_P; marker=:circle, markersize=10, color=(:black, 0), strokecolor=color, strokewidth=2)
    #         # lines!(ax, [src_P, src_P .+ offset]; color=:white)
    #         lines!(ax, kinkline(src_P .+ offset, src_P); color=:black, linewidth = 1)
    #         # lines!(ax, kinkline(src_P .+ offset, src_P); color=:black, linewidth=3)
    #         # lines!(ax, kinkline(src_P .+ offset, src_P); color)
    #         text!(ax, src_P .+ offset; text, align, color=:black, strokecolor=:black)
    #         # text!(ax, src_P .+ offset; text, align, color=:black, font=:bold, fontsize=18, strokecolor=:black, strokewidth=2)
    #         # text!(ax, src_P .+ offset; text, align, color, font=:bold, fontsize=18)
    #         # translate!(sc1, 0, 0, 99)
    #         translate!(sc2, 0, 0, 200)
    #     end

    #     contours[irow] = plt

    #     myhidexdecorations!(ax, irow < 3)

    #     xlims!(ax, limits[1])
    #     ylims!(ax, limits[2])

    # end
    # cb = Colorbar(fig[:,2], contours[1], label=rich(Γup, " at seafloor (yr)"))
    # cb.height = Relative(0.6)


    # rowgap!(fig.layout, 10)
    # colgap!(fig.layout, 10)

    # # save plot
    # outputfile = joinpath(outputdir, "reemergence_time_at_seafloor_meanmaxmin_AU_$(time_window).png")
    # @info "Saving mean reemergence time at sea floor as image file:\n  $(outputfile)"
    # save(outputfile, fig)


# end