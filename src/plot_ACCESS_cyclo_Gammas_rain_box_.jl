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
# try
#     using CairoMakie
# catch
#     using CairoMakie
# end
# using GeoMakie
# using SwarmMakie
# using Interpolations
# using OceanBasins
# using NaNStatistics
# using Format

# include("plotting_functions.jl")

# model = "ACCESS-ESM1-5"

# experiment = "historical"
# # experiment = "piControl"

# time_window = "Jan1990-Dec1999"


# # Gadi directory for input files
# # inputdirfun(member) = "/scratch/xv83/TMIP/data/$model/$experiment/all members/$(time_window)"
# inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/all_members/$(time_window)/cyclomonth"
# outputdir = inputdir
# mkpath(inputdir)


# gridinputdir = "/scratch/xv83/TMIP/data/$model/$experiment/r1i1p1f1/$(time_window)"
# areacello_ds = open_dataset(joinpath(gridinputdir, "areacello.nc"))
# volcello_ds = open_dataset(joinpath(gridinputdir, "volcello.nc"))
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
# (; lon_vertices, lat_vertices, lon, lat, zt, v3D,) = gridmetrics
# lev = zt
# # Make indices
# indices = makeindices(gridmetrics.v3D)
# (; wet3D, N) = indices



# Γdown = rich("Γ", superscript("↓"))
# Γup = rich("Γ", superscript("↑"))




# # # Plot zonal averages


# @info "Loading age/adjoint age lazily"
# Γinyr3D_ds = open_dataset(joinpath(inputdir, "age_timemean.nc"))
# Γoutyr3D_ds = open_dataset(joinpath(inputdir, "adjointage_timemean.nc"))

# @info "Loading in memory"
# Γinyr3D = readcubedata(Γinyr3D_ds.age_timemean)
# Γoutyr3D = readcubedata(Γoutyr3D_ds.adjointage_timemean)

# basin_keys = (:GBL, :ATL, :PAC, :IND, :SO)
# basin_strs = ("Global", "Atlantic", "Pacific", "Indian", "Southern Ocean")
# basin_functions = ((lat, lon, O) -> trues(size(lat)), isatlantic, ispacific, isindian, (lat, lon, O) -> lat .≤ -40)
# basin_values = (reshape(f(lat[:], lon[:], OCEANS), size(lat)) for f in basin_functions)
# basins = (; (basin_keys .=> basin_values)...)


# Nmembers_Γin = length(Γinyr3D.member)
# data_Γin = reduce(vcat, average(Γinyr3D, gridmetrics; mask, dims=(1,2,3)) for (basin, mask) in pairs(basins))
# category_Γin = reduce(vcat, fill(string(basin), Nmembers_Γin) for (basin, mask) in pairs(basins))
# Nmembers_Γout = length(Γoutyr3D.member)
# data_Γout = reduce(vcat, average(Γoutyr3D, gridmetrics; mask, dims=(1,2,3)) for (basin, mask) in pairs(basins))
# category_Γout = reduce(vcat, fill(string(basin), Nmembers_Γout) for (basin, mask) in pairs(basins))

raincloudioptions = (
    boxplot_nudge = -0.3,
    boxplot_width = 0.3,
    jitter_width = 0.3,
    clouds = nothing,
    hist_bins = 0:5:3000,
    cloud_width = 0.3,
    gap = 0.1,
)

fig = Figure()
ax = Axis(fig[1, 1];
    # ylabel="basin",
    xticks=(1:length(basins), collect(basin_strs)),
    ylabel="Ideal mean age (yr)",
)
rainclouds!(ax, category_Γin, data_Γin; raincloudioptions...)
ylims!(ax, (0, 100ceil(maximum([data_Γin; data_Γout])/100)))
ax = Axis(fig[1, 2];
    # ylabel="basin",
    xticks=(1:length(basins), collect(basin_strs)),
    ylabel="Reemergence time (yr)",
)
rainclouds!(ax, category_Γout, data_Γout; raincloudioptions...)
ylims!(ax, (0, 100ceil(maximum([data_Γin; data_Γout])/100)))
fig
# save plot
outputfile = joinpath(outputdir, "Gammas_rainclouds.png")
@info "Saving ideal age rainclouds as image file:\n  $(outputfile)"
save(outputfile, fig)






raincloudioptions = (
    boxplot_width = 0.2,
    jitter_width = 0.1,
    clouds = nothing,
    hist_bins = 0:5:3000,
    # cloud_width = 0.3,
    center_boxplot = false,
    boxplot_nudge = -0.1,
    side_nudge = -0.3,
    # gap = 0.1,
    markersize = 3,
    orientation = :horizontal,
)

fig = Figure(size=(900, 600))
ax = Axis(fig[1, 1];
    # ylabel="basin",
    yticks=(1:length(basins), collect(basin_strs)),
    # ylabel="Ideal mean age/reemergence time (yr)",
    xlabel="years",
)
rainclouds!(ax, category_Γin, data_Γin; raincloudioptions...)
rainclouds!(ax, category_Γout, data_Γout; raincloudioptions..., side=:right)
xlims!(ax, (0, 100ceil(maximum([data_Γin; data_Γout])/100)))

text = rich("ideal age")
text!(ax, 20 + maximum(data_Γin[category_Γin .== "GBL"]), 1.15; text, align=(:left,:center), color=Makie.wong_colors()[1])
text = rich("reemergence time")
text!(ax, 20 + maximum(data_Γout[category_Γout .== "GBL"]), 0.85; text, align=(:left,:center), color=Makie.wong_colors()[2])
fig
# save plot
outputfile = joinpath(outputdir, "Gammas_rainclouds_v2.png")
@info "Saving ideal age rainclouds as image file:\n  $(outputfile)"
save(outputfile, fig)






# same but with only seafloor values

# SEAFLOORMASK = seafloormask(wet3D)
# data_Γin_seafloor = reduce(vcat, average(Γinyr3D, gridmetrics; mask = mask .* SEAFLOORMASK, dims=(1,2,3)) for (basin, mask) in pairs(basins))
# data_Γout_seafloor = reduce(vcat, average(Γoutyr3D, gridmetrics; mask = mask .* SEAFLOORMASK, dims=(1,2,3)) for (basin, mask) in pairs(basins))


fig = Figure(size=(900, 600))
ax = Axis(fig[1, 1];
    # ylabel="basin",
    yticks=(1:length(basins), collect(basin_strs)),
    # ylabel="Ideal mean age/reemergence time (yr)",
    xlabel="years",
)
rainclouds!(ax, category_Γin, data_Γin_seafloor; raincloudioptions...)
rainclouds!(ax, category_Γout, data_Γout_seafloor; raincloudioptions..., side=:right)
xlims!(ax, (0, 100ceil(maximum([data_Γin_seafloor; data_Γout_seafloor])/100)))

text = rich("ideal age at seafloor")
text!(ax, -20 + minimum(data_Γin_seafloor[category_Γin .== "GBL"]), 1.15; text, align=(:right,:center), color=Makie.wong_colors()[1])
text = rich("reemergence time at seafloor")
text!(ax, -20 + minimum(data_Γout_seafloor[category_Γout .== "GBL"]), 0.85; text, align=(:right,:center), color=Makie.wong_colors()[2])
fig
# save plot
outputfile = joinpath(outputdir, "Gammas_rainclouds_v2_seafloor.png")
@info "Saving ideal age rainclouds as image file:\n  $(outputfile)"
save(outputfile, fig)









# data_Γin_seafloor = reduce(vcat, seafloorvalue1D(x3D, wet3D) for x3D in eachslice(Γinyr3D, dims=:member))
# data_Γout_seafloor = reduce(vcat, seafloorvalue1D(x3D, wet3D) for x3D in eachslice(Γoutyr3D, dims=:member))
#





# begin
#     @show data = let
#         # dataGBL[:, :basin] .= 1
#         dataPAC[:, :basin] .= 1
#         dataIND[:, :basin] .= 2
#         dataATL[:, :basin] .= 3
#         # vcat(dataGBL, dataATL, dataPAC, dataIND)
#         vcat(dataATL, dataPAC, dataIND)
#     end

#     # Remove ACCESS1-0 r2i1p1 from the data because age = 0???
#     delete!(data, findall((data.model .== "ACCESS1-0") .& (data.member .== "r2i1p1")))

#     fig = Figure(size=(800, 300))
#     emptyax = Axis(fig[1, 1];
#         # ylabel="Basin",
#         # yticks = (1:length(basin_strs) + 1, ["Global", basin_strs...]),
#         yticks=(1:length(basin_strs), collect(basin_strs)),
#         xticksvisible=false, yticksvisible=false,
#         xticklabelsvisible=false,
#         # xticksize = 0, yticksize = 0,
#         xgridvisible = false, ygridvisible = false,
#         leftspinevisible = false, rightspinevisible = false,
#         topspinevisible = false, bottomspinevisible = false,
#         # yticks=(0.5:1:3.5, fill("", 4)),
#         # xlabel="Basin mean of ideal mean age (yr)",
#     )
#     ax = Axis(fig[1, 1];
#         # ylabel="Basin",
#         # yticks = (1:length(basin_strs) + 1, ["Global", basin_strs...]),
#         # yticks=(1:length(basin_strs), collect(basin_strs)),
#         yticks=(0.5:1:3.5, fill("", 4)),
#         yticksvisible=false,
#         xticks=0:100:1500, xtrimspine=true,
#         xlabel="Basin mean of ideal mean age (yr)",
#         xgridvisible = true, ygridvisible = false,
#         leftspinevisible = false, rightspinevisible = false,
#         topspinevisible = false, bottomspinevisible = true,
#     )
#     linkaxes!(emptyax, ax)
#     model_index(model) = findfirst(isequal(model), models)
#     cmap = cgrad(:seaborn_colorblind, alpha=1, categorical=true)
#     dodge_gap = 0.1
#     gap = 0.2
#     # think white layer to spearate basins
#     [hspan!(ax, y - gap/2, y + gap/2, color=:white) for y in 0.5:1:3.5]
#     boxplot!(ax, data.basin, data.mean_age;
#         dodge = model_index.(data.model),
#         orientation = :horizontal,
#         label=data.model,
#         markersize=4,
#         dodge_gap,
#         gap,
#         # color = (:blue, 0.2),
#         color=cmap[model_index.(data.model)],
#         # colormap = cgrad(:Archambault, alpha = 0.3),
#     )
#     # rainclouds!(ax, model_index.(dataGBL.model), dataGBL.mean_age;
#     #     # boxplot_width = 0.1,
#     #     # boxplot_nudge = 0.1,
#     #     clouds = hist,
#     #     center_boxplot = true,
#     #     # cloud_width = 0.1,
#     #     # gap = 0.5,
#     # )
#     # beeswarm!(ax, model_index.(dataGBL.model), dataGBL.mean_age;
#     #     markersize = 5,
#     #     orientation = :horizontal,
#     #     color = (:black, 0.3),
#     # )
#     ylims!(ax, (0.5, 3.5))
#     xlims!(ax, (10, nothing))
#     n_dodge = length(models)
#     @show dodge_width = Makie.scale_width(dodge_gap, n_dodge)
#     width = 1 - gap
#     for model in models
#         # Only label the Indian
#         # Pacific index is 1
#         x = minimum(data.mean_age[(data.basin .== 1) .& (data.model .== model)]) - 20
#         # x = maximum(data.mean_age[(data.basin .== 1)]) + 20
#         @show shift = Makie.shift_dodge(model_index.(model), dodge_width, dodge_gap)
#         y = 1 + width * shift
#         # y = 3 + Makie.shift_dodge(model_index.(model), dodge_width, 0.0\)
#         text!(ax, x, y; text=model, align=(:right, :center), color=cmap[model_index.(model)])
#     end

#     fig
# end
# # save plot
# outputfile = joinpath(TMIPDIR, "extra", "ideal_age_rainclouds.png")
# @info "Saving ideal age rainclouds as image file:\n  $(outputfile)"
# save(outputfile, fig)























