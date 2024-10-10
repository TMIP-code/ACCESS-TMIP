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
using SwarmMakie
using Interpolations
using OceanBasins
using NaNStatistics
using Format

include("plotting_functions.jl")

models = ["ACCESS1-0", "ACCESS1-3", "ACCESS-ESM1-5", "ACCESS-CM2"]

CMIP5_models = ["ACCESS1-0", "ACCESS1-3"]

experiment = "historical"
# experiment = "piControl"

time_window = "Jan1990-Dec1999"
# time_window = "Jan1071-Dec1100" # <- last 30 years of ACCESS-ESM1-5 piControl
# time_window = "Jan1420-Dec1449" # <- last 30 years of ACCESS-CM2 piControl


# Directory for data
if gethostname() == "benoits-MacBook-Pro.local"
    TMIPDIR = "/Users/benoitpasquier/Data/TMIP"
else # on Gadi. TODO: change this to the correct path
    TMIPDIR = "/scratch/xv83/TMIP"
end
DATADIR = joinpath(TMIPDIR, "data")

requiredvariables = ["ideal_mean_age", "volcello"]



basin_keys = (:PAC, :IND, :ATL)
basin_strs = ("Pacific", "Indian", "Atlantic")
basin_functions = (ispacific, isindian, isatlantic)

dataGBL = DataFrame(model=String[], member=String[], mean_age=Float64[])
dataATL = DataFrame(model=String[], member=String[], mean_age=Float64[])
dataPAC = DataFrame(model=String[], member=String[], mean_age=Float64[])
dataIND = DataFrame(model=String[], member=String[], mean_age=Float64[])
# Fetch the ideal mean age + volume and compute the global mean of the ideal mean age
for model in models

    members = readdir(joinpath(DATADIR, "$model/$experiment"))
    members = [m for m in members if m ≠ ".DS_Store"]

    # sort members by r, i, p[, f]
    memmber_regex = model ∈ CMIP5_models ? r"r(\d+)i(\d+)p(\d+)" : r"r(\d+)i(\d+)p(\d+)f(\d+)"
    parse_member(member) = parse.(Int, match(memmber_regex, member).captures)
    members = sort(members, by=x -> parse_member(x))

    println("$model")

    for member in members

        inputdir = joinpath(DATADIR, "$model/$experiment/$member/$(time_window)")
        println("  $member")

        # Exit loop (`continue`) if any of the required variables is missing
        hasrequireddata(variable_name) = isfile(joinpath(inputdir, "$variable_name.nc"))
        all(hasrequireddata(v) for v in requiredvariables) || continue

        # Load age and volume
        age = open_dataset(joinpath(inputdir, "ideal_mean_age.nc")).age |> Array
        volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
        areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))
        mlotst_ds = open_dataset(joinpath(inputdir, "mlotst.nc"))
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

        modelgrid = makemodelgrid(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
        (; lon, lat) = modelgrid

        # Save the global mean of the ideal mean age
        _FillValue = get(volcello.properties, "_FillValue", NaN)
        volecllo = volcello |> Array
        volcello = replace(volcello, missing => NaN, _FillValue => NaN)
        mean_age = nansum(age .* volcello) ./ nansum(volcello)
        push!(dataGBL, (model, member, mean_age))
        # println(" $mean_age $(nansum(age)) $(nansum(volcello))")

        basin_values = (reshape(f(lat[:], lon[:], OCEANS), size(lat)) for f in basin_functions)
        basins = (; (basin_keys .=> basin_values)...)

        ATL_mean_age = nansum(age .* volcello .* basins.ATL) ./ nansum(volcello .* basins.ATL)
        PAC_mean_age = nansum(age .* volcello .* basins.PAC) ./ nansum(volcello .* basins.PAC)
        IND_mean_age = nansum(age .* volcello .* basins.IND) ./ nansum(volcello .* basins.IND)
        push!(dataATL, (model, member, ATL_mean_age))
        push!(dataPAC, (model, member, PAC_mean_age))
        push!(dataIND, (model, member, IND_mean_age))

    end

end

fig = Figure()
ax = Axis(fig[1, 1];
    xlabel="Model",
    xticks=(1:length(models), models),
    ylabel="Global mean ideal mean age (yr)",
)
model_index(model) = findfirst(isequal(model), models)
# boxplot!(ax, model_index.(dataGBL.model), dataGBL.mean_age)
rainclouds!(ax, model_index.(dataGBL.model), dataGBL.mean_age;
    # boxplot_width = 0.1,
    # boxplot_nudge = 0.1,
    clouds=hist,
    center_boxplot=true,
    # cloud_width = 0.1,
    # gap = 0.5,
)
ylims!(ax, (0, nothing))
fig
# save plot
outputfile = joinpath(TMIPDIR, "extra", "ideal_age_rainclouds.png")
@info "Saving ideal age rainclouds as image file:\n  $(outputfile)"
save(outputfile, fig)



begin
    @show data = let
        # dataGBL[:, :basin] .= 1
        dataPAC[:, :basin] .= 1
        dataIND[:, :basin] .= 2
        dataATL[:, :basin] .= 3
        # vcat(dataGBL, dataATL, dataPAC, dataIND)
        vcat(dataATL, dataPAC, dataIND)
    end

    # Remove ACCESS1-0 r2i1p1 from the data because age = 0???
    delete!(data, findall((data.model .== "ACCESS1-0") .& (data.member .== "r2i1p1")))

    fig = Figure(size=(800, 300))
    emptyax = Axis(fig[1, 1];
        # ylabel="Basin",
        # yticks = (1:length(basin_strs) + 1, ["Global", basin_strs...]),
        yticks=(1:length(basin_strs), collect(basin_strs)),
        xticksvisible=false, yticksvisible=false,
        xticklabelsvisible=false,
        # xticksize = 0, yticksize = 0,
        xgridvisible = false, ygridvisible = false,
        leftspinevisible = false, rightspinevisible = false,
        topspinevisible = false, bottomspinevisible = false,
        # yticks=(0.5:1:3.5, fill("", 4)),
        # xlabel="Basin mean of ideal mean age (yr)",
    )
    ax = Axis(fig[1, 1];
        # ylabel="Basin",
        # yticks = (1:length(basin_strs) + 1, ["Global", basin_strs...]),
        # yticks=(1:length(basin_strs), collect(basin_strs)),
        yticks=(0.5:1:3.5, fill("", 4)),
        yticksvisible=false,
        xticks=0:100:1500, xtrimspine=true,
        xlabel="Basin mean of ideal mean age (yr)",
        xgridvisible = true, ygridvisible = false,
        leftspinevisible = false, rightspinevisible = false,
        topspinevisible = false, bottomspinevisible = true,
    )
    linkaxes!(emptyax, ax)
    model_index(model) = findfirst(isequal(model), models)
    cmap = cgrad(:seaborn_colorblind, alpha=1, categorical=true)
    dodge_gap = 0.1
    gap = 0.2
    # think white layer to spearate basins
    [hspan!(ax, y - gap/2, y + gap/2, color=:white) for y in 0.5:1:3.5]
    boxplot!(ax, data.basin, data.mean_age;
        dodge = model_index.(data.model),
        orientation = :horizontal,
        label=data.model,
        markersize=4,
        dodge_gap,
        gap,
        # color = (:blue, 0.2),
        color=cmap[model_index.(data.model)],
        # colormap = cgrad(:Archambault, alpha = 0.3),
    )
    # rainclouds!(ax, model_index.(dataGBL.model), dataGBL.mean_age;
    #     # boxplot_width = 0.1,
    #     # boxplot_nudge = 0.1,
    #     clouds = hist,
    #     center_boxplot = true,
    #     # cloud_width = 0.1,
    #     # gap = 0.5,
    # )
    # beeswarm!(ax, model_index.(dataGBL.model), dataGBL.mean_age;
    #     markersize = 5,
    #     orientation = :horizontal,
    #     color = (:black, 0.3),
    # )
    ylims!(ax, (0.5, 3.5))
    xlims!(ax, (10, nothing))
    n_dodge = length(models)
    @show dodge_width = Makie.scale_width(dodge_gap, n_dodge)
    width = 1 - gap
    for model in models
        # Only label the Indian
        # Pacific index is 1
        x = minimum(data.mean_age[(data.basin .== 1) .& (data.model .== model)]) - 20
        # x = maximum(data.mean_age[(data.basin .== 1)]) + 20
        @show shift = Makie.shift_dodge(model_index.(model), dodge_width, dodge_gap)
        y = 1 + width * shift
        # y = 3 + Makie.shift_dodge(model_index.(model), dodge_width, 0.0\)
        text!(ax, x, y; text=model, align=(:right, :center), color=cmap[model_index.(model)])
    end

    fig
end
# save plot
outputfile = joinpath(TMIPDIR, "extra", "ideal_age_rainclouds.png")
@info "Saving ideal age rainclouds as image file:\n  $(outputfile)"
save(outputfile, fig)























