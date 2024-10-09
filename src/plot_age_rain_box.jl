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

models = ["ACCESS1-3", "ACCESS-ESM1-5", "ACCESS-CM2"]

CMIP5_models = ["ACCESS1-3"]

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



basin_keys = (:ATL, :PAC, :IND)
basin_strs = ("Atlantic", "Pacific", "Indian")
basin_functions = (isatlantic, ispacific, isindian)

dataGBL = DataFrame(model=String[], mean_age=Float64[])
dataATL = DataFrame(model=String[], mean_age=Float64[])
dataPAC = DataFrame(model=String[], mean_age=Float64[])
dataIND = DataFrame(model=String[], mean_age=Float64[])
# Fetch the ideal mean age + volume and compute the global mean of the ideal mean age
for model in models

    CMIP_version = model ∈ CMIP5_models ? "CMIP5" : "CMIP6"

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
        if model ∈ ("ACCESS-ESM1-5", "ACCESS-CM2")
            lon_vertices = readcubedata(volcello_ds.lon_verticies) # xmip issue: https://github.com/jbusecke/xMIP/issues/369
            lat_vertices = readcubedata(volcello_ds.lat_verticies) # xmip issue: https://github.com/jbusecke/xMIP/issues/369
        elseif model ∈ ("ACCESS1-3", "ACCESS1-0")
            lon_vertices = readcubedata(volcello_ds.lon_vertices) # no xmip so default key
            lat_vertices = readcubedata(volcello_ds.lat_vertices) # no xmip so default key
        else
            error("Need to hardcode the lon_vertices and lat_vertices keys for $model in the tests")
        end

        modelgrid = makemodelgrid(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
        (; lon, lat) = modelgrid

        # Save the global mean of the ideal mean age
        _FillValue = get(volcello.properties, "_FillValue", NaN)
        volecllo = volcello |> Array
        volcello = replace(volcello, missing => NaN, _FillValue => NaN)
        mean_age = nansum(age .* volcello) ./ nansum(volcello)
        push!(dataGBL, (model, mean_age))
        # println(" $mean_age $(nansum(age)) $(nansum(volcello))")

        basin_values = (reshape(f(lat[:], lon[:], OCEANS), size(lat)) for f in basin_functions)
        basins = (; (basin_keys .=> basin_values)...)

        ATL_mean_age = nansum(age .* volcello .* basins.ATL) ./ nansum(volcello .* basins.ATL)
        PAC_mean_age = nansum(age .* volcello .* basins.PAC) ./ nansum(volcello .* basins.PAC)
        IND_mean_age = nansum(age .* volcello .* basins.IND) ./ nansum(volcello .* basins.IND)
        push!(dataATL, (model, ATL_mean_age))
        push!(dataPAC, (model, PAC_mean_age))
        push!(dataIND, (model, IND_mean_age))

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
        dataATL[:, :basin] .= 1
        dataPAC[:, :basin] .= 2
        dataIND[:, :basin] .= 3
        # vcat(dataGBL, dataATL, dataPAC, dataIND)
        vcat(dataATL, dataPAC, dataIND)
    end

    fig = Figure(size=(800, 300))
    ax = Axis(fig[1, 1];
        ylabel="Basin",
        # yticks = (1:length(basin_strs) + 1, ["Global", basin_strs...]),
        yticks=(1:length(basin_strs), collect(basin_strs)),
        xlabel="Global mean ideal mean age (yr)",
    )
    model_index(model) = findfirst(isequal(model), models)
    cmap = cgrad(:seaborn_colorblind, alpha=1, categorical=true)
    dodge_gap = 0.1
    gap = 0.2
    boxplot!(ax, data.basin, data.mean_age;
        dodge = model_index.(data.model),
        orientation = :horizontal,
        label=data.model,
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
    xlims!(ax, (0, nothing))
    n_dodge = length(models)
    @show dodge_width = Makie.scale_width(dodge_gap, n_dodge)
    width = 1 - gap
    for model in models
        # Only label the Indian
        x = minimum(dataIND.mean_age[model.==dataATL.model]) - 20
        # Indian index is 3
        @show shift = Makie.shift_dodge(model_index.(model), dodge_width, dodge_gap)
        y = 3 + width * shift
        # y = 3 + Makie.shift_dodge(model_index.(model), dodge_width, 0.0\)
        text!(ax, x, y; text=model, align=(:right, :center), color=cmap[model_index.(model)])
    end
    fig
end
# save plot
outputfile = joinpath(TMIPDIR, "extra", "ideal_age_rainclouds.png")
@info "Saving ideal age rainclouds as image file:\n  $(outputfile)"
save(outputfile, fig)























