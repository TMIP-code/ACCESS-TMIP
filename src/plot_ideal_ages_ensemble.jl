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
using CairoMakie
using GeoMakie
using Interpolations
using OceanBasins
using NaNStatistics
using CairoMakie.Colors

include("plotting_functions.jl")

models = ("ACCESS-ESM1-5", "ACCESS-CM2", "ACCESS1-3")
CMIP_versions = ("CMIP6", "CMIP6", "CMIP5")
colors = cgrad([cgrad(:tableau_colorblind, categorical = true)[1:2]; :black])
# colors = cgrad(:okabe_ito, categorical = true)[[3,2,1]]
colors = [(c, 0.5) for c in colors]

experiment = "historical"

time_window = "Jan1990-Dec1999"

parse_member(member, member_regex) = parse.(Int, match(member_regex, member).captures)

Γinyr1Ds = Dict()
zts = Dict()
agessc1Ds = Dict()

for (model, CMIP_version) in zip(models, CMIP_versions)

    # Gadi directory for input files
    inputdirfun(member) = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)"

    requiredvariables = ["mlotst", "volcello", "areacello", "ideal_mean_age", "agessc"]
    hasrequireddata(member, variable_name) = isfile(joinpath(inputdirfun(member), "$variable_name.nc"))
    hasrequireddata(member) = all(variable_name -> hasrequireddata(member, variable_name), requiredvariables)

    @info "Processing $model $experiment $(time_window)"

    Γinyr1Ds[model] = Dict()
    zts[model] = Dict()
    agessc1Ds[model] = Dict()

    # find all members for which the inputdir contains umo.nc, vmo.nc, mlotst.nc, volcello.nc, and areacello.nc
    members = readdir("/scratch/xv83/TMIP/data/$model/$experiment")

    # sort members by r, i, p[, f]

    member_regex = CMIP_version == "CMIP6" ? r"r(\d+)i(\d+)p(\d+)f(\d+)" : r"r(\d+)i(\d+)p(\d+)"
    members = sort(members, by = x -> parse_member(x, member_regex))
    dataavailability = DataFrame(
        :member => members,
        :has_it_all => hasrequireddata.(members),
        [Symbol(var) => [hasrequireddata(member, var) for member in members] for var in requiredvariables]...,
    )
    show(dataavailability; allrows = true)
    println()



    for member in members[dataavailability.has_it_all]
    # member = last(members)

        @info "  member $member"

        Γinyr1Ds[model][member] = Dict()
        agessc1Ds[model][member] = Dict()

        inputdir = inputdirfun(member)

        # Load umo, vmo, mlotst, volcello, and areacello
        mlotst_ds = open_dataset(joinpath(inputdir, "mlotst.nc"))
        volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
        areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))
        agessc_ds = open_dataset(joinpath(inputdir, "agessc.nc"))

        mlotst = mlotst_ds["mlotst"] |> Array{Float64}

        # Make makegridmetrics
        gridmetrics = makegridmetrics(; areacello_ds, volcello_ds, mlotst_ds)

        # Make indices
        indices = makeindices(gridmetrics.v3D)

        # unpack model grid
        (; lon, lat, zt, v3D,) = gridmetrics

        zts[model][member] = zt

        Γin_ds = open_dataset(joinpath(inputdir, "ideal_mean_age.nc"))

        Γinyr3D = Γin_ds["age"] |> Array{Float64}
        agessc3D = agessc_ds["agessc"] |> Array{Float64}

        # save horizontal average
        basin_keys = (:ATL, :PAC, :IND)
        basin_strs = ("Atlantic", "Pacific", "Indian")
        basin_functions = (isatlantic, ispacific, isindian)
        basin_values = (reshape(f(lat[:], lon[:], OCEANS), size(lat)) for f in basin_functions)
        basins = (; (basin_keys .=> basin_values)...)

        for (basin_key, mask) in pairs(basins)
            Γinyr1Ds[model][member][basin_key] = horizontalaverage(Γinyr3D, gridmetrics; mask)
            agessc1Ds[model][member][basin_key] = horizontalaverage(agessc3D, gridmetrics; mask)
        end

    end
end


# # Plot Γ↓ basin profiles
# basin_keys = (:ATL, :PAC, :IND)
# basin_strs = ("Atlantic", "Pacific", "Indian")

fig = Figure(size = (1000, 1000), fontsize = 18)
axs = Array{Any,2}(undef, (2, 3))
plts = Array{Any,2}(undef, (2, 3))

zlim = (6000, 0)

for (icol, (basin_key, basin_str)) in enumerate(zip(basin_keys, basin_strs))

    for (irow, x3Ds) in enumerate((Γinyr1Ds, agessc1Ds))

        ax = Axis(fig[irow, icol],
            xgridvisible = true, ygridvisible = true,
            ylabel = "depth (m)",
            xlabel = "ideal mean age (yr)",
        )

        for (imodel, (model, color)) in enumerate(zip(models, colors))

            for member in keys(x3Ds[model])

                x1D = x3Ds[model][member][basin_key]
                zt = zts[model][member]

                lines!(ax, x1D, zt; linewidth = 3, color, label = model)

            end

        end

        ylims!(ax, zlim)

        hidexdecorations!(ax,
            label = irow < 2, ticklabels = irow < 2,
            ticks = irow < 2, grid = false)
        hideydecorations!(ax,
            label = icol > 1, ticklabels = icol > 1,
            ticks = icol > 1, grid = false)


        axs[irow, icol] = ax

        if icol == 1 && irow == 1
            axislegend(ax, merge = true, unique = true, framevisible = false)
        end

        if icol == 1
            if irow == 1
                text = "Steady state from 1990s transport matrix"
            else
                text = "Mean 1990s agessc variable"
            end
            Label(fig[irow, 0]; text, tellheight = false, rotation = π/2)
        end

    end

    linkxaxes!(axs[:, icol]...)

    Label(fig[0, icol], basin_str, tellwidth = false)

end

Label(fig[-1, 1:3], "Ideal mean age: Steady-state from transport matrices vs agessc (1990s)", tellwidth = false, fontsize = 24)

# save plot
inputdir = joinpath("/scratch/xv83/TMIP/extra/")
mkpath(inputdir)
outputfile = joinpath(inputdir, "ideal_age_all_basin_profiles.png")
@info "Saving ideal mean age profiles as image file:\n  $(outputfile)"
save(outputfile, fig)



