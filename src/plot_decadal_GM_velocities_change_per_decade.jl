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
using Unitful: s, yr, m, kg, K, °C
# Must run twice v1.11.1 https://github.com/JuliaLang/julia/issues/56216
try
    using CairoMakie
catch
    using CairoMakie
end
using GeoMakie
using PairPlots
using Interpolations
using OceanBasins
using NaNStatistics
using Format
using GibbsSeaWater
using Statistics
using Makie.StructArrays
using ProgressMeter


include("plotting_functions.jl")

vars = ("tx_trans_gm", "ty_trans_gm", "tx_trans_submeso", "ty_trans_submeso")
# vars = ("tx_trans_gm", "ty_trans_gm",)
Nvars = length(vars)

datadir = "/g/data/xv83/TMIP/data"

modeldir = joinpath(datadir, "ACCESS-ESM1-5")

# experiments = ["historical", "ssp370"]
# exp_prefix = ["HI", "SSP-370"]
# experiments = ["historical"]
# exp_prefix = ["HI"]
experiments = ["ssp370"]
exp_prefix = ["SSP-370"]

# members = 1:40
# members = 1:20
# members = 3:20
members = 21:40
# members = [1; 3:8]
# members = [1]

colormap = cgrad([:white; collect(cgrad(:managua))])


for (experiment, exp_prefix) in zip(experiments, exp_prefix)

    println("$experiment")
    expdir = joinpath(modeldir, experiment)
    decades = (experiment == "historical") ? (1850:10:2010) : (2010:10:2090)
    Ndecades = length(decades)

    for member in members

        println("  $member")

        CSIRO_member = "$exp_prefix-$(format(member + 4, width=2, zeropadding=true))"
        inputdir = joinpath(expdir, CSIRO_member)

        fig = Figure(size = (200Nvars, 200(Ndecades - 1)))

        @time for (icol, var) in enumerate(vars)

            decade_str = "$(first(decades))s"
            fname = joinpath(inputdir, "month_$(var)_$(decade_str).nc")
            !isfile(fname) && continue
            println("    $fname")
            previous_x = nanmean(replace(ncread(fname, var) |> Array, -1f20 => NaN), dims = 4) |> vec

            for (irow, decade) in enumerate(decades[2:end])

                decade_str = "$(decade)s"
                fname = joinpath(inputdir, "month_$(var)_$(decade_str).nc")
                x = nanmean(replace(ncread(fname, var) |> Array, -1f20 => NaN), dims = 4) |> vec

                ax = Axis(fig[irow, icol], aspect = 1)
                # scatter!(aprevious_x, x, x)

                points = StructArray{Point2f}((previous_x, x))
                ds = datashader!(ax, points; colormap, async = false)
                ablines!(ax, 0, 1, color = (:black, 0.1), linewidth = 10)
                vlines!(ax, 0, color = (:black, 0.1), linewidth = 10)
                hlines!(ax, 0, color = (:black, 0.1), linewidth = 10)

                hidedecorations!(ax)

                (icol == 1) && Label(fig[irow, 0]; text = decade_str, tellheight = false)

                previous_x = x

            end

            Label(fig[0, icol]; text = var, tellwidth = false)

        end

        # save plot
        outputfile = joinpath(inputdir, "decadal_GM_velocities_change.png")
        @info "Saving image file:\n  $(outputfile)"
        save(outputfile, fig)

    end

end





