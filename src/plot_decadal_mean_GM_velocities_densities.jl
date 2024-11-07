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
using CairoMakie
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

# members = 1:40
members = [1; 3:8]
# members = [1]

# decades = 1850:10:2010
# decades = 1850:10:1870
decades = 1990:10:1990
Ndecades = length(decades)


for member in members

    CMIP6_member = "r$(member)i1p1f1"

    inputdir = "/scratch/xv83/TMIP/data/ACCESS-ESM1-5/historical/$CMIP6_member/Jan1990-Dec1999/"

    fig = Figure(size = (200Nvars, 200Ndecades))

    for (icol, var) in enumerate(vars)

        for (irow, decade) in enumerate(decades)

            decade_str = "$(decade)s"
            fname = joinpath(inputdir, "$(var).nc")
            x3D = ncread(fname, var) |> Array
            x3Dnans = isnan.(x3D)
            x3Dzeros = iszero.(x3D)
            x3Dfillvalues = x3D .== -1f20

            x = x3D[.!x3Dnans .& .!x3Dzeros .& .!x3Dfillvalues]

            ax = Axis(fig[irow, icol], aspect = 1)
            density!(ax, x)

            (icol > 1) && linkaxes!(ax, content(fig[1, icol]))

            hidedecorations!(ax)

            (icol == 1) && Label(fig[irow, 0]; text = decade_str, tellheight = false)

            text!(ax, 0, 1; text="$(sum(x3Dnans)) nans", align = (:left, :top), space = :relative)
            text!(ax, 0, 0.85; text="$(sum(x3Dzeros)) zeros", align = (:left, :top), space = :relative)
            text!(ax, 0, 0.7; text="$(sum(x3Dfillvalues)) fillvalues", align = (:left, :top), space = :relative)

        end

        Label(fig[0, icol]; text = var, tellwidth = false)

    end

    # save plot
    outputfile = joinpath(inputdir, "decadal_mean_GM_velocities_densities.png")
    @info "Saving image file:\n  $(outputfile)"
    save(outputfile, fig)

end








