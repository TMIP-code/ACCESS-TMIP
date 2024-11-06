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

# vars = ("tx_trans_gm", "ty_trans_gm", "tx_trans_submeso", "ty_trans_submeso")
vars = ("tx_trans_gm", "ty_trans_gm",)
Nvars = length(vars)

# members = 1:40
# members = [1; 3:8]
members = [1]

# decades = 1850:10:2010
decades = 1850:10:1870
Ndecades = length(decades)


@showprogress desc = "members" for member in members

    CSIRO_member = "HI-$(format(member + 4, width = 2, zeropadding = true))"

    inputdir = "/g/data/xv83/TMIP/data/ACCESS-ESM1-5/$CSIRO_member/"

    fig = Figure(size = (200Nvars, 200(Ndecades - 1)))

    @showprogress desc = "$member vars" for (icol, var) in enumerate(vars)

        @showprogress desc = "$var decades" for (irow, decade) in enumerate(decades)

            decade_str = "$(decade)s"
            fname = joinpath(inputdir, "month_$(var)_$(decade_str).nc")
            x4D = ncread(fname, var) |> Array
            x = x4D[x4D .≠ -1f20]

            ax = Axis(fig[irow, icol], aspect = 1)
            density!(ax, x)

            (icol > 1) && linkaxes!(ax, content(fig[1, icol]))
            hidedecorations!(ax)

            (icol == 1) && Label(fig[irow, 0]; text = decade_str, tellheight = false)

        end

        Label(fig[0, icol]; text = var, tellwidth = false)

    end

    # save plot
    outputfile = joinpath(inputdir, "decadal_GM_velocities_densities.png")
    @info "Saving image file:\n  $(outputfile)"
    save(outputfile, fig)

end








