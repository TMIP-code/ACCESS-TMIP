# # qsub -I -P y99 -q express -l mem=47GB -l storage=scratch/gh0+scratch/xv83+scratch/p66 -l walltime=01:00:00 -l ncpus=12
# # qsub -I -P y99 -l mem=47GB -l storage=scratch/gh0+scratch/xv83+scratch/p66 -l walltime=01:00:00 -l ncpus=12


# # TODO: Plot MOC overtutning for all files in a separate script
# # because I keep seeing some weird stuff (different \rho ranges, AABW min locations, etc.)
# # Probably a 3x6 panels for 3 resolutios and cycles* (4 + 2 for 0.1deg...)


# using Pkg
# Pkg.activate(".")
# Pkg.instantiate()

# using OceanTransportMatrixBuilder
# using NetCDF
# using YAXArrays
# using YAXArrayBase
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
# using Interpolations
# using OceanBasins
# using Statistics
# using NaNStatistics
# using StatsBase
# using FileIO
# using Contour
# using GeometryBasics
# using GeometryOps
# using LibGEOS
# using DataFrames
# # using LaTeXStrings
# using Format
# using KernelDensity
# using RollingWindowArrays
# using Dates
# using DataStructures

# include("plotting_functions.jl")

# #############
# # Load data #
# #############


# files = [
#     # "/scratch/y99/TMIP/data/ACCESS-OM2-1/1deg_jra55_iaf_omip2_cycle1/psi_tot.nc"
#     # "/scratch/y99/TMIP/data/ACCESS-OM2-1/1deg_jra55_iaf_omip2_cycle2/psi_tot.nc"
#     # "/scratch/y99/TMIP/data/ACCESS-OM2-1/1deg_jra55_iaf_omip2_cycle3/psi_tot.nc"
#     # "/scratch/y99/TMIP/data/ACCESS-OM2-1/1deg_jra55_iaf_omip2_cycle4/psi_tot.nc"
#     # "/scratch/y99/TMIP/data/ACCESS-OM2-1/1deg_jra55_iaf_omip2_cycle5/psi_tot.nc"
#     # "/scratch/y99/TMIP/data/ACCESS-OM2-1/1deg_jra55_iaf_omip2_cycle6/psi_tot.nc"
#     # "/scratch/y99/TMIP/data/ACCESS-OM2-025/025deg_jra55_iaf_omip2_cycle1/psi_tot.nc"
#     # "/scratch/y99/TMIP/data/ACCESS-OM2-025/025deg_jra55_iaf_omip2_cycle2/psi_tot.nc"
#     # "/scratch/y99/TMIP/data/ACCESS-OM2-025/025deg_jra55_iaf_omip2_cycle3/psi_tot.nc"
#     # "/scratch/y99/TMIP/data/ACCESS-OM2-025/025deg_jra55_iaf_omip2_cycle4/psi_tot.nc"
#     # "/scratch/y99/TMIP/data/ACCESS-OM2-025/025deg_jra55_iaf_omip2_cycle5/psi_tot.nc"
#     # "/scratch/y99/TMIP/data/ACCESS-OM2-025/025deg_jra55_iaf_omip2_cycle6/psi_tot.nc"
#     "/scratch/y99/TMIP/data/ACCESS-OM2-01/01deg_jra55v140_iaf/psi_tot.nc"
#     "/scratch/y99/TMIP/data/ACCESS-OM2-01/01deg_jra55v140_iaf_cycle2/psi_tot.nc"
#     "/scratch/y99/TMIP/data/ACCESS-OM2-01/01deg_jra55v140_iaf_cycle3/psi_tot.nc"
#     "/scratch/y99/TMIP/data/ACCESS-OM2-01/01deg_jra55v140_iaf_cycle4/psi_tot.nc"
#     # "/scratch/y99/TMIP/data/ACCESS-OM2-01/01deg_jra55v140_iaf_cycle4_jra55v150_extension/psi_tot.nc"
#     # "/scratch/y99/TMIP/data/ACCESS-OM2-01/01deg_jra55v150_iaf_cycle1/psi_tot.nc"
# ]
# model_str(file) = split(file, "/")[6]
# resolution_str(file) = split(model_str(file), "-")[3]
# jra55version_str(file) = occursin("jra55v150", file) ? "jra55v150" : "jra55v140" # CHECK that it's v1.4 when unspecified
# cycle_match(file) = match(r"cycle(\d+)", file)
# cycle_str(file) = (c = cycle_match(file); isnothing(c) ? 1 : parse(Int, c.captures[1]))
# files_df = DataFrame()
# id_str(file) = split(file, "/")[7]
# for file in files
#     model = model_str(file)
#     resolution = resolution_str(file)
#     cycle = cycle_str(file)
#     jra_version = jra55version_str(file)
#     id = id_str(file)
#     push!(files_df, (; id, file, model, resolution, cycle, jra_version))
# end
# @show(files_df)

# psis = OrderedDict{String, Any}()
# ρs = OrderedDict{String, Any}()
# lats = OrderedDict{String, Any}()
# times = OrderedDict{String, Any}()
# meanpsis = OrderedDict{String, Any}()


# for row in eachrow(files_df)
#     (; id, file, model, resolution, cycle, jra_version) = row
#     @info "Loading $model $jra_version cycle $cycle overturning streamfunction"
#     ds = open_dataset(file)
#     psis[id] = set(readcubedata(ds.psi_tot), :time => Ti)
#     meanpsis[id] = dropdims(mean(psis[id]; dims = Ti); dims = Ti)
#     ρs[id] = ds.potrho
#     lats[id] = ds.grid_yu_ocean
#     times[id] = ds.time
# end

# ########
# # Plot #
# ########

# colors = Makie.wong_colors()
resolutions = unique(files_df.resolution)
resolution° = Dict(
    "1" => "1°",
    "025" => "0.25°",
    "01" => "0.1°",
)
cycles = 1:4
Nrows = length(resolutions)
Ncols = length(cycles)

fig = Figure(;
    size = (400 * Ncols, 450 * Nrows),
    fonts = (; regular = "Dejavu"),
)
# g = fig[1, 1] = GridLayout(Ncols, Nrows)

cfs = Array{Any, 2}(undef, Nrows, Ncols)


ρmin, ρmax = extrema(ρs["01deg_jra55v140_iaf"].val)
latmin, latmax = extrema(lats["01deg_jra55v140_iaf"].val)
ρmin -= eps(ρmin)
# yscale = ReversibleScale(ρ -> exp(ρ - ρmin), x -> log(x) + ρmin, limits = (ρmin, ρmax))
yscale = ReversibleScale(ρ -> (ρ - ρmin)^4, x -> x^(1 / 4) + ρmin, limits = (ρmin, ρmax))

levels = -24:2:24
colormap = :delta
limits = (-90, +90, 1035, 1037.3)

for (ijrow, row) in enumerate(eachrow(files_df))

    (; id, file, model, resolution, cycle, jra_version) = row
    icol, irow = CartesianIndices((Ncols, Nrows))[ijrow].I

    # ((icol == 1) && (irow == 1)) || ((icol == 2) && (irow == 3)) || continue # for testing

    ax = Axis(
        fig[irow, icol];
        xlabel = "Latitude",
        ylabel = "Density (kg/m³)",
        yreversed = true,
        # limits = (extrema(lats[id].val)..., ρmin, ρmax),
        limits,
        yscale = yscale,
        # xticks = -80:20:80,
        xticks = -60:30:60,
        # yticks = round.([ρ for ρ in ρbox if ρ ≉ ρmin], digits = 1),
        # xticks = round.(Int, [l for l in latbox if l ≉ latmin]),
        xtickformat = ytickformat,
        # title = id,
    )

    lat = lats[id].val

    ρ = ρs[id].val

    # Skim some of weirdness at large densities
    iρ = ρ .<= limits[4]
    ρ = ρ[iρ]
    ψ = meanpsis[id].data[:, iρ] / 1.0e9 # convert to Sv

    co = cfs[irow, icol] = contourf!(
        ax, lat, ρ, ψ;
        levels,
        colormap,
        extendlow = :auto,
        extendhigh = :auto,
    )
    translate!(co, 0, 0, -100)

    text!(ax, 0, 0; text = id, align = (:left, :bottom), space = :relative, offset = (5, 3))

    # contour!(ax, lat, ρ, ψ; levels = levels[(end ÷ 2 + 1):end], color = :black, linewidth = 0.5)
    # contour!(ax, lat, ρ, ψ; levels = levels[1:(end ÷ 2)], color = :black, linewidth = 0.5, linestyle = :dash)

    myhidexdecorations!(ax, irow < Nrows)
    myhideydecorations!(ax, icol > 1)

    (icol == 1) && Label(fig[irow, 0], "$(resolution°[resolution]) ACCESS-OM2", tellheight = false, rotation = π / 2, font = :bold)
    (irow == 1) && Label(fig[0, icol], rich("Cycle $cycle", cycle > 4 ? superscript("*") : ""), tellwidth = false, font = :bold)

    # add a red box over those 0.1deg runs that are not part of the main 6 cycles...
    # (resolution == "01") && (icol > 4) && Box(fig[irow, icol]; color = (Pattern('/'), 0.3), strokecolor = :transparent)
    lp = Makie.LinePattern(; direction = Makie.Vec2f(1, 1), width = 1, tilesize = (12, 12), linecolor = (:gray, 0.5), backgroundcolor = (:white, 0.5))
    (resolution == "01") && (icol > 4) && Box(fig[irow, icol]; color = lp, strokecolor = :transparent)

end

linkxaxes!([x for x in fig.content if x isa Axis]...)

cb = Colorbar(
    fig[Nrows + 1, 1:Ncols], cfs[1, 1];
    vertical = false, flipaxis = false,
    label = "Oveturning Streamfunction (Sv)",
    tickformat = divergingcbarticklabelformat,
)
cb.width = Relative(1 / 4)

# Box(fig[Nrows, 1:Ncols], color = (:green, 0.2))
# Box(fig[Nrows + 1, 1:Ncols], color = (:red, 0.2))

fig

imagefile = "/scratch/y99/TMIP/data/OM2-01_MOC_4cycles.png"
save(imagefile, fig)
@info "Saved AABW timeseries plot to $imagefile"
# imagefile = "/scratch/y99/TMIP/data/OM2-01_MOC_4cycles.pdf"
# save(imagefile, fig)
# @info "Saved AABW timeseries plot to $imagefile"
