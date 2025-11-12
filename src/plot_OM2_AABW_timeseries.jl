# qsub -I -P y99 -q express -l mem=47GB -l storage=scratch/gh0+scratch/xv83+scratch/p66 -l walltime=01:00:00 -l ncpus=12

using Pkg
Pkg.activate(".")
Pkg.instantiate()

using OceanTransportMatrixBuilder
using NetCDF
using YAXArrays
using YAXArrayBase
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
using GeometryBasics
using GeometryOps
using LibGEOS
# using LaTeXStrings
using Format
using KernelDensity
using RollingWindowArrays
using Dates
using DataStructures

include("plotting_functions.jl")

#############
# Load data #
#############

files = OrderedDict(
    "1" => "/scratch/y99/TMIP/data/ACCESS-OM2-1/1deg_jra55_iaf_omip2_cycle6/psi_tot.nc", # 1 degree
    "025" => "/scratch/y99/TMIP/data/ACCESS-OM2-025/025deg_jra55_iaf_omip2_cycle6/psi_tot.nc", # 0.25 degree
    "01" => "/scratch/y99/TMIP/data/ACCESS-OM2-01/01deg_jra55v140_iaf_cycle4/psi_tot.nc", # 0.1 degree
)


psis = OrderedDict{String, Any}()
ρs = OrderedDict{String, Any}()
lats = OrderedDict{String, Any}()
times = OrderedDict{String, Any}()
AABWs = OrderedDict{String, Any}()
idxs = OrderedDict{String, Any}()
meanpsis = OrderedDict{String, Any}()
latselectors = OrderedDict(
    # "any lat" => Colon(),
    "lat ≤ 50°S" => Where(≤(-50)),
    "lat ≤ 40°S" => Where(≤(-40)),
    "50°S ≤ lat ≤ 40°S" => -50 .. -40,
    "lat ≈ 40°S" => Near(-40),
)
ρselectors = OrderedDict(
    # "any ρ" => Colon(),
    "ρ ≥ 1036 kg/m³" => Where(≥(1036.0)),
    "ρ ≥ 1036.8 kg/m³" => Where(≥(1036.8)),
)


for (res, file) in pairs(files)
    @info """Loading OM2-$res overturning streamfunction
      file: $file"""
    ds = open_dataset(file)
    psis[res] = set(readcubedata(ds.psi_tot), :time => Ti)
    meanpsis[res] = dropdims(mean(psis[res]; dims = Ti); dims = Ti)
    ρs[res] = ds.potrho
    lats[res] = ds.grid_yu_ocean
    times[res] = ds.time
    println("  Computing AABW transport metrics:")
    AABW = OrderedDict{String, Any}()
    idx = OrderedDict{String, Any}()
    psi = psis[res]
    (grid_yu_ocean, potrho) = otherdims(psi, (Ti,))
    for (latstr, latselector) in pairs(latselectors)
        println(" - $latstr:")
        AABW_sub = OrderedDict{String, Any}()
        idx_sub = OrderedDict{String, Any}()
        for (ρstr, ρselector) in pairs(ρselectors)
            # AABW_sub[ρstr] = minimum(psi[latselector(lats[res]), ρselector(ρs[res]), :], dims = (:grid_yu_ocean, :potrho)) |> vec
            psi_subdomain = psi[DimensionalData.rebuild(grid_yu_ocean, latselector), DimensionalData.rebuild(potrho, ρselector)]
            AABW_sub[ρstr] = minimum(psi_subdomain, dims = (grid_yu_ocean, potrho))
            idx_sub[ρstr] = argmin(psi_subdomain.data, dims = 1:(length(size(psi_subdomain)) - 1))
            println("   - $ρstr (size: $(size(AABW_sub[ρstr]))")
        end
        AABW[latstr] = AABW_sub
        idx[latstr] = idx_sub
    end
    AABWs[res] = AABW
    idxs[res] = idx
end

# colors = Makie.wong_colors()
colors = cgrad(:tab10, categorical = true)

fig = Figure(size = (1500, 750))
g = fig[1,1] = GridLayout(length(ρselectors), length(latselectors))
axs = Array{Any, 2}(undef, (length(ρselectors), length(latselectors)))

for (icol, (latstr, latselector)) in enumerate(pairs(latselectors))
    for (irow, (ρstr, ρselector)) in enumerate(pairs(ρselectors))
        ax = axs[irow, icol] = Axis(
            g[irow, icol];
            ylabel = "AABW transport (Sv)",
            # limits = (DateTime(1900), DateTime(2030), nothing, nothing),
            # xticks = DateTime.(1900:20:2030),
            # xticks = 1900:20:2030,
        )
        # ax, plt =
        # ax, plt = CairoMakie.lines(fig[irow, icol], [DateTime(1990)], [10.0]; axis=(; xticks = DateTime.(1900:20:2030), ylabel = "AABW transport (Sv)"))
        plts = Vector{Any}(undef, length(files))
        for (ires, res) in enumerate(keys(files))
            AABW = AABWs[res][latstr][ρstr].data |> vec
            time = times[res].val
            # rolling averages
            AABW_decade = fill(NaN, size(AABW))
            AABW_rolling = mean.(rolling(AABW, 120; center = true))
            offset = AABW_rolling.offsets[1]
            AABW_decade[offset .+ (0:(length(AABW_rolling.parent) - 1))] = AABW_rolling.parent
            AABW_year = fill(NaN, size(AABW))
            AABW_rolling = mean.(rolling(AABW, 12; center = true))
            offset = AABW_rolling.offsets[1]
            AABW_year[offset .+ (0:(length(AABW_rolling.parent) - 1))] = AABW_rolling.parent
            plts[ires] = lines!(ax, DateTime.(time), -AABW / 1.0e9; color = colors[ires], alpha = 0.1) # Sv
            lines!(ax, DateTime.(time), -AABW_year / 1.0e9; color = colors[ires], alpha = 0.3) # Sv
            lines!(ax, DateTime.(time), -AABW_decade / 1.0e9; color = colors[ires], label = res) # Sv
        end
        (icol == irow == 1) && axislegend(ax, position = :lt)
        ax.xticks = DateTimeTicks(6)
        myhidexdecorations!(ax, irow < length(ρselectors))
        # myhideydecorations!(ax, icol > 1)

        (icol == 1) && Label(g[irow, 0], ρstr, tellheight = false, rotation = π / 2)

        # Add inset plot of mean overturning streamfunction + AABW box for minimum
        res = "025"
        ρmin, ρmax = extrema(ρs["025"].val)
        ρmin -= eps(ρmin)
        # yscale = ReversibleScale(ρ -> exp(ρ - ρmin), x -> log(x) + ρmin, limits = (ρmin, ρmax))
        yscale = ReversibleScale(ρ -> (ρ - ρmin)^4, x -> x^(1 / 4) + ρmin, limits = (ρmin, ρmax))
        # yscale = ReversibleScale(ρ -> (ρ - ρmin)^2, x -> sqrt(x) + ρmin, limits = (ρmin, ρmax))
        ax_inset = Axis(
            g[irow, icol];
            width = Relative(0.4),
            height = Relative(0.2),
            halign = 0.9,
            valign = 0.9,
            yreversed = true,
            limits = (extrema(lats[res].val)..., ρmin, ρmax),
            yscale = yscale,
            yticks = [1030, 1033, 1035, 1036, 1037],
        )
        contourf!(
            ax_inset, lats[res].val, ρs[res].val, meanpsis[res].data / 1.0e9; # Sv
            levels = -24:2:24,
            colormap = :delta,
            extendlow = :auto,
            extendhigh = :auto,
        )
        (grid_yu_ocean, potrho) = dims(meanpsis[res])
        latbox = grid_yu_ocean[collect(extrema(DimensionalData.selectindices(grid_yu_ocean, latselector)))].val
        ρbox = potrho[collect(extrema(DimensionalData.selectindices(potrho, ρselector)))].val
        lines!(ax_inset, latbox[[1, 2, 2, 1, 1]], ρbox[[1, 1, 2, 2, 1]]; color = :red)
    end
    Label(g[0, icol], latstr, tellwidth = false)
end

# linkaxes!(axs[:]...)
rowgap!(g, 0.0)
colgap!(g, 0.0)

fig

imagefile = "/scratch/y99/TMIP/data/AABW_timeseries.png"
save(imagefile, fig)
@info "Saved AABW timeseries plot to $imagefile"
imagefile = "/scratch/y99/TMIP/data/AABW_timeseries.pdf"
save(imagefile, fig)
@info "Saved AABW timeseries plot to $imagefile"
