# qsub -I -P y99 -q express -l mem=47GB -l storage=scratch/gh0+scratch/xv83+scratch/p66 -l walltime=01:00:00 -l ncpus=12
# qsub -I -P y99 -l mem=47GB -l storage=scratch/gh0+scratch/xv83+scratch/p66 -l walltime=01:00:00 -l ncpus=12

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
using DataFrames
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


files = [
    "/scratch/y99/TMIP/data/ACCESS-OM2-1/1deg_jra55_iaf_omip2_cycle1/psi_tot.nc"
    "/scratch/y99/TMIP/data/ACCESS-OM2-1/1deg_jra55_iaf_omip2_cycle2/psi_tot.nc"
    "/scratch/y99/TMIP/data/ACCESS-OM2-1/1deg_jra55_iaf_omip2_cycle3/psi_tot.nc"
    "/scratch/y99/TMIP/data/ACCESS-OM2-1/1deg_jra55_iaf_omip2_cycle4/psi_tot.nc"
    "/scratch/y99/TMIP/data/ACCESS-OM2-1/1deg_jra55_iaf_omip2_cycle5/psi_tot.nc"
    "/scratch/y99/TMIP/data/ACCESS-OM2-1/1deg_jra55_iaf_omip2_cycle6/psi_tot.nc"
    "/scratch/y99/TMIP/data/ACCESS-OM2-025/025deg_jra55_iaf_omip2_cycle1/psi_tot.nc"
    "/scratch/y99/TMIP/data/ACCESS-OM2-025/025deg_jra55_iaf_omip2_cycle2/psi_tot.nc"
    "/scratch/y99/TMIP/data/ACCESS-OM2-025/025deg_jra55_iaf_omip2_cycle3/psi_tot.nc"
    "/scratch/y99/TMIP/data/ACCESS-OM2-025/025deg_jra55_iaf_omip2_cycle4/psi_tot.nc"
    "/scratch/y99/TMIP/data/ACCESS-OM2-025/025deg_jra55_iaf_omip2_cycle5/psi_tot.nc"
    "/scratch/y99/TMIP/data/ACCESS-OM2-025/025deg_jra55_iaf_omip2_cycle6/psi_tot.nc"
    "/scratch/y99/TMIP/data/ACCESS-OM2-01/01deg_jra55v140_iaf/psi_tot.nc"
    "/scratch/y99/TMIP/data/ACCESS-OM2-01/01deg_jra55v140_iaf_cycle2/psi_tot.nc"
    "/scratch/y99/TMIP/data/ACCESS-OM2-01/01deg_jra55v140_iaf_cycle3/psi_tot.nc"
    "/scratch/y99/TMIP/data/ACCESS-OM2-01/01deg_jra55v140_iaf_cycle4/psi_tot.nc"
    "/scratch/y99/TMIP/data/ACCESS-OM2-01/01deg_jra55v140_iaf_cycle4/psi_tot.nc"
    "/scratch/y99/TMIP/data/ACCESS-OM2-01/01deg_jra55v140_iaf_cycle4_jra55v150_extension/psi_tot.nc"
    "/scratch/y99/TMIP/data/ACCESS-OM2-01/01deg_jra55v150_iaf_cycle1/psi_tot.nc"
]
model_str(file) = split(file, "/")[6]
resolution_str(file) = split(model_str(file), "-")[3]
jra55version_str(file) = occursin("jra55v150", file) ? "jra55v150" : "jra55v140" # CHECK that it's v1.4 when unspecified
cycle_match(file) = match(r"cycle(\d+)", file)
cycle_str(file) = (c = cycle_match(file); isnothing(c) ? 1 : parse(Int, c.captures[1]))
files_df = DataFrame()
id_str(model, cycle, jra) = "$(model)_$(jra)_cycle$(cycle)"
for file in files
    model = model_str(file)
    resolution = resolution_str(file)
    cycle = cycle_str(file)
    jra_version = jra55version_str(file)
    id = id_str(model, cycle, jra_version)

    push!(files_df, (; id, file, model, resolution, cycle, jra_version))
end
@show(files_df)

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
    # "lat ≤ 40°S" => Where(≤(-40)),
    "50°S ≤ lat ≤ 40°S" => -50 .. -40,
    "lat ≈ 40°S" => Near(-40),
    "lat ≈ 0°" => Near(0),
    "lat ≈ 20°N" => Near(20),
)
ρselectors = OrderedDict(
    # "any ρ" => Colon(),
    "ρ ≥ 1036 kg/m³" => Where(≥(1036.0)),
    # "ρ ≥ 1036.8 kg/m³" => Where(≥(1036.8)),
)


for row in eachrow(files_df)
    (; id, file, model, resolution, cycle, jra_version) = row
    @info "Loading $model $jra_version cycle $cycle overturning streamfunction"
    ds = open_dataset(file)
    psis[id] = set(readcubedata(ds.psi_tot), :time => Ti)
    meanpsis[id] = dropdims(mean(psis[id]; dims = Ti); dims = Ti)
    ρs[id] = ds.potrho
    lats[id] = ds.grid_yu_ocean
    times[id] = ds.time
    println("  Computing AABW transport metrics:")
    AABW = OrderedDict{String, Any}()
    idx = OrderedDict{String, Any}()
    psi = psis[id]
    (grid_yu_ocean, potrho) = otherdims(psi, (Ti,))
    for (latstr, latselector) in pairs(latselectors)
        println(" - $latstr:")
        AABW_sub = OrderedDict{String, Any}()
        idx_sub = OrderedDict{String, Any}()
        for (ρstr, ρselector) in pairs(ρselectors)
            # AABW_sub[ρstr] = minimum(psi[latselector(lats[id]), ρselector(ρs[id]), :], dims = (:grid_yu_ocean, :potrho)) |> vec
            psi_subdomain = psi[DimensionalData.rebuild(grid_yu_ocean, latselector), DimensionalData.rebuild(potrho, ρselector)]
            AABW_sub[ρstr] = minimum(psi_subdomain, dims = (grid_yu_ocean, potrho))
            idx_sub[ρstr] = argmin(psi_subdomain.data, dims = 1:(length(size(psi_subdomain)) - 1))
            println("   - $ρstr (size: $(size(AABW_sub[ρstr]))")
        end
        AABW[latstr] = AABW_sub
        idx[latstr] = idx_sub
    end
    AABWs[id] = AABW
    idxs[id] = idx
end

########
# Plot #
########

# colors = Makie.wong_colors()
resolutions = unique(files_df.resolution)
resolution° = Dict(
    "1" => "1°",
    "025" => "0.25°",
    "01" => "0.1°",
)
colors = Dict((res => c) for (res, c) in zip(resolutions, cgrad(:tab10, categorical = true)))

fig = Figure(;
    size = (400 * length(ρselectors), 300 * length(latselectors)),
    fonts = (; regular = "Dejavu"),
)
g = fig[1, 1] = GridLayout(length(latselectors), length(ρselectors))

for (icol, (ρstr, ρselector)) in enumerate(pairs(ρselectors))
    for (irow, (latstr, latselector)) in enumerate(pairs(latselectors))
        axis = (
            ylabel = "AABW transport (Sv)",
            limits = (nothing, nothing, 0, nothing),
            xticks = DateTime.(1900:10:2030),
            xtickformat = "yyyy",
        )
        # ax, plt = CairoMakie.lines(fig[irow, icol], [DateTime(1990)], [10.0]; axis=(; xticks = DateTime.(1900:20:2030), ylabel = "AABW transport (Sv)"))
        for row in eachrow(files_df)
            (; id, file, model, resolution, cycle, jra_version) = row
            AABW = AABWs[id][latstr][ρstr].data |> vec
            time = times[id].val
            # rolling averages
            AABW_decade = fill(NaN, size(AABW))
            AABW_rolling = mean.(rolling(AABW, 120; center = true))
            offset = AABW_rolling.offsets[1]
            AABW_decade[offset .+ (0:(length(AABW_rolling.parent) - 1))] = AABW_rolling.parent
            AABW_year = fill(NaN, size(AABW))
            AABW_rolling = mean.(rolling(AABW, 12; center = true))
            offset = AABW_rolling.offsets[1]
            AABW_year[offset .+ (0:(length(AABW_rolling.parent) - 1))] = AABW_rolling.parent
            color = colors[resolution]
            panel = g[irow, icol]
            # @show l = ax.layout
            # @show contents(ax)
            if isempty(contents(panel))
                CairoMakie.lines(panel, DateTime.(time), -AABW_year / 1.0e9; color, alpha = 0.3, axis) # Sv
            else
                lines!(panel, DateTime.(time), -AABW_year / 1.0e9; color, alpha = 0.3) # Sv
            end
            ax = content(panel)
            lines!(ax, DateTime.(time), -AABW_decade / 1.0e9; color, label = resolution°[resolution]) # Sv
            inan = .!isnan.(AABW_year)
            textoptions = (; text = string(cycle), color, fontsize = 8)
            text!(ax, first(DateTime.(time[inan])), first(-AABW_year[inan] / 1.0e9); textoptions..., align = (:right, :center), offset = (-5, 0))
            text!(ax, last(DateTime.(time[inan])), last(-AABW_year[inan] / 1.0e9); textoptions..., align = (:center, :center), offset = (+5, 0))
            # text!(ax, 0.5, 0.5; text = "foo", space = :relative)
        end
        ax = content(g[irow, icol])
        (icol == irow == 1) && axislegend(ax, position = :lb, merge = true)
        # ax.xticks = DateTimeTicks(6)

        myhidexdecorations!(ax, irow < length(latselectors))
        myhideydecorations!(ax, icol > 1)

        (icol == 1) && Label(g[irow, 0], latstr, tellheight = false, rotation = π / 2)

        # Add mean overturning streamfunction + AABW box for minimum
        for row in eachrow(files_df)
            (; id, file, model, resolution, cycle, jra_version) = row
            if resolution == "01" && cycle == 4
                break
            end
        end
        res = "01"
        idx = findfirst((files_df.resolution .== res) .& (files_df.cycle .== 4))
        @show id = files_df.id[idx]
        ρmin, ρmax = extrema(ρs[id].val)
        latmin, latmax = extrema(lats[id].val)
        ρmin -= eps(ρmin)
        # yscale = ReversibleScale(ρ -> exp(ρ - ρmin), x -> log(x) + ρmin, limits = (ρmin, ρmax))
        yscale = ReversibleScale(ρ -> (ρ - ρmin)^4, x -> x^(1 / 4) + ρmin, limits = (ρmin, ρmax))
        # yscale = ReversibleScale(ρ -> (ρ - ρmin)^2, x -> sqrt(x) + ρmin, limits = (ρmin, ρmax))
        (grid_yu_ocean, potrho) = dims(meanpsis[id])
        latbox = grid_yu_ocean[collect(extrema(DimensionalData.selectindices(grid_yu_ocean, latselector)))].val
        ρbox = potrho[collect(extrema(DimensionalData.selectindices(potrho, ρselector)))].val
        ax_inset = Axis(
            g[irow, icol];
            width = Relative(0.6),
            height = Relative(0.3),
            halign = 0.95,
            valign = 0.95,
            yreversed = true,
            limits = (extrema(lats[id].val)..., ρmin, ρmax),
            yscale = yscale,
            yticks = round.([ρ for ρ in ρbox if ρ ≉ ρmin], digits = 1),
            xticks = round.(Int, [l for l in latbox if l ≉ latmin]),
            alignmode = Outside(),
            xticklabelsize = 8,
            yticklabelsize = 8,
        )
        cf = contourf!(
            ax_inset, lats[id].val, ρs[id].val, meanpsis[id].data / 1.0e9; # Sv
            levels = -24:2:24,
            colormap = :delta,
            extendlow = :auto,
            extendhigh = :auto,
        )
        Colorbar(g[irow, icol], cf;
            width = Relative(0.6),
            height = Relative(0.3),
            halign = 0.95,
            valign = 0.95,)
        lines!(ax_inset, latbox[[1, 2, 2, 1, 1]], ρbox[[1, 1, 2, 2, 1]]; color = :red)
    end
    Label(g[0, icol], ρstr, tellwidth = false)
end

# linkxaxes!([panel.content for panel in g.content if panel.content isa Axis]...)
# rowgap!(g, 0.0)
# colgap!(g, 0.0)

# box1 = Box(g[2, 2, Makie.GridLayoutBase.Outer()], cornerradius = 0, color = (:tomato, 0.5), strokecolor = :transparent)
# box2 = Box(g[1, 1, Makie.GridLayoutBase.Outer()], cornerradius = 0, color = (:teal, 0.5), strokecolor = :transparent)

fig

imagefile = "/scratch/y99/TMIP/data/AABW_timeseries_allcycles.png"
save(imagefile, fig)
@info "Saved AABW timeseries plot to $imagefile"
imagefile = "/scratch/y99/TMIP/data/AABW_timeseries_allcycles.pdf"
save(imagefile, fig)
@info "Saved AABW timeseries plot to $imagefile"
