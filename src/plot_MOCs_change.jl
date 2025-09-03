# qsub -I -P xv83 -q express -l mem=47GB -l storage=scratch/gh0+scratch/xv83 -l walltime=01:00:00 -l ncpus=12
# This is Fig. S7? in Pasquier et al. (GRL, 2025) # TODO: check figure number

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
using Unitful: s, yr, m, kg
using CairoMakie
using GeoMakie
using Interpolations
using OceanBasins
using Statistics
using StatsBase
using NaNStatistics
using Format
using GeometryBasics
using LibGEOS
using GeometryOps

include("plotting_functions.jl")

model = "ACCESS-ESM1-5"

time_windows = [
    "Jan1850-Dec1859",
    "Jan2030-Dec2039",
    "Jan2090-Dec2099",
]

experiments = [
    "historical",
    "ssp370",
    "ssp370",
]
experiments2 = [
    "historical",
    "SSP3-7.0",
    "SSP3-7.0",
]

members = ["r$(r)i1p1f1" for r in 1:40]

lumpby = "month"
months = 1:12


# Gadi directory for input files
fixedvarsinputdir = "/scratch/xv83/TMIP/data/$model"
# Load areacello and volcello for grid geometry
volcello_ds = open_dataset(joinpath(fixedvarsinputdir, "volcello.nc"))
areacello_ds = open_dataset(joinpath(fixedvarsinputdir, "areacello.nc"))


# Load fixed variables in memory
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

# Make makegridmetrics
gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
(; lon_vertices, lat_vertices) = gridmetrics

# Make indices
indices = makeindices(gridmetrics.v3D)


ϕsnorth = map(time_windows, experiments) do time_window, experiment

    ϕsnorth_ensemblemean = mean(map(members) do member

        inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)"
        cycloinputdir = joinpath(inputdir, "cyclomonth")
        umo_ds = open_dataset(joinpath(cycloinputdir, "umo.nc"))
        vmo_ds = open_dataset(joinpath(cycloinputdir, "vmo.nc"))
        ψᵢGM_ds = open_dataset(joinpath(cycloinputdir, "tx_trans_gm.nc"))
        ψⱼGM_ds = open_dataset(joinpath(cycloinputdir, "ty_trans_gm.nc"))
        ψᵢsubmeso_ds = open_dataset(joinpath(cycloinputdir, "tx_trans_submeso.nc"))
        ψⱼsubmeso_ds = open_dataset(joinpath(cycloinputdir, "ty_trans_submeso.nc"))
        mlotst_ds = open_dataset(joinpath(cycloinputdir, "mlotst.nc"))

        mean_days_in_month = umo_ds.mean_days_in_month |> Array
        w = Weights(mean_days_in_month)

        umo = dropdims(mean(readcubedata(umo_ds.umo), w; dims=:month), dims=:month)
        vmo = dropdims(mean(readcubedata(vmo_ds.vmo), w; dims=:month), dims=:month)

        ψᵢGM = dropdims(mean(readcubedata(ψᵢGM_ds.tx_trans_gm), w; dims=:month), dims=:month)
        ψⱼGM = dropdims(mean(readcubedata(ψⱼGM_ds.ty_trans_gm), w; dims=:month), dims=:month)
        ψᵢsubmeso = dropdims(mean(readcubedata(ψᵢsubmeso_ds.tx_trans_submeso), w; dims=:month), dims=:month)
        ψⱼsubmeso = dropdims(mean(readcubedata(ψⱼsubmeso_ds.ty_trans_submeso), w; dims=:month), dims=:month)

        # Replace missing values and convert to arrays
        # I think latest YAXArrays converts _FillValues to missing
        ψᵢGM = replace(ψᵢGM |> Array, missing => 0) .|> Float64
        ψⱼGM = replace(ψⱼGM |> Array, missing => 0) .|> Float64
        ψᵢsubmeso = replace(ψᵢsubmeso |> Array, missing => 0) .|> Float64
        ψⱼsubmeso = replace(ψⱼsubmeso |> Array, missing => 0) .|> Float64

        # Also remove missings in umo and vmo
        umo = replace(umo, missing => 0)
        vmo = replace(vmo, missing => 0)

        # Take the vertical diff of zonal/meridional transport diagnostics to get their mass transport
        (nx, ny, _) = size(ψᵢGM)
        ϕᵢGM = diff([fill(0.0, nx, ny, 1);;; ψᵢGM |> Array], dims=3)
        ϕⱼGM = diff([fill(0.0, nx, ny, 1);;; ψⱼGM |> Array], dims=3)
        ϕᵢsubmeso = diff([fill(0.0, nx, ny, 1);;; ψᵢsubmeso |> Array], dims=3)
        ϕⱼsubmeso = diff([fill(0.0, nx, ny, 1);;; ψⱼsubmeso |> Array], dims=3)

        # TODO fix incompatible dimensions betwewen umo and ϕᵢGM/ϕᵢsubmeso Dim{:i} and Dim{:xu_ocean}
        ϕ = let umo = umo + ϕᵢGM + ϕᵢsubmeso, vmo = vmo + ϕⱼGM + ϕⱼsubmeso
            facefluxesfrommasstransport(; umo, vmo, gridmetrics, indices)
        end

        return ϕ.north
    end)

end

# unpack model grid
(; lon, lat, zt, v3D,) = gridmetrics
lev = zt
# unpack indices
(; wet3D, N) = indices

# basins
basin_keys = (:ATL, :INDOPAC, :GBL)
basin_strs = ("Atlantic", "Indo-Pacific", "Global")
OCEANS = OceanBasins.oceanpolygons()
isatlnoSO(lat, lon, o) = isatlantic(lat, lon, o) .& (lat .≥ -30)
isindopacific(lat, lon, o) = (isindian(lat, lon, o) .| ispacific(lat, lon, o)) .& (lat .≥ -30)
isglobal(lat, lon, o) = trues(size(lat))
basin_functions = (isatlnoSO, isindopacific, isglobal)
basin_values = (reshape(f(lat[:], lon[:], OCEANS), size(lat)) for f in basin_functions)
basins = (; (basin_keys .=> basin_values)...)
basin_latlims_values = [clamp.(0 .* (-5, +5) .+ extrema(lat[.!isnan.(v3D[:,:,1]) .& basin[:,:,1]]), -80, 80) for basin in basins]
basin_latlims = (; (basin_keys .=> basin_latlims_values)...)



# Plot meridional overturning circulation for each basin

ρ = 1035.0    # density (kg/m^3)

levels = -24:2:24
colormap = cgrad(:curl, length(levels) + 1; categorical=true, rev=true)
extendlow = colormap[1]
extendhigh = colormap[end]
colormap = cgrad(colormap[2:end-1]; categorical=true)

fig = Figure(size = (1000, 250length(ϕsnorth)), fontsize = 18)
axs = Array{Any, 2}(undef, (length(ϕsnorth), length(basin_keys)))
contours = Array{Any, 2}(undef, (length(ϕsnorth), length(basin_keys)))
for (icol, (basin_key, basin)) in enumerate(pairs(basins))

    for (irow, (x3D, time_window, experiment)) in enumerate(zip(ϕsnorth, time_windows, experiments2))

        x2D = dropdims(reverse(nancumsum(reverse(nansum(basin .* x3D, dims = 1), dims=3), dims = 3), dims=3), dims = 1) # kg/s
        x2Dmask = zonalaverage(1, gridmetrics; mask = basin) .> 0
        x2D[.!x2Dmask] .= NaN

        # convert to Sv
        x2D = x2D / 1e6 / ρ # Sv

        local ax = Axis(fig[irow, icol],
            backgroundcolor = :lightgray,
            xgridvisible = true, ygridvisible = true,
            xgridcolor = (:black, 0.05), ygridcolor = (:black, 0.05),
            ylabel = "depth (m)")

        X = dropdims(maximum(lat, dims=1), dims=1)
        Y = zt
        Z = x2D
        co = contourf!(ax, X, Y, Z;
            levels,
            colormap,
            nan_color = :lightgray,
            extendlow,
            extendhigh,
        )
        translate!(co, 0, 0, -100)
        contours[irow, icol] = co

        xlim = basin_latlims[basin_key]
        # basin2 = LONGTEXT[basin]

        ax.yticks = (ztick, zticklabel)
        xticks = -90:30:90
        ax.xticks = (xticks, latticklabel.(xticks))
        ylims!(ax, zlim)
        # xlims!(ax, (-90, 90))
        xlims!(ax, xlim)


        hidexdecorations!(ax,
            label = irow < size(axs, 1), ticklabels = irow < size(axs, 1),
            ticks = irow < size(axs, 1), grid = false)
        hideydecorations!(ax,
            label = icol > 1, ticklabels = icol > 1,
            ticks = icol > 1, grid = false)

        (icol == 1) && Label(fig[irow, 0]; text = "$experiment $(time_window[4:7])s", rotation = π/2, tellheight = false)

        axs[irow, icol] = ax
    end

    Label(fig[0, icol], "$(basin_strs[icol]) MOC", tellwidth = false)


end

cb = Colorbar(fig[1:size(axs, 1), length(basin_keys) + 1], contours[1, 1];
    vertical = true, flipaxis = true,
    # ticks = (, cbarticklabelformat.(levels)),
    tickformat = x -> map(t -> replace(format("{:+d}", t), "-" => "−"), x),
    label = "MOC (Sv)",
    )
cb.height = Relative(0.666)

for (ibasin, xlims) in enumerate(basin_latlims)
    # Label(f[0, ibasin], LONGTEXT[basin], fontsize=20, tellwidth=false)
    colsize!(fig.layout, ibasin, Auto(xlims[2] - xlims[1]))
end

rowgap!(fig.layout, 15)
# rowgap!(fig.layout, 4, 15)
# # # rowgap!(fig.layout, 5, 10)
colgap!(fig.layout, 15)
# title = "$model ensemble-mean MOC"
# Label(fig[-1, 1:3]; text = title, fontsize = 20, tellwidth = false)

labels = [
    "a" "b" "c"
    "d" "e" "f"
    "g" "h" "i"
]

labeloptions = (
    font = :bold,
    align = (:left, :bottom),
    offset = (5, 2),
    space = :relative,
    fontsize = 24
)

for (ax, label) in zip(axs, labels)
    txt = text!(ax, 0, 0; text = label, labeloptions..., strokecolor = :white, strokewidth = 3)
    translate!(txt, 0, 0, 100)
    txt= text!(ax, 0, 0; text = label, labeloptions...)
    translate!(txt, 0, 0, 100)
end

fig

outputdir = "output/plots"

# save plot
outputfile = joinpath(outputdir, "$(model)_MOC_change_for_Reviewer2.png")
@info "Saving MOC as image file:\n  $(outputfile)"
save(outputfile, fig)
outputfile = joinpath(outputdir, "$(model)_MOC_change_for_Reviewer2.pdf")
@info "Saving MOC as image file:\n  $(outputfile)"
save(outputfile, fig)





