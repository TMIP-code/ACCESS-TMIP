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


# unpack model grid
(; lon, lat, zt, v3D,) = gridmetrics
lev = zt
maxlat = dropdims(maximum(lat, dims=1), dims=1)

# unpack indices
(; wet3D, N) = indices

# basins
basin_keys = (:ATL, :INDOPAC, :GBL, :SO)
basin_strs = ("Atlantic", "Indo-Pacific", "Global", "Southern Ocean")
OCEANS = OceanBasins.oceanpolygons()
isSO(lat, lon, o) = lat .< -30
isatlnoSO(lat, lon, o) = isatlantic(lat, lon, o) .& (lat .≥ -30)
isindopacific(lat, lon, o) = (isindian(lat, lon, o) .| ispacific(lat, lon, o)) .& (lat .≥ -30)
isglobal(lat, lon, o) = trues(size(lat))
basin_functions = (isatlnoSO, isindopacific, isglobal, isSO)
basin_values = (reshape(f(lat[:], lon[:], OCEANS), size(lat)) for f in basin_functions)
basins = (; (basin_keys .=> basin_values)...)
basin_latlims_values = [clamp.(0 .* (-5, +5) .+ extrema(lat[.!isnan.(v3D[:,:,1]) .& basin[:,:,1]]), -80, 80) for basin in basins]
basin_latlims = (; (basin_keys .=> basin_latlims_values)...)



# Compute MOC as max or min of stream function

ρ = 1035.0    # density (kg/m^3)

MOCdf = DataFrame(
    time_window = String[], experiment = String[], member = String[],
    AMOC = Float64[], AMOC_z = Float64[], AMOC_lat = Float64[],
    # IPMOC = Float64[], IPMOC_lat = Float64[], IPMOC_z = Float64[],
    USMOC = Float64[], USMOC_lat = Float64[], USMOC_z = Float64[],
    LSMOC = Float64[], LSMOC_lat = Float64[], LSMOC_z = Float64[],
)

for (time_window, experiment) in zip(time_windows, experiments)

    @show time_window

    for member in members

        @show member

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

        # AMOC
        basin = basins.ATL
        x2D = dropdims(reverse(nancumsum(reverse(nansum(basin .* ϕ.north, dims = 1), dims=3), dims = 3), dims=3), dims = 1) # kg/s
        x2Dmask = zonalaverage(1, gridmetrics; mask = basin) .> 0
        x2D[.!x2Dmask] .= 0
        x2D[:, 1:10] .= 0
        # convert to Sv
        x2D = x2D / 1e6 / ρ # Sv
        AMOCval, AMOCidx = findmin(x2D)

        # USMOC
        basin = basins.SO
        x2D = dropdims(reverse(nancumsum(reverse(nansum(basin .* ϕ.north, dims = 1), dims=3), dims = 3), dims=3), dims = 1) # kg/s
        x2Dmask = zonalaverage(1, gridmetrics; mask = basin) .> 0
        x2D[.!x2Dmask] .= 0
        x2D[:, 1:10] .= 0
        x2D[maxlat .≥ -30, :] .= 0

        # convert to Sv
        x2D = x2D / 1e6 / ρ # Sv
        USMOCval, USMOCidx = findmin(x2D)
        LSMOCval, LSMOCidx = findmax(x2D)

        push!(MOCdf, (time_window, experiment, member,
            AMOCval, zt[AMOCidx[2]], maxlat[AMOCidx[1]],
            USMOCval, zt[USMOCidx[2]], maxlat[USMOCidx[1]],
            LSMOCval, zt[LSMOCidx[2]], maxlat[LSMOCidx[1]]
        ))
    end

end


# Plot meridional overturning circulation for each basin

fig = Figure(size = (800, 500), fontsize = 18)

raincloudioptions = (
    # boxplot_nudge = -0.3,
    # boxplot_width = 0.3,
    # jitter_width = 0.3,
    clouds = violin,
    # hist_bins = 0:0.25:30,
    # cloud_width = 0.3,
    # gap = 0,
    markersize = 5,
    orientation = :horizontal,
)

axisoptions = (
    # ylabel="basin",
    xticks = 0:2:36,
    # limits = (12, 24, 0.5, 3.5),
    limits = (6, 34, 0.5, 3.5),
    yticks = (1:length(time_windows), ["$e\n$(t[4:7])s" for (e, t) in zip(reverse(experiments2), reverse(time_windows))]),
    yticklabelrotation = pi/2,
    yticksvisible = false,
    # yticklabelsvisible = false,
    ygridvisible = false,
    xtrimspine = true,
    rightspinevisible = false,
    topspinevisible = false,
    leftspinevisible = false,
)

textoptions = (
    align = (:center, :top),
    offset = (0, -25),
    fontsize = 12
)

textoptions2 = (
    align = (:center, :bottom),
    offset = (-50, +15),
    fontsize = 18,
    font = :bold,
)

# AMOC
ax = Axis(fig[1, 1]; axisoptions..., xlabel = "MOC (Sv)")
rp = rainclouds!(ax, reverse(MOCdf.time_window), reverse(-MOCdf.AMOC); raincloudioptions...)
sp = rp.plots[findfirst(x -> x isa Makie.Scatter, rp.plots)]
# for (i, (e, t)) in enumerate(zip(reverse(experiments2), reverse(time_windows)))
#     text!(ax, -mean(MOCdf.AMOC[MOCdf.time_window .== t]), i; text = "$e $(t[4:7])s", textoptions...)
# end
text!(ax, -mean(MOCdf.AMOC[MOCdf.time_window .== time_windows[1]]), 3; text = "AMOC", textoptions2..., color = sp.color)

# USMOC
USMOClimits = (6, 36, 0.5, 3.5)
# ax = Axis(fig[2, 1]; axisoptions..., xlabel = "upper-cell SMOC (Sv)", limits = USMOClimits)
rp = rainclouds!(ax, reverse(MOCdf.time_window), reverse(-MOCdf.USMOC); raincloudioptions...)
sp = rp.plots[findfirst(x -> x isa Makie.Scatter, rp.plots)]
# for (i, (e, t)) in enumerate(zip(reverse(experiments2), reverse(time_windows)))
#     text!(ax, -mean(MOCdf.USMOC[MOCdf.time_window .== t]), i; text = "$e $(t[4:7])s", textoptions...)
# end
text!(ax, -mean(MOCdf.USMOC[MOCdf.time_window .== time_windows[1]]), 3; text = "upper-cell\nSMOC", textoptions2..., color = sp.color, offset = (+50, 20))


# LSMOC
LSMOClimits = (6, 36, 0.5, 3.5)
rp = rainclouds!(ax, reverse(MOCdf.time_window), reverse(MOCdf.LSMOC); raincloudioptions...)
sp = rp.plots[findfirst(x -> x isa Makie.Scatter, rp.plots)]
# for (i, (e, t)) in enumerate(zip(reverse(experiments2), reverse(time_windows)))
#     text!(ax, mean(MOCdf.LSMOC[MOCdf.time_window .== t]), i; text = "$e $(t[4:7])s", textoptions...)
# end
text!(ax, mean(MOCdf.LSMOC[MOCdf.time_window .== time_windows[1]]), 3; text = "lower-cell\nSMOC", textoptions2..., color = sp.color, offset = (-50, 20))



fig

outputdir = "output/plots"

# save plot
outputfile = joinpath(outputdir, "$(model)_MOC_variability_for_Reviewer2.png")
@info "Saving MOC as image file:\n  $(outputfile)"
save(outputfile, fig)
outputfile = joinpath(outputdir, "$(model)_MOC_variability_for_Reviewer2.pdf")
@info "Saving MOC as image file:\n  $(outputfile)"
save(outputfile, fig)





