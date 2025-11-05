# # qsub -I -P xv83 -q express -l mem=47GB -l storage=scratch/gh0+scratch/xv83 -l walltime=01:00:00 -l ncpus=12
# # This is Fig. 1 in Pasquier et al. (GRL, 2025)

# using Pkg
# Pkg.activate(".")
# Pkg.instantiate()

# using OceanTransportMatrixBuilder
# using NetCDF
# using YAXArrays
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
# using Format

# model = "ACCESS-ESM1-5"

# time_window_1850s = "Jan1850-Dec1859"
# time_window_2030s = "Jan2030-Dec2039"
# time_window_2090s = "Jan2090-Dec2099"
# time_windows = ["1850s", "2030s", "2090s"]
# experiments = ["historical", "ssp370", "ssp370"]
# experiments2 = ["historical", "SSP3-7.0", "SSP3-7.0"]

# members = ["r$(r)i1p1f1" for r in 1:40]
# # members = ["r$(r)i1p1f1" for r in 1:3]

# # Load areacello and volcello for grid geometry
# fixedvarsinputdir = "/scratch/xv83/TMIP/data/$model"
# volcello_ds = open_dataset(joinpath(fixedvarsinputdir, "volcello.nc"))
# areacello_ds = open_dataset(joinpath(fixedvarsinputdir, "areacello.nc"))

# # Load fixed variables in memory
# areacello = readcubedata(areacello_ds.areacello)
# volcello = readcubedata(volcello_ds.volcello)
# lon = readcubedata(volcello_ds.lon)
# lat = readcubedata(volcello_ds.lat)
# lev = volcello_ds.lev
# # Identify the vertices keys (vary across CMIPs / models)
# volcello_keys = propertynames(volcello_ds)
# lon_vertices_key = volcello_keys[findfirst(x -> occursin("lon", x) & occursin("vert", x), string.(volcello_keys))]
# lat_vertices_key = volcello_keys[findfirst(x -> occursin("lat", x) & occursin("vert", x), string.(volcello_keys))]
# lon_vertices = readcubedata(getproperty(volcello_ds, lon_vertices_key))
# lat_vertices = readcubedata(getproperty(volcello_ds, lat_vertices_key))
# # Make makegridmetrics
# gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
# (; lon_vertices, lat_vertices, lon, lat, zt, v3D, thkcello, Z3D) = gridmetrics
# lev = zt
# # Make indices
# indices = makeindices(gridmetrics.v3D)
# (; wet3D, N) = indices

# experiment_dir = "/scratch/xv83/TMIP/data/$model/$(experiments[1])"
# mlotst_files_1850s = [joinpath(experiment_dir, member, time_window_1850s, "cyclomonth", "mlotst.nc") for member in members]
# mlotst_1850s_ds = open_mfdataset(DimArray(mlotst_files_1850s, Dim{:member}(members)))
# mlotst_1850s = readcubedata(mlotst_1850s_ds.mlotst)

# experiment_dir = "/scratch/xv83/TMIP/data/$model/$(experiments[2])"
# mlotst_files_2030s = [joinpath(experiment_dir, member, time_window_2030s, "cyclomonth", "mlotst.nc") for member in members]
# mlotst_2030s_ds = open_mfdataset(DimArray(mlotst_files_2030s, Dim{:member}(members)))
# mlotst_2030s = readcubedata(mlotst_2030s_ds.mlotst)

# experiment_dir = "/scratch/xv83/TMIP/data/$model/$(experiments[3])"
# mlotst_files_2090s = [joinpath(experiment_dir, member, time_window_2090s, "cyclomonth", "mlotst.nc") for member in members]
# mlotst_2090s_ds = open_mfdataset(DimArray(mlotst_files_2090s, Dim{:member}(members)))
# mlotst_2090s = readcubedata(mlotst_2090s_ds.mlotst)

# mlotst_1850s_yearlymax = dropdims(maximum(mlotst_1850s, dims = :month), dims = :month)
# mlotst_2030s_yearlymax = dropdims(maximum(mlotst_2030s, dims = :month), dims = :month)
# mlotst_2090s_yearlymax = dropdims(maximum(mlotst_2090s, dims = :month), dims = :month)

# mlotst_1850s_yearlymax_ensemblemean = dropdims(mean(mlotst_1850s_yearlymax, dims = :member), dims = :member)
# mlotst_1850s_yearlymax_ensemblemax = dropdims(maximum(mlotst_1850s_yearlymax, dims = :member), dims = :member)
# mlotst_1850s_yearlymax_ensemblemin = dropdims(minimum(mlotst_1850s_yearlymax, dims = :member), dims = :member)
# mlotst_1850s_yearlymax_ensemblerange = mlotst_1850s_yearlymax_ensemblemax - mlotst_1850s_yearlymax_ensemblemin

# mlotst_2030s_yearlymax_ensemblemean = dropdims(mean(mlotst_2030s_yearlymax, dims = :member), dims = :member)
# mlotst_2030s_yearlymax_ensemblemax = dropdims(maximum(mlotst_2030s_yearlymax, dims = :member), dims = :member)
# mlotst_2030s_yearlymax_ensemblemin = dropdims(minimum(mlotst_2030s_yearlymax, dims = :member), dims = :member)
# mlotst_2030s_yearlymax_ensemblerange = mlotst_2030s_yearlymax_ensemblemax - mlotst_2030s_yearlymax_ensemblemin

# mlotst_2090s_yearlymax_ensemblemean = dropdims(mean(mlotst_2090s_yearlymax, dims = :member), dims = :member)
# mlotst_2090s_yearlymax_ensemblemax = dropdims(maximum(mlotst_2090s_yearlymax, dims = :member), dims = :member)
# mlotst_2090s_yearlymax_ensemblemin = dropdims(minimum(mlotst_2090s_yearlymax, dims = :member), dims = :member)
# mlotst_2090s_yearlymax_ensemblerange = mlotst_2090s_yearlymax_ensemblemax - mlotst_2090s_yearlymax_ensemblemin

# MLD_ensemble_means = [mlotst_1850s_yearlymax_ensemblemean, mlotst_2030s_yearlymax_ensemblemean, mlotst_2090s_yearlymax_ensemblemean]
# MLD_ensemble_ranges = [mlotst_1850s_yearlymax_ensemblerange, mlotst_2030s_yearlymax_ensemblerange, mlotst_2090s_yearlymax_ensemblerange]


include("plotting_functions.jl")

usecontourf = false

axs = Array{Any, 2}(undef, (3, 2))
contours = Array{Any, 2}(undef, (3, 2))
nrows, ncols = size(axs)

fig = Figure(size = (ncols * 500, nrows * 250 + 100), fontsize = 18)

yticks = -60:30:60
xticks = -120:60:(120 + 360)

# myscale = ReversibleScale(
#     x -> sign(x) * log10(abs(x / 5) + 1),
#     x -> sign(x) * (exp10(abs(x)) - 1) * 5;
#     # x -> x,
#     # x -> x;
#     limits=(0f0, 3f0),
#     name=:myscale
# )

# levels = [0, 50, 70, 100, 140, 200, 300, 500, 700, 1000, 1400, 2000, 3000]
levels = [0, 50, 100, 200, 500, 1000, 2000]
colorscale = mk_piecewise_linear(levels)

colorrange = extrema(levels)
# pseudocolorrange = myscale.(colorrange)
colormap = cgrad(:thermal, length(levels); categorical = true)
extendlow = lowclip = nothing
extendhigh = highclip = colormap[end]
colormap = cgrad(colormap[1:(end - 1)], categorical = true)

# colormap = cgrad(:tol_ylorbr, length(levels); categorical=true)
# lowclip = nothing
# highclip = colormap[end]
# colormap = cgrad(colormap[1:end-1], categorical=true)
# # pseudologlevels = myscale.(levels)

for (irow, (MLD_ensemble_mean, MLD_ensemble_range)) in enumerate(zip(MLD_ensemble_means, MLD_ensemble_ranges))
    # Plot ensemble mean
    icol = 1
    axs[irow, icol] = ax = Axis(fig[irow, icol]; yticks, xticks, xtickformat, ytickformat, aspect = DataAspect())
    contours[irow, icol] = if usecontourf
        plotcontourfmap!(ax, MLD_ensemble_mean, gridmetrics; levels, colormap, extendhigh, colorscale)
    else
        plotmap!(ax, MLD_ensemble_mean, gridmetrics; colorrange, colormap, highclip, colorscale) # <- need to fix wrapping longitude for contour levels
    end
    myhidexdecorations!(ax, irow < nrows)
    myhideydecorations!(ax, icol > 1)

    # Plot ensemble range
    icol = 2
    axs[irow, icol] = ax = Axis(fig[irow, icol]; yticks, xticks, xtickformat, ytickformat, aspect = DataAspect())
    contours[irow, icol] = if usecontourf
        plotcontourfmap!(ax, MLD_ensemble_range, gridmetrics; levels, colormap, colorscale, extendhigh)
    else
        plotmap!(ax, MLD_ensemble_range, gridmetrics; colorrange, colormap, highclip, colorscale) # <- need to fix wrapping longitude for contour levels
    end
    myhidexdecorations!(ax, irow < nrows)
    myhideydecorations!(ax, icol > 1)

    Label(fig[irow, 0]; text = "$(experiments2[irow]) $(time_windows[irow])", rotation = π / 2, tellheight = false)
end


label = "ensemble mean MLD (m)"
# cb = Colorbar(fig[nrows + 1, 1], contours[1, 1]; label, vertical = false, flipaxis = false, ticks = levels)
# cb.width = Relative(2/3)

cb = Colorbar(
    fig[nrows + 1, 1];
    limits = (1, length(levels)),
    label,
    vertical = false,
    flipaxis = false,
    colormap,
    highclip,
    ticks = (1:length(levels), string.(levels)),
)
cb.width = Relative(2 / 3)

label = "ensemble range MLD (m)"
cb = Colorbar(
    fig[nrows + 1, 2];
    limits = (1, length(levels)),
    label,
    vertical = false,
    flipaxis = false,
    colormap,
    highclip,
    ticks = (1:length(levels), string.(levels)),
)
cb.width = Relative(2 / 3)
# cb = Colorbar(fig[nrows + 1, 2], contours[1, 2]; label, vertical = false, flipaxis = false, ticks = levels, scale = colorscale)
# cb.width = Relative(2/3)

# column labels
# Label(fig[0, 1]; text = "ensemble mean", tellwidth = false)
# Label(fig[0, 2]; text = "ensemble range (internal variability)", tellwidth = false)

labels = [
    "a" "b"
    "c" "d"
    "e" "f"
]

labeloptions = (
    font = :bold,
    align = (:left, :bottom),
    offset = (5, 2),
    space = :relative,
    fontsize = 24,
)

for (ax, label) in zip(axs, labels)
    txt = text!(ax, 0, 0; text = "$label", labeloptions..., strokecolor = :white, strokewidth = 3)
    translate!(txt, 0, 0, 100)
    txt = text!(ax, 0, 0; text = "$label", labeloptions...)
    translate!(txt, 0, 0, 100)
end

# Label(fig[0, 1:2]; text = "$(time_window[4:7])s Seafloor Sequestration Efficiency ($(length(members)) members)$(yearly_str2)", fontsize = 24, tellwidth = false)
rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)

colsize!(fig.layout, 1, Aspect(1, 2.0))
colsize!(fig.layout, 2, Aspect(1, 2.0))

# save plot
suffix = usecontourf ? "_ctrf" : ""

resize_to_layout!(fig)

outputdir = "/scratch/xv83/TMIP/data/$model/$(experiments[2])/all_members"


outputfile = joinpath(outputdir, "MLDs.png")
@info "Saving MLD image file:\n  $(outputfile)"
save(outputfile, fig)
outputfile = joinpath(outputdir, "MLDs.pdf")
@info "Saving MLD image file:\n  $(outputfile)"
save(outputfile, fig)
