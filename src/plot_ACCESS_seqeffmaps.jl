# qsub -I -P xv83 -l mem=64GB -l storage=scratch/gh0+scratch/xv83 -l walltime=02:00:00 -l ncpus=12

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

include("plotting_functions.jl")


model = "ACCESS-ESM1-5"

# for time_window in ["Jan1850-Dec1859", "Jan1990-Dec1999", "Jan2030-Dec2039", "Jan2090-Dec2099"]
time_window = "Jan2030-Dec2039"
    experiment = parse(Int, time_window[4:7]) ≤ 2010 ? "historical" : "ssp370"
    # experiment = "historical"
    # time_window = "Jan1850-Dec1859"
    # time_window = "Jan1990-Dec1999"
    # experiment = "ssp370"
    # time_window = "Jan2030-Dec2039"
    # time_window = "Jan2090-Dec2099"

    # members = ["r$(r)i1p1f1" for r in 1:40]
    members = ["r$(r)i1p1f1" for r in 1:3]

    # Gadi directory for input files
    # inputdirfun(member) = "/scratch/xv83/TMIP/data/$model/$experiment/all members/$(time_window)"
    inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/all_members/$(time_window)/cyclomonth"
    outputdir = inputdir
    mkpath(inputdir)



    mlotst_files = ["/scratch/xv83/TMIP/data/$model/$experiment/$member/$time_window/cyclomonth/mlotst.nc" for member in members]
    mlotst_ds = open_mfdataset(DimArray(mlotst_files, Dim{:member}(members)))
    mlotst = readcubedata(mlotst_ds.mlotst)

    mlotst_yearmax_ensemblemean = dropdims(maximum(maximum(mlotst, dims=:month), dims=:member), dims=(:month, :member)) .-
        dropdims(minimum(maximum(mlotst, dims=:month), dims=:member), dims=(:month, :member))

    # Load areacello and volcello for grid geometry
    fixedvarsinputdir = "/scratch/xv83/TMIP/data/$model"
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
    (; lon_vertices, lat_vertices, lon, lat, zt, v3D, thkcello, Z3D) = gridmetrics
    lev = zt
    # Make indices
    indices = makeindices(gridmetrics.v3D)
    (; wet3D, N) = indices



    ℰ_files = ["/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)/calE.nc" for member in members]
    ℰ_ds = open_mfdataset(DimArray(ℰ_files, Dim{:member}(members)))
    ℰ = readcubedata(ℰ_ds.calE)



    years = ℰ_ds.Ti |> Array
    iyears =

    ℰ_ensemblemean = dropdims(mean(ℰ, dims = 4), dims = 4)
    ℰ_ensemblemax = dropdims(maximum(ℰ, dims = 4), dims = 4)
    ℰ_ensemblemin = dropdims(minimum(ℰ, dims = 4), dims = 4)
    ℰ_ensemblerange = ℰ_ensemblemax - ℰ_ensemblemin



    axs = Array{Any,2}(undef, (3, 2))
    contours = Array{Any,2}(undef, (3, 2))
    nrows, ncols = size(axs)

    fig = Figure(size = (ncols * 500, nrows * 250 + 100), fontsize = 18)

    yticks = -60:30:60
    xticks = -120:60:120

    for (irow, year) in enumerate([100, 300, 1000])

        iyear = year + 1

        # Plot ensemble mean
        icol = 1
        levels = 0:10:100
        colormap = cgrad(:Zissou1Continuous, 10, categorical = true, rev = true)
        # colormap = cgrad(:Hiroshige, 10, categorical = true, rev = true)
        colorrange = extrema(levels)

        axs[irow, icol] = ax = Axis(fig[irow, icol]; yticks, xticks, xtickformat, ytickformat)

        contours[irow, icol] = plotmap!(ax, 100 * ℰ_ensemblemean[:, :, iyear], gridmetrics; colorrange, colormap) # <- need to fix wrapping longitude for contour levels

        myhidexdecorations!(ax, irow < nrows)
        myhideydecorations!(ax, icol > 1)

        # Plot ensemble range
        icol = 2
        levels = 0:5:25
        colormap = cgrad(:linear_wyor_100_45_c55_n256, 5, categorical = true)
        colorrange = extrema(levels)

        axs[irow, icol] = ax = Axis(fig[irow, icol]; yticks, xticks, xtickformat, ytickformat)

        contours[irow, icol] = plotmap!(ax, 100 * ℰ_ensemblerange[:, :, iyear], gridmetrics; colorrange, colormap) # <- need to fix wrapping longitude for contour levels

        myhidexdecorations!(ax, irow < nrows)
        myhideydecorations!(ax, icol > 1)

        Label(fig[irow, 0]; text = "after $year years", rotation = π/2, tellheight = false)

    end


    label = rich("ensemble mean ℰ (%)")
    cb = Colorbar(fig[nrows + 1, 1], contours[1, 1]; label, vertical = false, flipaxis = false, ticks = 0:20:100)
    cb.width = Relative(2/3)

    label = rich("ensemble range, max ℰ − min ℰ (%)")
    cb = Colorbar(fig[nrows + 1, 2], contours[1, 2]; label, vertical = false, flipaxis = false, ticks = 0:5:25)
    cb.width = Relative(2/3)

    # column labels
    # Label(fig[0, 1]; text = "ensemble mean", tellwidth = false)
    # Label(fig[0, 2]; text = "ensemble range (internal variability)", tellwidth = false)

    labels = [
        "a" "d"
        "b" "e"
        "c" "f"
    ]

    for (ax, label) in zip(axs, labels)
        txt = text!(ax, 0, 1; text = "$label", labeloptions..., strokecolor = :white, strokewidth = 3)
        translate!(txt, 0, 0, 100)
        txt = text!(ax, 0, 1; text = "$label", labeloptions...)
        translate!(txt, 0, 0, 100)
    end

    Label(fig[0, 1:2]; text = "$(time_window[4:7])s Seafloor Sequestration Efficiency ($(length(members)) members)", fontsize = 24, tellwidth = false)
    # rowgap!(fig.layout, 10)
    # colgap!(fig.layout, 10)

    # save plot
    outputfile = joinpath(outputdir, "seqeff_$(time_window).png")
    @info "Saving seqeff on sea floor as image file:\n  $(outputfile)"
    save(outputfile, fig)
    outputfile = joinpath(outputdir, "seqeff_$(time_window).pdf")
    @info "Saving seqeff on sea floor as image file:\n  $(outputfile)"
    save(outputfile, fig)

