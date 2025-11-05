# qsub -I -P xv83 -l mem=47GB -l storage=scratch/gh0+scratch/xv83 -l walltime=02:00:00 -l ncpus=12

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
using GeometryBasics
using GeometryOps
using LibGEOS
using Format
# using LaTeXStrings

model = "ACCESS-ESM1-5"

time_window1 = "Jan2030-Dec2039"
time_window2 = "Jan2090-Dec2099"
experiment1 = parse(Int, time_window1[4:7]) ≤ 2010 ? "historical" : "ssp370"
experiment2 = parse(Int, time_window2[4:7]) ≤ 2010 ? "historical" : "ssp370"


members = ["r$(r)i1p1f1" for r in 1:40]
# members = ["r$(r)i1p1f1" for r in 1:3]

# Gadi directory for input files
# inputdirfun(member) = "/scratch/xv83/TMIP/data/$model/$experiment/all members/$(time_window)"
inputdir1 = "/scratch/xv83/TMIP/data/$model/$experiment1/all_members/$(time_window1)/cyclomonth"
inputdir2 = "/scratch/xv83/TMIP/data/$model/$experiment2/all_members/$(time_window2)/cyclomonth"
outputdir = inputdir2
mkpath(inputdir2)


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

# Load \Gamma out
κVdeep = 3.0e-5 # m^2/s
κVML = 1.0      # m^2/s
κH = 300.0      # m^2/s
κVdeep_str = "kVdeep" * format(κVdeep, conversion = "e")
κVML_str = "kVML" * format(κVML, conversion = "e")
κH_str = "kH" * format(κH, conversion = "d")
upwind = false
upwind_str = upwind ? "" : "_centered"
upwind_str2 = upwind ? "upwind" : "centered"
yearly = true
yearly_str = yearly ? "_yearly" : ""
yearly_str2 = yearly ? "(yearly)" : ""

# τ = timescales for which we plot ℰ(τ)
τs = [100, 300, 1000]

if yearly
    # Read 2030s files
    ℰ_files = ["/scratch/xv83/TMIP/data/$model/$experiment1/$member/$(time_window1)/seqeff$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str).nc" for member in members]
    ℰ_ds = open_mfdataset(DimArray(ℰ_files, Dim{:member}(members)))
    # Load the data for the selected timescales only (to avoid using too much memory)
    ℰ = readcubedata(ℰ_ds.seqeff[Ti = At(τs)])
    years = ℰ_ds.Ti |> Array
    ℰ_ensemblemean1 = dropdims(mean(ℰ, dims = :member), dims = :member)
    # Read 2090s files
    ℰ_files = ["/scratch/xv83/TMIP/data/$model/$experiment2/$member/$(time_window2)/seqeff$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str).nc" for member in members]
    ℰ_ds = open_mfdataset(DimArray(ℰ_files, Dim{:member}(members)))
    ℰ = readcubedata(ℰ_ds.seqeff[Ti = At(τs)])
    ℰ_ensemblemean2 = dropdims(mean(ℰ, dims = :member), dims = :member)
    ℰ_diff = ℰ_ensemblemean2 - ℰ_ensemblemean1
else
    # TODO update when periodic adjoint propagator runs are finished
    # Read 2030s files
    ℰ_files = ["/scratch/xv83/TMIP/data/$model/$experiment1/$member/$(time_window1)/calE$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str).nc" for member in members]
    ℰ_ds = open_mfdataset(DimArray(ℰ_files, Dim{:member}(members)))
    ℰ = readcubedata(ℰ_ds.calE)
    years = ℰ_ds.Ti |> Array
    ℰ_ensemblemean1 = dropdims(mean(ℰ, dims = 4), dims = 4)
    # Read 2090s files
    ℰ_files = ["/scratch/xv83/TMIP/data/$model/$experiment2/$member/$(time_window2)/calE$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str).nc" for member in members]
    ℰ_ds = open_mfdataset(DimArray(ℰ_files, Dim{:member}(members)))
    ℰ = readcubedata(ℰ_ds.calE)
    ℰ_ensemblemean2 = dropdims(mean(ℰ, dims = 4), dims = 4)
    ℰ_diff = ℰ_ensemblemean2 - ℰ_ensemblemean1
end


include("plotting_functions.jl")

usecontourf = false

axs = Array{Any, 2}(undef, (3, 2))
contours = Array{Any, 2}(undef, (3, 2))
nrows, ncols = size(axs)

fig = Figure(size = (ncols * 500, nrows * 250 + 100), fontsize = 18)

yticks = -60:30:60
xticks = -120:60:(120 + 360)

for (irow, year) in enumerate([100, 300, 1000])

    iyear = year + 1

    # Plot ensemble mean
    icol = 1
    levels = 0:10:100
    colormap = cgrad(:Zissou1Continuous, length(levels) - 1, categorical = true, rev = true)
    # colormap = cgrad(:Hiroshige, 10, categorical = true, rev = true)
    colorrange = extrema(levels)

    axs[irow, icol] = ax = Axis(fig[irow, icol]; yticks, xticks, xtickformat, ytickformat)

    contours[irow, icol] = if usecontourf
        plotcontourfmap!(ax, 100 * ℰ_ensemblemean2[:, :, irow], gridmetrics; levels, colormap)
    else
        plotmap!(ax, 100 * ℰ_ensemblemean2[:, :, irow], gridmetrics; colorrange, colormap) # <- need to fix wrapping longitude for contour levels
    end

    myhidexdecorations!(ax, irow < nrows)
    myhideydecorations!(ax, icol > 1)

    # Plot difference
    icol = 2
    levels = -50:10:50
    colormap = cgrad(cgrad(:tol_bu_rd, length(levels), categorical = true)[[1:(end ÷ 2 + 1); (end ÷ 2 + 1):end]], categorical = true)
    highclip = colormap[end]
    lowclip = colormap[1]
    colormap = cgrad(colormap[2:(end - 1)], categorical = true)
    colorrange = extrema(levels)

    axs[irow, icol] = ax = Axis(fig[irow, icol]; yticks, xticks, xtickformat, ytickformat)

    contours[irow, icol] = if usecontourf
        plotcontourfmap!(ax, 100 * ℰ_diff[:, :, irow], gridmetrics; levels, colormap)
    else
        contours[irow, icol] = plotmap!(ax, 100 * ℰ_diff[:, :, irow], gridmetrics; colorrange, colormap, highclip, lowclip) # <- need to fix wrapping longitude for contour levels
    end

    myhidexdecorations!(ax, irow < nrows)
    myhideydecorations!(ax, icol > 1)

    Label(fig[irow, 0]; text = "τ = $year years", rotation = π / 2, tellheight = false)

end

# ℰstr = rich("ℰ", superscript("—"ℰ̅ ‾))
ℰstr = rich("ℰ", rich("‾", offset = (-0.5, 0.15)))
# τstr = rich("τ", font = :italic)
𝒓 = rich("r", font = :bold_italic)
ℰfun = rich(ℰstr, "(", 𝒓, ", τ)")
label = rich("ensemble mean ", ℰfun, " (%)")
cb = Colorbar(fig[nrows + 1, 1], contours[1, 1]; label, vertical = false, flipaxis = false, ticks = 0:20:100)
cb.width = Relative(2 / 3)

label = rich("mean 2090s ", ℰfun, " − mean 2030s ", ℰfun, " (%)")
cb = Colorbar(fig[nrows + 1, 2], contours[1, 2]; label, vertical = false, flipaxis = false, ticks = -50:20:50, tickformat = divergingcbarticklabelformat)
cb.width = Relative(2 / 3)

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

Label(fig[0, 1:2]; text = "Climate Change Effect on Seafloor Sequestration Efficiency ($(length(members)) members)$(yearly_str2)", fontsize = 24, tellwidth = false)
rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)

# save plot
suffix = usecontourf ? "_ctrf" : ""

outputfile = joinpath(outputdir, "seqeff_diff$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str)_$(time_window2)$(suffix).png")
@info "Saving seqeff on sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)
outputfile = joinpath(outputdir, "seqeff_diff$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str)_$(time_window2)$(suffix).pdf")
@info "Saving seqeff on sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)
