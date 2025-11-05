# qsub -I -P xv83 -q express -l mem=47GB -l storage=scratch/gh0+scratch/xv83 -l walltime=01:00:00 -l ncpus=12
# This is Fig. 1 in Pasquier et al. (GRL, 2025)

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

model = "ACCESS-ESM1-5"

time_window = "Jan2030-Dec2039"
experiment = parse(Int, time_window[4:7]) ≤ 2010 ? "historical" : "ssp370"

members = ["r$(r)i1p1f1" for r in 1:40]
# members = ["r$(r)i1p1f1" for r in 1:3]

# Gadi directory for input files
# inputdirfun(member) = "/scratch/xv83/TMIP/data/$model/$experiment/all members/$(time_window)"
inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/all_members/$(time_window)/cyclomonth"
outputdir = inputdir
mkpath(inputdir)


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

# Preferred diffusivities
κVdeep = 3.0e-5 # m^2/s
κVML = 1.0      # m^2/s
κH = 300.0      # m^2/s
κVdeep_str = "kVdeep" * format(κVdeep, conversion = "e")
κVML_str = "kVML" * format(κVML, conversion = "e")
κH_str = "kH" * format(κH, conversion = "d")
upwind = false
upwind_str = upwind ? "" : "_centered"
upwind_str2 = upwind ? "upwind" : "centered"

# Use yearly time-stepped simulations or monthly ones?
yearly = false
yearly_str = yearly ? "_yearly" : ""
yearly_str2 = yearly ? "(yearly)" : ""

# τ = timescales for which we plot ℰ(τ)
τs = [100, 300, 1000]

if yearly
    ℰ_files = ["/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)/seqeff$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str).nc" for member in members]
    ℰ_ds = open_mfdataset(DimArray(ℰ_files, Dim{:member}(members)))
    # Load the data for the selected timescales only (to avoid using too much memory)
    ℰ = readcubedata(ℰ_ds.seqeff[Ti = At(τs)])
else
    ℰ_files = ["/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)/calE$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str).nc" for member in members]
    ℰ_ds = open_mfdataset(DimArray(ℰ_files, Dim{:member}(members)))
    ℰ = readcubedata(ℰ_ds.calE[Ti = At(τs)])
end


years = ℰ_ds.Ti |> Array

ℰ_ensemblemean = dropdims(mean(ℰ, dims = 4), dims = 4)
ℰ_ensemblemax = dropdims(maximum(ℰ, dims = 4), dims = 4)
ℰ_ensemblemin = dropdims(minimum(ℰ, dims = 4), dims = 4)
ℰ_ensemblerange = ℰ_ensemblemax - ℰ_ensemblemin


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

    axs[irow, icol] = ax = Axis(fig[irow, icol]; yticks, xticks, xtickformat, ytickformat, aspect = DataAspect())

    contours[irow, icol] = if usecontourf
        plotcontourfmap!(ax, 100 * ℰ_ensemblemean[:, :, irow], gridmetrics; levels, colormap)
    else
        plotmap!(ax, 100 * ℰ_ensemblemean[:, :, irow], gridmetrics; colorrange, colormap) # <- need to fix wrapping longitude for contour levels
    end

    myhidexdecorations!(ax, irow < nrows)
    myhideydecorations!(ax, icol > 1)

    # Plot ensemble range
    icol = 2
    levels = 0:10:60
    colormap = cgrad(:tol_ylorbr, length(levels), categorical = true)
    highclip = colormap[end]
    colormap = cgrad(colormap[1:(end - 1)], categorical = true)
    colorrange = extrema(levels)

    axs[irow, icol] = ax = Axis(fig[irow, icol]; yticks, xticks, xtickformat, ytickformat, aspect = DataAspect())

    contours[irow, icol] = if usecontourf
        plotcontourfmap!(ax, 100 * ℰ_ensemblerange[:, :, irow], gridmetrics; levels, colormap)
    else
        contours[irow, icol] = plotmap!(ax, 100 * ℰ_ensemblerange[:, :, irow], gridmetrics; colorrange, colormap, highclip) # <- need to fix wrapping longitude for contour levels
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
label = rich("ensemble mean $(time_window[4:7])s sequestration efficiency, ", ℰfun, " (%)")
cb = Colorbar(fig[nrows + 1, 1], contours[1, 1]; label, vertical = false, flipaxis = false, ticks = 0:20:100)
cb.width = Relative(2 / 3)

label = rich("ensemble range, max ", ℰfun, " − min ", ℰfun, " (%)")
cb = Colorbar(fig[nrows + 1, 2], contours[1, 2]; label, vertical = false, flipaxis = false, ticks = 0:10:60)
cb.width = Relative(2 / 3)

# column labels
# Label(fig[0, 1]; text = "ensemble mean", tellwidth = false)
# Label(fig[0, 2]; text = "ensemble range (internal variability)", tellwidth = false)

labels = [
    "a" "d"
    "b" "e"
    "c" "f"
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

outputfile = joinpath(outputdir, "seqeff$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str)_$(time_window)$(suffix).png")
@info "Saving seqeff on sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)
outputfile = joinpath(outputdir, "seqeff$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str)_$(time_window)$(suffix).pdf")
@info "Saving seqeff on sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)

# Save the data to be uploaded with paper
metadata = Dict(
    "description" => "Sequestration efficiency as plotted in Fig. 1 in Pasquier et al. (2025)",
    "model" => model,
    "experiment" => experiment,
    "time window" => time_window,
    "unit" => "",
    "Ti unit" => "yr",
)
cube4D = DimensionalData.rebuild(
    areacello_ds["areacello"];
    data = [ℰ_ensemblemean;;;; ℰ_ensemblerange],
    dims = (dims(ℰ_ensemblemean)..., dims(DimArray(ones(2), Dim{:statistic}(["ensemble mean", "ensemble range"])))[1]),
    metadata = metadata,
)
arrays = Dict(:seqeff => cube4D, :lat => readcubedata(volcello_ds.lat), :lon => readcubedata(volcello_ds.lon))
ds = Dataset(; properties = metadata, arrays...)

# Save to netCDF file
outputfile = joinpath(inputdir, "Pasquier_etal_GRL_2025_Fig1_data.nc")
@info "Saving mean sequestration efficiency as netCDF file:\n  $(outputfile)"
# ds_chunked = setchunks(ds, (x = 60, y = 60, Ti = length(ds.Ti)))
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)
