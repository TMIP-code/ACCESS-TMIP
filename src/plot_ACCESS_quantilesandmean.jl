# qsub -I -P xv83 -q express -l mem=47GB -l storage=scratch/gh0+scratch/xv83 -l walltime=01:00:00 -l ncpus=12
# This is Fig. 2 in Pasquier et al. (GRL, 2025)

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

members_mean = ["r$(r)i1p1f1" for r in 1:40]
members = ["r$(r)i1p1f1" for r in 1:40]

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
κVdeep_str = "kVdeep" * format(κVdeep, conversion="e")
κVML_str = "kVML" * format(κVML, conversion="e")
κH_str = "kH" * format(κH, conversion="d")
upwind = false
upwind_str = upwind ? "" : "_centered"
upwind_str2 = upwind ? "upwind" : "centered"

# Use yearly time-stepped simulations or monthly ones?
yearly = false
yearly_str = yearly ? "_yearly" : ""
yearly_str2 = yearly ? "(yearly)" : ""

# little helper function to get the year given a quantile
function yearatquantile(ℰ, ℰlevel)
    isnan(ℰ[1]) && return NaN
    out = findfirst(ℰ .< ℰlevel)
    isnothing(out) ? maximum(years) : Float64(out)
end

# To avoid loading and carrying 100s of GB of data around,
# preprocess each member before and only save the 2D data needed for plots.
varname = yearly ? "seqeff" : "calE"
ℰ_file0 = "/scratch/xv83/TMIP/data/$model/$experiment/$(first(members))/$(time_window)/$(varname)$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str).nc"
ℰ_ds0 = open_dataset(ℰ_file0)
years = ℰ_ds0.Ti |> Array
τℰ1050_ensemble = reduce((x,y) -> cat(x, y, dims = 4), map(members) do member
    @info "loading $member ℰ"
    ℰ_file = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)/$(varname)$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str).nc"
    ℰ_ds = open_dataset(ℰ_file)
    ℰ = readcubedata(ℰ_ds[varname])
    τℰ10 = map(
        ts -> yearatquantile(ts, 0.9),
        view(ℰ, i, j, :) for i in 1:size(ℰ,1), j in 1:size(ℰ,2)
    )
    τℰ50 = map(
        ts -> yearatquantile(ts, 0.5),
        view(ℰ, i, j, :) for i in 1:size(ℰ,1), j in 1:size(ℰ,2)
    )
    [τℰ10;;; τℰ50]
end)
τℰ1050_ensemblemean = dropdims(mean(τℰ1050_ensemble, dims = 4), dims = 4)
τℰ1050_ensemblemax = dropdims(maximum(τℰ1050_ensemble, dims = 4), dims = 4)
τℰ1050_ensemblemin = dropdims(minimum(τℰ1050_ensemble, dims = 4), dims = 4)
τℰ1050_ensemblerange = τℰ1050_ensemblemax - τℰ1050_ensemblemin

include("plotting_functions.jl") # load seafloorvalue function

Γout_ensemble = reduce((x, y) -> cat(x, y, dims = 3), map(members) do member
    @info "loading $member Γ†"
    Γout_file = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)/cyclomonth/reemergence_time$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str).nc"
    Γoutyr4D_ds = open_dataset(Γout_file)
    Γoutyr3D = dropdims(mean(readcubedata(Γoutyr4D_ds.adjointage), dims = Ti), dims = Ti)
    seafloorvalue(Γoutyr3D, wet3D, gridmetrics)
end)

Γout_ensemblemean = dropdims(mean(Γout_ensemble, dims = 3), dims = 3)
Γout_ensemblemax = dropdims(maximum(Γout_ensemble, dims = 3), dims = 3)
Γout_ensemblemin = dropdims(minimum(Γout_ensemble, dims = 3), dims = 3)
Γout_ensemblerange = Γout_ensemblemax - Γout_ensemblemin




usecontourf = false

axs = Array{Any,2}(undef, (3, 2))
contours = Array{Any,2}(undef, (3, 2))
nrows, ncols = size(axs)

fig = Figure(size = (ncols * 500, nrows * 250 + 100), fontsize = 18)

yticks = -60:30:60
xticks = -120:60:120 + 360

datamean = (Γout_ensemblemean, τℰ1050_ensemblemean[:,:,2] .|> Float64, τℰ1050_ensemblemean[:,:,1] .|> Float64)
datarange = (Γout_ensemblerange, τℰ1050_ensemblerange[:,:,2] .|> Float64, τℰ1050_ensemblerange[:,:,1] .|> Float64)
𝒓 = rich("r", font = :bold_italic)
Γstr = rich("Γ", superscript("†"), rich("‾", offset = (-0.55, 0.25)), rich("‾", offset = (-0.85, 0.25)))
Γfun = rich(Γstr, rich("(", 𝒓, ")", offset = (0.4, 0)))
# Qstr = rich("Q", font = :italic)
ℰstr = rich("ℰ", rich("‾", offset = (-0.5, 0.15)))
# Q10 = rich(Qstr, "(0.1)")
# Q50 = rich(Qstr, "(0.5)")
# rowlabels = (rich("Mean time, ", Γfun), rich("Median time (", ℰstr, " = 50 %)"), rich("10th percentile time (", ℰstr, " = 90 %)"))
rowlabels = (rich("Mean time, ", Γfun), rich("Median time (", ℰstr, " = 50 %)"), rich("10th %ile time (", ℰstr, " = 90 %)"))
# rowlabels = (rich("Mean time, ", Γfun), rich("Median time (", ℰstr, " = 50 %)"), rich("10 % time (", ℰstr, " = 90 %)"))


for (irow, (x2Dmean, x2Drange, text)) in enumerate(zip(datamean, datarange, rowlabels))

    # Plot ensemble mean
    icol = 1
    levels = 0:200:3000
    colormap = cgrad(:viridis, length(levels), categorical = true)
    highclip = colormap[end]
    colormap = cgrad(colormap[1:end-1], categorical = true)
    colorrange = extrema(levels)

    axs[irow, icol] = ax = Axis(fig[irow, icol]; yticks, xticks, xtickformat, ytickformat, aspect = DataAspect())

    contours[irow, icol] = if usecontourf
        plotcontourfmap!(ax, x2Dmean, gridmetrics; levels, colormap)
    else
        plotmap!(ax, x2Dmean, gridmetrics; colorrange, colormap, highclip) # <- need to fix wrapping longitude for contour levels
    end

    myhidexdecorations!(ax, irow < nrows)
    myhideydecorations!(ax, icol > 1)

    # Plot ensemble range
    icol = 2
    levels = 0:100:1000
    colormap = cgrad(:amp, length(levels), categorical = true)
    highclip = colormap[end]
    colormap = cgrad(colormap[1:end-1], categorical = true)
    colorrange = extrema(levels)

    axs[irow, icol] = ax = Axis(fig[irow, icol]; yticks, xticks, xtickformat, ytickformat, aspect = DataAspect())

    contours[irow, icol] = if usecontourf
        plotcontourfmap!(ax, x2Drange, gridmetrics; levels, colormap)
    else
        plotmap!(ax, x2Drange, gridmetrics; colorrange, colormap, highclip) # <- need to fix wrapping longitude for contour levels
    end

    myhidexdecorations!(ax, irow < nrows)
    myhideydecorations!(ax, icol > 1)

    Label(fig[irow, 0]; text, rotation = π/2, tellheight = false)

end


label = rich("ensemble mean $(time_window[4:7])s characteristic timescales (years)")
cb = Colorbar(fig[nrows + 1, 1], contours[1, 1]; label, vertical = false, flipaxis = false, ticks = 0:1000:3000)
cb.width = Relative(2/3)

label = rich("ensemble range, max − min (years)")
cb = Colorbar(fig[nrows + 1, 2], contours[1, 2]; label, vertical = false, flipaxis = false, ticks = 0:200:1000)
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

# Label(fig[0, 1:2]; text = "$(time_window[4:7])s Seafloor Reemergence Time ($(length(members)) members)", fontsize = 24, tellwidth = false)
# Label(fig[0, 1:2]; text = "$(time_window[4:7])s Characteristic Timescales of Reemergence ($(length(members)) members)$(yearly_str2)", fontsize = 24, tellwidth = false)
rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)

colsize!(fig.layout, 1, Aspect(1, 2.0))
colsize!(fig.layout, 2, Aspect(1, 2.0))
resize_to_layout!(fig)

# save plot
suffix = usecontourf ? "_ctrf" : ""


outputfile = joinpath(outputdir, "reemergencetime$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str)_$(time_window)$(suffix).png")
@info "Saving reemergencetime on sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)
outputfile = joinpath(outputdir, "reemergencetime$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)$(yearly_str)_$(time_window)$(suffix).pdf")
@info "Saving reemergencetime on sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)



# datamean = (Γout_ensemblemean, τℰ50_ensemblemean, τℰ10_ensemblemean)
ikeep = .!isnan.(datamean[1]) .& (seafloorvalue(Z3D, wet3D) .> 3000)
data = datamean[2][ikeep] ./ datamean[1][ikeep] .- 1
weights = Weights(gridmetrics.area2D[ikeep])
fig, ax, plt = hist(data; weights, bins = -1:0.05:0)
outputfile = joinpath(outputdir, "median_over_mean_histogram_$(time_window)$(suffix).png")
@info "Saving reemergencetime on sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)
@show mean(data, weights)
@show std(data, weights)
@show quantile(data, weights, 0:0.1:1)


ikeep = .!isnan.(datamean[1]) .& (seafloorvalue(Z3D, wet3D) .> 3000)
data = datamean[3][ikeep] ./ datamean[2][ikeep] .- 1
weights = Weights(gridmetrics.area2D[ikeep])
fig, ax, plt = hist(data; weights, bins = -1:0.05:1)
outputfile = joinpath(outputdir, "tenthpercentile_over_median_histogram_$(time_window)$(suffix).png")
@info "Saving reemergencetime on sea floor as image file:\n  $(outputfile)"
@show save(outputfile, fig)
@show mean(data, weights)
@show quantile(data, weights, 0:0.1:1)



ikeep = .!isnan.(datamean[1]) .& (seafloorvalue(Z3D, wet3D) .> 3000)
data = datamean[1][ikeep]
weights = Weights(gridmetrics.area2D[ikeep])
fig, ax, plt = hist(data; weights, bins = 0:100:3000)
outputfile = joinpath(outputdir, "mean_histogram_$(time_window)$(suffix).png")
@info "Saving reemergencetime on sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)
@show mean(data, weights)
@show std(data, weights)
@show quantile(data, weights, 0:0.1:1)

ikeep = .!isnan.(datamean[1]) .& (seafloorvalue(Z3D, wet3D) .> 3000)
data = datamean[2][ikeep]
weights = Weights(gridmetrics.area2D[ikeep])
fig, ax, plt = hist(data; weights, bins = 0:100:3000)
outputfile = joinpath(outputdir, "median_histogram_$(time_window)$(suffix).png")
@info "Saving reemergencetime on sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)
mean(data, weights)
quantile(data, weights, 0:0.1:1)

ikeep = .!isnan.(datamean[1]) .& (seafloorvalue(Z3D, wet3D) .> 3000)
data = datamean[3][ikeep]
weights = Weights(gridmetrics.area2D[ikeep])
fig, ax, plt = hist(data; weights, bins = 0:100:3000)
outputfile = joinpath(outputdir, "tenthpercentile_histogram_$(time_window)$(suffix).png")
@info "Saving as image file:\n  $(outputfile)"
save(outputfile, fig)
@show mean(data, weights)
@show std(data, weights)
@show quantile(data, weights, 0:0.1:1)




ikeep = .!isnan.(datamean[1]) .& (seafloorvalue(Z3D, wet3D) .> 3000)
data = datarange[1][ikeep] ./ datamean[1][ikeep]
weights = Weights(gridmetrics.area2D[ikeep])
fig, ax, plt = hist(data; weights, bins = 0:0.1:2)
outputfile = joinpath(outputdir, "rangemeanage_over_meanage_$(time_window)$(suffix).png")
@info "Saving as image file:\n  $(outputfile)"
save(outputfile, fig)
@show mean(data, weights)
@show std(data, weights)
@show quantile(data, weights, 0:0.1:1)



# Save the data to be uploaded with paper
metadata = Dict(
    "description" => "Characteristic timescales as plotted in Fig. 2 in Pasquier et al. (2025)",
    "model" => model,
    "experiment" => experiment,
    "time window" => time_window,
    "unit" => "years",
)
cube4D = DimensionalData.rebuild(areacello_ds["areacello"];
    data = [cat(datamean..., dims=3);;;; cat(datarange..., dims=3)],
    dims = (
        dims(readcubedata(volcello_ds.lat))...,
        dims(DimArray(ones(3), Dim{:timescale}(["mean", "median", "10th percentile"])))[1],
        dims(DimArray(ones(2), Dim{:statistic}(["ensemble mean", "ensemble range"])))[1],
    ),
    metadata = metadata,
)
arrays = Dict(:timescales => cube4D, :lat => readcubedata(volcello_ds.lat), :lon => readcubedata(volcello_ds.lon))
ds = Dataset(; properties = metadata, arrays...)

# Save to netCDF file
outputfile = joinpath(inputdir, "Pasquier_etal_GRL_2025_Fig2_data.nc")
@info "Saving characteristic timescales as netCDF file:\n  $(outputfile)"
# ds_chunked = setchunks(ds, (x = 60, y = 60, Ti = length(ds.Ti)))
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)


