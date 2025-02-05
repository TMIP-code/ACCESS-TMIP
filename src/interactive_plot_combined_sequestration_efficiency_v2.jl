
using Pkg
Pkg.activate(".")
Pkg.instantiate()


ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
using OceanTransportMatrixBuilder
using NetCDF
using YAXArrays
using DataFrames
using DimensionalData
using SparseArrays
using LinearAlgebra
using Unitful
using Unitful: s, yr, d
using Statistics
using Format
using Dates
using FileIO
using ProgressMeter
using GLMakie
using GeoMakie
using OceanBasins
using NaNStatistics
using GeometryBasics
using LibGEOS
using GeometryOps


include("plotting_functions.jl")


# script options
@show model = "ACCESS-ESM1-5"
if isempty(ARGS)
    member = "r1i1p1f1"
    # experiment = "historical"
    # time_window = "Jan1850-Dec1859"
    # time_window = "Jan1990-Dec1999"
    experiment = "ssp370"
    time_window = "Jan2030-Dec2039"
    # time_window = "Jan2090-Dec2099"
    WRITEDATA = "true"
else
    experiment, member, time_window, finalmonth, WRITEDATA = ARGS
end
WRITEDATA = parse(Bool, WRITEDATA)
@show experiment
@show member
@show time_window


# Load areacello and volcello for grid geometry
fixedvarsinputdir = "/Users/benoitpasquier/Data/TMIP/data/$model"
volcello_ds = open_dataset(joinpath(fixedvarsinputdir, "volcello.nc"))
areacello_ds = open_dataset(joinpath(fixedvarsinputdir, "areacello.nc"))

# Gadi directory for input files
lumpby = "month"
inputdir = "/Users/benoitpasquier/Data/TMIP/data/$model/$experiment/$member/$(time_window)"
cycloinputdir = joinpath(inputdir, "cyclo$lumpby")



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
(; lon_vertices, lat_vertices, v3D) = gridmetrics
lon2 = mod.(gridmetrics.lon .- 80, 360) .+ 80

# Make indices
indices = makeindices(v3D)
(; N, wet3D) = indices

# Load netCDF file
inputfile = joinpath(inputdir, "calE.nc")
@info "Loading mean sequestration efficiency as netCDF file:\n  $(inputfile)"
ds = open_dataset(inputfile)
ℰ = readcubedata(ds.calE).data

years = ds.Ti |> Array



# Figure
fig = Figure(size = (1200, 800))

ℰlevels = 0:0.01:1

ga = fig[1, 1] = GridLayout()
sg = SliderGrid(
    ga[1, 1],
    (label = "Year", range = [y for y in years if mod(y, 10) == 0], format = "{:1d}", startvalue = 100),
    (label = "Efficiency", range = ℰlevels, format = "{:.1f}", startvalue = 0.5),
    (label = "Longitude", range = 80:80+360, format = "{:.1f}°E", startvalue = 80 + 180),
    (label = "Latitude", range = -90:90, format = "{:.1f}°N", startvalue = 0),
    tellwidth = false,
    tellheight = false)

check_injection_timeseries = Observable(false)

subga = ga[2, 1] = GridLayout(tellheight = false, tellwidth = false)
btn = Checkbox(subga[1, 1], checked = false, tellheight = false)
Label(subga[1, 2], "Toggle injection time series check", halign = :left)
alpha = @lift($(btn.checked) ? 1.0 : 0)



sliderobservables = [s.value for s in sg.sliders]
sliderP = lift(sliderobservables[3:4]...) do lon, lat
    (lon, lat)
end
sliderij = lift(sliderP) do sliderP
    Tuple(argmin(map(P -> norm(P .- sliderP), zip(lon2, lat))))
end
ℰij = lift(sliderij) do sliderij
    i, j = sliderij
    view(ℰ, i, j, :)
end
ℰk = lift(sliderij) do sliderij
    i, j = sliderij
    k = findlast(view(wet3D, i, j, :))
    isnothing(k) ? 0 : k
end


# time series
titled = @lift("$($(sliderP)): $($(sliderij)) level $($(ℰk))")
ax = Axis(fig[2, 2], xlabel = "years", ylabel = "sequestration efficiency",
    limits = (0, years[end], 0, 1),
    yticks = 0:0.2:1,
    title = titled)
Makie.deactivate_interaction!(ax, :rectanglezoom)
yearℰpoint = select_point(ax; marker = :circle, markersize = 10, color = :red)
on(yearℰpoint) do yearℰpoint
    year, ℰ = yearℰpoint
    set_close_to!(sg.sliders[1], year)
    set_close_to!(sg.sliders[2], ℰ)
end

vlines!(ax, sliderobservables[1], color = :black, linestyle = :dash)
hlines!(ax, sliderobservables[2], color = :black, linestyle = :dash)
# vlines!(ax, lift(P -> P[1], yearℰpoint), color = :red, linestyle = :dash)
# hlines!(ax, lift(P -> P[2], yearℰpoint), color = :red, linestyle = :dash)
# scatter!(ax, sliderobservables[1:2]...; marker = :circle, markersize = 10, color = :red)

# plot the time series
lines!(ax, years, ℰij; color = :blue)

# plot the time series of the injection locations as a check
src_Ps = (
    Karratha = (115.45849390000001,-16.56466979999999), # Carnarvon Basin?" North West of Australia
    Portland = (141.73529860000008,-38.93477809999996), # Otway Basin" South West of Melbourne (West of Tas)
    Marlo = (149.05333500000006, -38.25798499999996), # "Shark 1" Gippsland Basin" South East (East of Tas)
)
map(keys(src_Ps), values(src_Ps)) do srcname, src_P
    src_i, src_j = Tuple(argmin(map(P -> norm(P .- src_P), zip(lon, lat))))
    lines!(ax, years, ℰ[src_i, src_j, :]; alpha, color = :red)
    # Compare to actual data
    datainputdir = joinpath(fixedvarsinputdir, experiment, member, time_window)
    outputfile = joinpath(datainputdir, "$(srcname)_timeseries_injected.jld2")
    @info "Loading injection time series file:\n  $(outputfile)"
    umass = load(outputfile, "umass")
    src_mass = load(outputfile, "src_mass")
    lines!(ax, years[1:length(umass) - 4], umass[5:length(umass)] / src_mass; color = :black, linestyle = :dash, alpha)
end


# map of sequestration efficiency gievn a year
gb = fig[2, 1] = GridLayout()
titleb = @lift("Fraction of C sequestered after $($(sliderobservables[1])) years")
axb = Axis(gb[1, 1], xlabel = "lon", ylabel = "lat",
    limits = (80, 80 + 360, -90, 90),
    title = titleb)
Makie.deactivate_interaction!(axb, :rectanglezoom)
lonlatPb = select_point(axb; marker = :circle, markersize = 10, color = :red)
on(lonlatPb) do lonlatPb
    lon, lat = lonlatPb
    set_close_to!(sg.sliders[3], lon)
    set_close_to!(sg.sliders[4], lat)
end
viewatyear(ℰ, year) = view(ℰ, :, :, year)
gbdata = lift(sliderobservables[1]) do year
    viewatyear(ℰ, year)
end
colormapb = cgrad(:viridis, 10, categorical = true)
imgb = surface!(axb, lon2, lat.data, gbdata, colormap = colormapb, colorrange = (0, 1), shading = NoShading)
cb = Colorbar(gb[1, 2], imgb)
scb = scatter!(axb, sliderP, color = :red, markersize = 10, strokecolor = :black)
translate!(scb, 0, 0, 10)

# map of duration given sequestration quantile
gc = fig[1, 2] = GridLayout()
titlec = @lift("Duration of sequestration at $($(sliderobservables[2])) efficiency")
axc = Axis(gc[1, 1], xlabel = "lon", ylabel = "lat",
    limits = (80, 80 + 360, -90, 90),
    title = titlec)
Makie.deactivate_interaction!(axc, :rectanglezoom)
lonlatPc = select_point(axc; marker = :circle, markersize = 10, color = :red)
on(lonlatPc) do lonlatPc
    lon, lat = lonlatPc
    set_close_to!(sg.sliders[3], lon)
    set_close_to!(sg.sliders[4], lat)
end
function yearatquantile(ℰ, ℰlevel)
    isnan(ℰ[1]) && return NaN
    out = findfirst(ℰ .< ℰlevel)
    isnothing(out) ? maximum(years) : Float64(out)
end
gcdata2 = Array{Float64}(undef, size(ℰ,1), size(ℰ,2), length(ℰlevels))
for (ilevel, ℰlevel) in enumerate(ℰlevels)
    gcdata2[:, :, ilevel] .= map(
        ts -> yearatquantile(ts, ℰlevel),
		view(ℰ, i, j, :) for i in 1:size(ℰ,1), j in 1:size(ℰ,2)
    )
end
gcdata = lift(sg.sliders[2].selected_index) do ilevel
    view(gcdata2, :, :, ilevel)
end
colormapc = cgrad(:magma, 10, categorical = true)
imgc = surface!(axc, lon2, lat.data, gcdata, colormap = colormapc, colorrange = (0, 1000), shading = NoShading)
scc = scatter!(axc, sliderP, color = :red, markersize = 10, strokecolor = :black)
translate!(scc, 0, 0, 1100)
cc = Colorbar(gc[1, 2], imgc, label = "years of sequestration")

fig





