
# using Pkg
# Pkg.activate(".")
# Pkg.instantiate()


# ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
# using OceanTransportMatrixBuilder
# using NetCDF
# using YAXArrays
# using DataFrames
# using DimensionalData
# using SparseArrays
# using LinearAlgebra
# using Unitful
# using Unitful: s, yr, d
# using Statistics
# using Format
# using Dates
# using FileIO
# using ProgressMeter
# using GLMakie
# using GeoMakie
# using OceanBasins
# using NaNStatistics



# include("plotting_functions.jl")


# # script options
# @show model = "ACCESS-ESM1-5"
# if isempty(ARGS)
#     member = "r1i1p1f1"
#     # experiment = "historical"
#     # time_window = "Jan1850-Dec1859"
#     # time_window = "Jan1990-Dec1999"
#     experiment = "ssp370"
#     time_window = "Jan2030-Dec2039"
#     # time_window = "Jan2090-Dec2099"
#     WRITEDATA = "true"
# else
#     experiment, member, time_window, finalmonth, WRITEDATA = ARGS
# end
# WRITEDATA = parse(Bool, WRITEDATA)
# @show experiment
# @show member
# @show time_window


# # Load areacello and volcello for grid geometry
# fixedvarsinputdir = "/Users/benoitpasquier/Data/TMIP/data/$model"
# volcello_ds = open_dataset(joinpath(fixedvarsinputdir, "volcello.nc"))
# areacello_ds = open_dataset(joinpath(fixedvarsinputdir, "areacello.nc"))

# # Gadi directory for input files
# lumpby = "month"
# inputdir = "/Users/benoitpasquier/Data/TMIP/data/$model/$experiment/$member/$(time_window)"
# cycloinputdir = joinpath(inputdir, "cyclo$lumpby")



# # Load fixed variables in memory
# areacello = readcubedata(areacello_ds.areacello)
# volcello = readcubedata(volcello_ds.volcello)
# lon = readcubedata(volcello_ds.lon)
# lon2 = mod.(gridmetrics.lon .- 80, 360) .+ 80
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
# (; lon_vertices, lat_vertices, v3D) = gridmetrics

# # Make indices
# indices = makeindices(v3D)
# (; N, wet3D) = indices

# # Load netCDF file
# inputfile = joinpath(inputdir, "calE.nc")
# @info "Loading mean sequestration efficiency as netCDF file:\n  $(inputfile)"
# ds = open_dataset(inputfile)
# ℰ = readcubedata(ds.calE).data

# years = ds.Ti |> Array



# Figure
fig = Figure(size = (1200, 800))


sg = SliderGrid(
    fig[1, 1],
    (label = "Year", range = years, format = "{:1d}", startvalue = 100),
    (label = "Efficiency", range = 0:0.1:1, format = "{:.1f}", startvalue = 0.5),
    (label = "Longitude", range = 80:80+360, format = "{:.1f}°E", startvalue = 360),
    (label = "Latitude", range = -90:90, format = "{:.1f}°N", startvalue = 0),
    tellwidth = false,
    tellheight = false)

sliderobservables = [s.value for s in sg.sliders]
sliderP = lift(sliderobservables[3:4]...) do lon, lat
    (lon, lat)
end
sliderij = lift(sliderP) do sliderP
    Tuple(argmin(map(P -> norm(P .- sliderP), zip(lon, lat))))
end
ℰij = lift(sliderij) do sliderij
    i, j = sliderij
    view(ℰ, i, j, :)
end


# time series
ax = Axis(fig[2, 2], xlabel = "years", ylabel = "sequestration efficiency",
    limits = (0, years[end], 0, 1))
Makie.deactivate_interaction!(ax, :rectanglezoom)
lines!(ax, years, ℰij; color = :blue)
yearℰpoint = select_point(ax; marker = :circle, markersize = 10, color = :red)
# vlines!(ax, sliderobservables[1], color = :black, linestyle = :dash)
# hlines!(ax, sliderobservables[2], color = :black, linestyle = :dash)
vlines!(ax, lift(P -> P[1], yearℰpoint), color = :black, linestyle = :dash)
hlines!(ax, lift(P -> P[2], yearℰpoint), color = :black, linestyle = :dash)
scatter!(ax, yearℰpoint)
# scatter!(ax, yearℰpoint; marker = :circle, markersize = 10, color = :red)


# map of sequestration efficiency gievn a year
gb = fig[2, 1] = GridLayout()
axb = Axis(gb[1, 1], xlabel = "lon", ylabel = "lat",
    limits = (80, 80 + 360, -90, 90))
Makie.deactivate_interaction!(axb, :rectanglezoom)
viewatyear(ℰ, year) = view(ℰ, :, :, year)
gbdata = lift(sliderobservables[1]) do year
    viewatyear(ℰ, year)
end
imgb = surface!(axb, lon2, lat.data, gbdata, colormap = :viridis, colorrange = (0, 1), shading = NoShading)
cb = Colorbar(gb[1, 2], imgb)
scb = scatter!(axb, sliderP, color = :red, markersize = 10, strokecolor = :black)
translate!(scb, 0, 0, 100)

# map of duration given sequestration quantile
gc = fig[1, 2] = GridLayout()
axc = Axis(gc[1, 1], xlabel = "lon", ylabel = "lat",
    limits = (80, 80 + 360, -90, 90))
Makie.deactivate_interaction!(axc, :rectanglezoom)
function yearatquantile(ℰ, ℰlevel)
    isnan(ℰ[1]) && return NaN
    out = findfirst(ℰ .< ℰlevel)
    isnothing(out) ? maximum(years) : Float64(out)
end
gcdata = lift(sliderobservables[2]) do ℰlevel
    map(
        ts -> yearatquantile(ts, ℰlevel),
		view(ℰ, i, j, :) for i in 1:size(ℰ,1), j in 1:size(ℰ,2)
    )
end
imgc = surface!(axc, lon2, lat.data, gcdata, colormap = :magma, colorrange = (0, 1000), shading = NoShading)
scc = scatter!(axc, sliderP, color = :red, markersize = 10, strokecolor = :black)
translate!(scc, 0, 0, 100)
cc = Colorbar(gc[1, 2], imgc, label = "years of sequestration")

fig





