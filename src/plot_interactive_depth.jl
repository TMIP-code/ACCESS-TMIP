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
# using GLMakie
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

# include("plotting_functions.jl")

file = if length(ARGS) < 1
    "/Users/z3319805/Data/TMIP/data/ACCESS-OM2-025/steady_age_kVdeep3e-05_kH75_kVML1e+00.nc"
else
    ARGS[1]
end

ds = open_dataset(file)

lon = ds.lon |> Array
lon2 = mod.(lon .- 80, 360) .+ 80

lat = ds.lat |> Array
age3D = ds.age |> Array
lev = ds.lev |> Array

Z0 = zeros(size(age3D)[1:2])

fig = Figure()

slider_levindex = Slider(fig[1, 2], range = reverse(1:length(lev)), startvalue = 20, horizontal = false)


age2D = lift(slider_levindex.value) do k
    view(age3D, :, :, k)
end

title = @lift("Age at depth $(round(Int, lev[$(slider_levindex.value)]))m")

ax = Axis(fig[1, 1], title = title, xlabel = "Longitude", ylabel = "Latitude")

sf = surface!(ax, lon2, lat, Z0, color = age2D, colorrange = (0, 1500), shading = NoShading, interpolate = true)
cb = Colorbar(fig[2, 1], sf; vertical = false, flipaxis = false)

fig