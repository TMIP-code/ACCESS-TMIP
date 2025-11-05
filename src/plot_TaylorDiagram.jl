# # qsub -I -P xv83 -l mem=64GB -l storage=scratch/gh0+scratch/xv83 -l walltime=02:00:00 -l ncpus=12

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
# # using LaTeXStrings
# using Format

# model = "ACCESS-ESM1-5"

# time_window = "Jan1850-Dec1859"
# experiment = "historical"
# member = "r20i1p1f1"

# # Gadi directory for input files
# # inputdirfun(member) = "/scratch/xv83/TMIP/data/$model/$experiment/all members/$(time_window)"
# inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/$(member)/$(time_window)/cyclomonth"
# outputdir = inputdir

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


# # Matrix ages for varied diffusivities
# @info "Loading age computed from matrices with different diffusivities"
# κVdeeps = [1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4]
# κHs = [5, 15, 50, 150, 500, 1500, 5000]
# κs = Iterators.product(κVdeeps, κHs)
# model_data = map(κs) do κ
#     κVdeep, κH = κ
#     κVdeep_str = "kVdeep" * format(κVdeep, conversion="e")
#     κH_str = "kH" * format(κH, conversion="d")
#     f = "ideal_mean_age_$(κVdeep_str)_$(κH_str).nc"
#     age_ds = open_dataset(joinpath(inputdir, f))
#     age3D = (age_ds.age[Ti=1].data + age_ds.age[Ti=12].data) / 2
#     age3D[wet3D]
# end


# @info "Loading Anderson Acceleration age"
# # Anderson Accelerated age
# AAfile = "/scratch/xv83/bp3051/access-esm/archive/andersonacceleration_test-n10-5415f621/age_output/ocean_age.res_0035.nc"
# obs_ds = open_dataset(AAfile)
# obs_data = obs_ds.age_global[Time=1].data[wet3D]


# # Taylor diagram function that returns all the required values
# # notation taken from the original Taylor paper
# # TODO check that the identity holds when using weights!
# function taylordiagramvalues(f, r, args...)

#     # STDs and means
#     σf = std(f, args...; corrected = false)
#     σr = std(r, args...; corrected = false)
#     f̄ = mean(f, args...)
#     r̄ = mean(r, args...)

#     # Correlation coefficient
#     R = cor([f r], args...)[2]

#     # Root Mean Square Difference
#     E = sqrt(mean((f .- r) .^ 2, args...))

#     # Bias
#     Ē = f̄ - r̄

#     # Centered Root Mean Square Difference
#     E′ = sqrt(mean(((f .- f̄) - (r .- r̄)) .^ 2, args...))

#     # Full Mean Square Difference
#     E² = E′^2 + Ē^2

#     # Normalized values (maybe that needs to be a kwarg)
#     Ê′ = E′ / σr
#     σ̂f = σf / σr
#     σ̂r = 1.0

#     return (; σr, σf, R, E, Ē, E′, E², Ê′, σ̂f, σ̂r)
# end

# # Calculate the Taylor diagram values
# w = weights(v3D[wet3D])
# TDvals = [taylordiagramvalues(data, obs_data, w) for data in model_data]
# σr = TDvals[1].σr
# σmax = 2TDvals[1].σr

# Do the actual plotting now
# First, construct the figure and a polar axis on the first quadrant
fig = Figure(size = (800, 400))

# Corrticks for Taylor diagram
# corrticks = [-1; -0.99; -0.95; -0.9:0.1:-0.7; -0.6:0.2:0.6; 0.7:0.1:0.9; 0.95; 0.99; 1.0]
corrticks = [0:0.2:0.6; 0.7:0.1:0.9; 0.95; 0.99; 1.0]
function myformat(corrtick)
    isinteger(corrtick) && (corrtick = Int(corrtick))
    str = string(corrtick)
    return replace(str, "-" => "−")
end

ax = PolarAxis(
    fig[1, 1];
    thetalimits = (0, π / 2), # first quadrant only
    thetagridcolor = (:black, 0.5),
    thetagridstyle = :dot,
    thetaticks = (acos.(corrticks), myformat.(corrticks)),
    # thetaminorticks,
    rlimits = (0, σmax),
    rgridcolor = cgrad(:Archambault, categorical = true)[1],
    rticklabelcolor = cgrad(:Archambault, categorical = true)[1],
    rgridstyle = :dash,
    rticks = 0:σr:σmax,
)

# Create isolines of root centered mean squared difference (E′) and labels
levels = (0.5:0.5:4) .* σr
rgrid = (0:0.01:1) .* σmax
θgrid = (0:0.01:π)
E′fun(σf, σr, R) = sqrt(σf^2 + σr^2 - 2 * σf * σr * R)
# E′grid = [sqrt(r^2 + σr^2 - 2 * σr * r * cos(θ)) for θ in θgrid, r in rgrid]
E′grid = [E′fun(r, σr, cos(θ)) for θ in θgrid, r in rgrid]
# labelformatter(E′s) = map(E′ -> rich("$(E′/σr)", rich(" σ", subscript("ref"))), E′s)
labelformatter(E′s) = map(E′ -> "$(format(round(10E′ / σr) / 10, stripzeros = true))σᵣ", E′s)
contour!(
    ax, θgrid, rgrid, E′grid;
    levels,
    labels = true,
    labelformatter,
    color = cgrad(:Archambault, categorical = true)[3]
)

# skill score isolines
R₀ = 0.9940801590730742 # maximum correlation obtainable from ensemble
S(σf, σr, R) = 4 * (1 + R) / ((σf / σr + σr / σf)^2 * (1 + R₀))
# S(σf, σr, R) = 4 * (1 + R)^4 / ((σf/σr + σr/σf)^2 * (1 + R₀)^4)
Sgrid = [S(r, σr, cos(θ)) for θ in θgrid, r in rgrid]
Slevels = [0:0.1:0.9; 0.95; 0.99]
contour!(
    ax, θgrid, rgrid, Sgrid;
    levels = Slevels,
    labels = true,
    color = cgrad(:Archambault, categorical = true)[4]
)

# Now, plot the actual data
σfs = [vals.σf for vals in TDvals]
Rs = [vals.R for vals in TDvals]

# Transform data to Cartesian space?
"A transformation function that goes from correlation and standard deviation to the Taylor plot's Cartesian space."
xy_from_R_and_σ(R, σ) = Point2(σ * R, sqrt(σ^2 - (σ * R)^2))
Ps = xy_from_R_and_σ.(Rs, σfs)
x, y = collect.(zip(xy_from_R_and_σ.(Rs, σfs)...) |> collect)
X = reshape(x, length(κVdeeps), length(κHs)) # <- not sure that was necessary
Y = reshape(y, length(κVdeeps), length(κHs))
# Above I used R = 1 and σr for the reference point
# Plot all matrix ages
Xcol = reduce(vcat, [col; NaN] for col in eachcol(X))
Ycol = reduce(vcat, [col; NaN] for col in eachcol(Y))
Xrow = reduce(vcat, [row; NaN] for row in eachrow(X))
Yrow = reduce(vcat, [row; NaN] for row in eachrow(Y))
transformation = Transformation(ax.scene.transformation; transform_func = identity)
scatterlines!(
    ax, [Ps[1, :]; Ps[:, end]; Ps[end, end:-1:1]; Ps[end:-1:1, 1]];
    color = :blue,
    markersize = 3,
    linewidth = 1,
    transformation,
)
text!(ax, Ps[1, 1]; text = "min", transformation, align = (:center, :center), fontsize = 6, offset = (0, 5))
text!(ax, Ps[end, 1]; text = "max Vdeep", transformation, align = (:center, :center), fontsize = 6, offset = (0, 5))
text!(ax, Ps[end, end]; text = "max", transformation, align = (:center, :center), fontsize = 6, offset = (0, 5))
text!(ax, Ps[1, end]; text = "max H", transformation, align = (:center, :center), fontsize = 6, offset = (0, 5))
# scatter!(ax, x[:], y[:];
#     color = [TDval.Ē for TDval in TDvals][:],
#     colorrange = (-200, 200),
#     colormap = :berlin,
#     marker = :cross,
#     transformation = Transformation(ax.scene.transformation; transform_func = identity)
# )

# Plot reference (AA age)
scatter!(
    ax, Point2(xy_from_R_and_σ(1, σr));
    color = :black,
    transformation,
)

outputfile = joinpath(outputdir, "Taylor_diagram.png")
@info "Saving image file:\n  $(outputfile)"
save(outputfile, fig)
