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

# # Gadi directory for input files
# # inputdirfun(member) = "/scratch/xv83/TMIP/data/$model/$experiment/all members/$(time_window)"
# inputdir(member) = "/scratch/xv83/TMIP/data/$model/$experiment/$(member)/$(time_window)/cyclomonth"
# outputdir = "/scratch/xv83/TMIP/data/$model/$experiment/all_members/$(time_window)/cyclomonth"
# mkpath(outputdir)

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






# # Matrix cyclo ages for varied members
# @info "Loading age computed from matrices of different members"
# κVML = 1e-7   # m^2/s
# κVdeep = 1e-7 # m^2/s
# κH = 5        # m^2/s
# κVdeep_str = "kVdeep" * format(κVdeep, conversion="e")
# κH_str = "kH" * format(κH, conversion="d")
# κVML_str = "kVML" * format(κVML, conversion="e")
# members = ["r$(i)i1p1f1" for i in 1:40]
# model_data = map(members) do member
#     f = "ideal_mean_age_$(κVdeep_str)_$(κH_str)_$(κVML_str).nc"
#     age_ds = open_dataset(joinpath(inputdir(member), f))
#     # age3D = (age_ds.age[Ti=1].data + age_ds.age[Ti=12].data) / 2
#     # age3D = age_ds.age[Ti=1].data
#     age3D = dropdims(mean(age_ds.age, dims=:Ti), dims=:Ti).data
#     age3D[wet3D]
# end

# @info "Loading Anderson Acceleration age"
# # Anderson Accelerated age
# AAdir = "/scratch/xv83/bp3051/access-esm/archive/andersonacceleration_test-n10-5415f621/age_output"
# AAfile = joinpath(AAdir, "/ocean_age.res_0035.nc")
# obs_ds = open_dataset(AAfile)
# obs_data = obs_ds.age_global[Time=1].data[wet3D]

# iters = 0:35
# obs_data2 = map(iters) do iter
#     AAfile = joinpath(AAdir, "ocean_age.res_$(format(iter; width = 4, zeropadding = true)).nc")
#     obs_ds = open_dataset(AAfile)
#     obs_data = obs_ds.age_global[Time=1].data[wet3D]
# end


# # Matrix steady ages for varied members
# @info "Loading steady age computed from matrices of different members"
# model_data2 = map(members) do member
#     f = "steady_state_ideal_mean_age_$(κVdeep_str)_$(κH_str)_$(κVML_str).nc"
#     age_ds = open_dataset(joinpath(inputdir(member), f))
#     age3D = age_ds.age.data
#     age3D[wet3D]
# end
# # Mean-flow steady ages for varied members
# @info "Loading steady age computed from mean-flow matrices of different members"
# model_data3 = map(members) do member
#     f = "steady_state_ideal_mean_age_$(κVdeep_str)_$(κH_str)_$(κVML_str)_meanflow.nc"
#     age_ds = open_dataset(joinpath(inputdir(member), f))
#     age3D = age_ds.age.data
#     age3D[wet3D]
# end


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
# TDvals2 = [taylordiagramvalues(data, obs_data, w) for data in model_data2]
# TDvals3 = [taylordiagramvalues(data, obs_data, w) for data in model_data3]
# TDvalsobs = [taylordiagramvalues(data, obs_data, w) for data in obs_data2]






# Do the actual plotting now
# First, construct the figure and a polar axis on the first quadrant
fig = Figure(size=(400, 400))

κVdeep_str = rich("κVdeep = ", format(κVdeep, conversion="e"), " m s", superscript("−2"))
κH_str = rich("κH = ", format(κH, conversion="d"), " m s", superscript("−2"))
κVML_str = rich("κVML = ", format(κVML, conversion="e"), " m s", superscript("−2"))
title = rich(κH_str, ", ", κVdeep_str, ", ", κVML_str)

# Corrticks for Taylor diagram
# corrticks = [-1; -0.99; -0.95; -0.9:0.1:-0.7; -0.6:0.2:0.6; 0.7:0.1:0.9; 0.95; 0.99; 1.0]
corrticks = [0:0.2:0.6; 0.7:0.1:0.9; 0.95; 0.99; 1.0]
function myformat(corrtick)
    isinteger(corrtick) && (corrtick = Int(corrtick))
    str = string(corrtick)
    return replace(str, "-" => "−")
end
rtickformat(rs) = map(r -> "$(format(round(10r/σr)/10, stripzeros = true))σᵣ", rs)
ax = PolarAxis(fig[1, 1];
    thetalimits = (0, π/2), # first quadrant only
    thetagridcolor = (:black, 0.5),
    thetagridstyle = :dot,
    thetaticks = (acos.(corrticks), myformat.(corrticks)),
    # thetaminorticks,
    rlimits = (0, σmax),
    rgridcolor = cgrad(:Archambault, categorical = true)[1],
    rticklabelcolor = cgrad(:Archambault, categorical = true)[1],
    rgridstyle = :dash,
    rticks = (0:σr/2:σmax),
    rtickformat,
)

# Create isolines of root centered mean squared difference (E′) and labels
levels = (0.5:0.5:4) .* σr
rgrid = (0:0.01:1) .* σmax
θgrid = (0:0.01:π)
E′fun(σf, σr, R) = sqrt(σf^2 + σr^2 - 2 * σf * σr * R)
# E′grid = [sqrt(r^2 + σr^2 - 2 * σr * r * cos(θ)) for θ in θgrid, r in rgrid]
E′grid = [E′fun(r, σr, cos(θ)) for θ in θgrid, r in rgrid]
# labelformatter(E′s) = map(E′ -> rich("$(E′/σr)", rich(" σ", subscript("ref"))), E′s)
labelformatter(E′s) = map(E′ -> "$(format(round(10E′/σr)/10, stripzeros = true))σᵣ", E′s)
contour!(ax, θgrid, rgrid, E′grid;
    levels,
    labels = true,
    labelformatter,
    color = cgrad(:Archambault, categorical = true)[3]
)

# # skill score isolines
# R₀ = 1 # maximum correlation obtainable. Don't think I'll need this but may be useful
# S(σf, σr, R) = 4 * (1 + R) / ((σf/σr + σr/σf)^2 * (1 + R₀))
# # S(σf, σr, R) = 4 * (1 + R)^4 / ((σf/σr + σr/σf)^2 * (1 + R₀)^4)
# Sgrid = [S(r, σr, cos(θ)) for θ in θgrid, r in rgrid]
# Slevels = [0:0.1:0.9; 0.95; 0.99]
# contour!(ax, θgrid, rgrid, Sgrid;
#     levels = Slevels,
#     labels = true,
#     color = cgrad(:Archambault, categorical = true)[4]
# )

# Now, plot the actual data
σfs = [vals.σf for vals in TDvals]
Rs = [vals.R for vals in TDvals]
σfs2 = [vals.σf for vals in TDvals2]
Rs2 = [vals.R for vals in TDvals2]
σfs3 = [vals.σf for vals in TDvals3]
Rs3 = [vals.R for vals in TDvals3]

# Transform data to Cartesian space?
"A transformation function that goes from correlation and standard deviation to the Taylor plot's Cartesian space."
xy_from_R_and_σ(R, σ) = Point2(σ * R, sqrt(σ^2 - (σ * R)^2))

# Plot all matrix ages
transformation = Transformation(ax.scene.transformation; transform_func = identity)
x, y = collect.(zip(xy_from_R_and_σ.(Rs, σfs)...) |> collect)
scatter!(ax, x, y;
    color = :red,
    markersize = 3,
    transformation,
)
offset = 100
txtline = [offset, offset/5]
fontsize = 10
lines!(ax, mean(x) .- txtline, mean(y) .+ txtline / 2; linewidth = 1, color = :black, transformation)
text!(ax, mean(x) - offset, mean(y) + offset / 2; text = "cyclostationary\nage", transformation, align = (:right, :bottom), offset = (0, 0), fontsize)

offset = 100
txtline = [offset, offset/5]
x2, y2 = collect.(zip(xy_from_R_and_σ.(Rs2, σfs2)...) |> collect)
scatter!(ax, x2, y2;
    color = :blue,
    markersize = 3,
    transformation,
)
lines!(ax, mean(x2) .+ txtline / 2, mean(y2) .+ txtline; linewidth = 1, color = :black, transformation)
text!(ax, mean(x2) + offset / 2, mean(y2) + offset; text = "steady age\n(mean matrices)", transformation, align = (:center, :bottom), offset = (0, 0), fontsize)


offset = 50
txtline = [offset, offset/2.5]
x3, y3 = collect.(zip(xy_from_R_and_σ.(Rs3, σfs3)...) |> collect)
scatter!(ax, x3, y3;
    color = :purple,
    markersize = 3,
    transformation,
)
lines!(ax, mean(x3) .+ txtline, mean(y3) .- txtline; linewidth = 1, color = :black, transformation)
text!(ax, mean(x3) + offset, mean(y3) - offset; text = "steady age\n(mean flow)", transformation, align = (:left, :top), offset = (0, 0), fontsize)
# text!(ax, x[end], y[end]; text = "max", transformation, align = (:center, :center), fontsize)
# scatter!(ax, x, y;
#     color = :red,
#     marker = :cross,
#     transformation = Transformation(ax.scene.transformation; transform_func = identity)
# )

# Add lines for AA Interation
σfsobs = [vals.σf for vals in TDvalsobs]
Rsobs = [vals.R for vals in TDvalsobs]
xobs, yobs = collect.(zip(xy_from_R_and_σ.(Rsobs, σfsobs)...) |> collect)
scatterlines!(ax, xobs, yobs;
    color = :black,
    linewidth = 1,
    markersize = 3,
    transformation,
)
text!(ax, xobs[1], yobs[1]; text = "AA sequence", transformation, align = (:left, :bottom), offset = (0, 3), fontsize)
lines!(ax, xobs[end] .+ txtline, yobs[end] .+ txtline; linewidth = 1, color = :black, transformation)
text!(ax, xobs[end] + offset, yobs[end] + offset; text = "AA age", transformation, align = (:left, :bottom), offset = (0, 0), fontsize)

# Add text for kappas
σκ = 1.5σr
Rκ = 0.4
xκ, yκ = xy_from_R_and_σ(Rκ, σκ)
text!(ax, xκ - 0.2σr, yκ - 0.2σr; text = κH_str, transformation, align = (:left, :bottom), fontsize)
text!(ax, xκ - 0.2σr, yκ - 0.2σr - offset; text = κVdeep_str, transformation, align = (:left, :bottom), fontsize)
text!(ax, xκ - 0.2σr, yκ - 0.2σr - 2offset; text = κVML_str, transformation, align = (:left, :bottom), fontsize)


# Plot reference (AA age)
scatter!(ax, Point2(xy_from_R_and_σ(1, σr));
    color = :black,
    transformation,
)

outputfile = joinpath(outputdir, "Taylor_diagram_members.png")
@info "Saving seqeff on sea floor as image file:\n  $(outputfile)"
save(outputfile, fig)


# bins = 0:50:2000
# # fig, ax, plt = hist(obs_data[Z3D[wet3D] .> 500]; label = "obs", bins, color = (:black, 0.2))
# # hist!(ax, model_data[1][Z3D[wet3D] .> 500]; label = "min", bins, color = (:red, 0.2))
# # hist!(ax, model_data[end][Z3D[wet3D] .> 500]; label = "max", bins, color = (:blue, 0.2))
# fig, ax, plt = density(obs_data[Z3D[wet3D] .> 500]; label = "obs", color = (:black, 0.2))
# density!(ax, model_data[1][Z3D[wet3D] .> 500]; label = "min", color = (:red, 0.2))
# density!(ax, model_data[end][Z3D[wet3D] .> 500]; label = "max", color = (:blue, 0.2))
# # axislegend(ax)
# outputfile = joinpath(outputdir, "hist_kVML.png")
# @info "Saving seqeff on sea floor as image file:\n  $(outputfile)"
# save(outputfile, fig)
