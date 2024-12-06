# qsub -I -P xv83 -q hugemem -l mem=360GB -l storage=scratch/gh0+scratch/xv83 -l walltime=02:00:00 -l ncpus=48

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
using LinearSolve
import Pardiso # import Pardiso instead of using (to avoid name clash?)
using NonlinearSolve
using ProgressMeter
try
    using CairoMakie
catch
    using CairoMakie
end
using GeoMakie
using OceanBasins
using NaNStatistics

include("plotting_functions.jl")

# The tracer equation is
#
#   ∂x(t)/∂t + T(t) x(t) = s(t) - Ω x(t)
#
# where Ω "relaxes" x to zero in the top layer.
#
# Applying Backward Euler time step gives
#
#   (I + δt M(t+δt)) x(t+δt) = x(t) + δt s(t+δt)
#
# where M(t) = T(t) + Ω.
#
# More succintly, if k represents the months (1..12)
#
#   (I + δt Mₖ₊₁) xₖ₊₁ = xₖ + δt sₖ₊₁          (1)
#
# Here δt = δt(k..k+1) is the time that separates
# the "center time" of climatological months k and k+1.
# So the δt that multiplies Mₖ is δ(k-1..k).
#
# Here we switch on s for 10 years at given locations
# and then switch it off to see how long it takes to disappear
# in the top layer.
#
# So if we have s on for the 2030s, we start the simulation with
#
#   x₀ = x(Dec 2029)
#




# Load matrix and grid metrics
@show model = "ACCESS-ESM1-5"
# @show experiment = ARGS[1]
# @show member = ARGS[2]
# @show time_window = ARGS[3]
@show experiment = "historical"
@show time_window = "Jan1850-Dec1859"
# @show time_window = "Jan1990-Dec1999"
# @show experiment = "ssp370"
# @show time_window = "Jan2030-Dec2039"
# @show time_window = "Jan2090-Dec2099"
@show member = "r1i1p1f1"

# Injection locations
# I Identify location with nearest towns to use the name for saving files
src_Ps = (
    Karratha = (115.45849390000001,-16.56466979999999), # Carnarvon Basin?" North West of Australia
    Portland = (141.73529860000008,-38.93477809999996), # Otway Basin" South West of Melbourne (West of Tas)
    Marlo    = (148.74669400000005,-38.6952          ), # Gippsland Basin" South East (East of Tas)
)
Nsrc = length(src_Ps)

lumpby = "month"
steps = 1:12
Nsteps = length(steps)
δt = ustrip(s, 1yr / Nsteps) # TODO maybe use exact mean number of days (more important for monthly because Feb)?


# Gadi directory for input files
inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)"
cycloinputdir = joinpath(inputdir, "cyclo$lumpby")
# Load areacello and volcello for grid geometry
volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))

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
(; lon_vertices, lat_vertices, v3D, ) = gridmetrics

# Make indices
indices = makeindices(v3D)
(; N, wet3D) = indices

issrf = let
    issrf3D = falses(size(wet3D))
    issrf3D[:,:,1] .= true
    issrf3D[wet3D]
end
Ω = sparse(Diagonal(Float64.(issrf)))


# Build matrices
# @time "building Ms" Ms = [
#     begin
#         inputfile = joinpath(cycloinputdir, "cyclo_matrix_$step.jld2")
#         @info "Loading matrix from $inputfile"
#         T = load(inputfile, "T")
#         T + Ω
#     end
#     for step in steps
# ]
@time "building Ms" Ms = map(steps) do step
    inputfile = joinpath(cycloinputdir, "cyclo_matrix_$step.jld2")
    @info "Loading matrix from $inputfile"
    T = load(inputfile, "T")
    T + Ω
end
@time "building mean_days_in_months" mean_days_in_months = map(steps) do step
    inputfile = joinpath(cycloinputdir, "cyclo_matrix_$step.jld2")
    @info "Loading mean_days_in_month from $inputfile"
    load(inputfile, "mean_days_in_month")
end
# So the δt that multiplies Mₖ is δ(k-1..k)
# which is 0.5 of the mean days in months k-1 and k
δts = map(eachindex(mean_days_in_months)) do k
    (mean_days_in_months[mod1(k - 1, 12)] + mean_days_in_months[k]) / 2
end

function initstepprob(A)
    prob = LinearProblem(A, ones(N, Nsrc))
    return init(prob, MKLPardisoIterate(; nprocs = 48), rtol = 1e-10)
end

p = (;
    δts,
    stepprobs = [initstepprob(I + δt * M) for (δt, M) in zip(δts, Ms)]
)

# Point sources for carbon injection
# for each source point we create a vector and we hcat them into a Nocn × Nsrc matrix
src_injections = reduce(hcat, map(src_Ps) do src_P
    src_i, src_j = Tuple(argmin(map(P -> norm(P .- src_P), zip(lon, lat))))
    src_k = findlast(wet3D[src_i, src_j,:])
    isnothing(src_k) && error("No seafloor at (lon,lat)=$src_P")
    src3D = fill(NaN, size(wet3D))
    src3D[wet3D] .= 0
    src3D[src_i, src_j, src_k] = 1
    src3D[wet3D]
end)

src_years = 10
src(y) = (y ≤ src_years) * src_injections
v = v3D[wet3D]
src_mass = v' * sum(src(y) * sum(δts) for y in 1:src_years)

# Plot
function plotandprint!(x3D, state, ksrc, year, text="")
    x3D[wet3D] = state[:, ksrc]
    fig = Figure(size=(600, 300))
    ax = Axis(fig[1, 1], title="Injected tracer (year $year: $text sequestered)")
    x2D = volumeintegral(x3D, gridmetrics; dim = 3)
    plt = plotmap!(ax, x2D, gridmetrics; colormap = :viridis)
    scatter!(ax, src_Ps[ksrc], marker=:star5, markersize=15, color=:red, strokecolor=:black, strokewidth=1)
    Colorbar(fig[1, 2], plt)

    # save plot
    str = keys(src_Ps)[ksrc]
    outputfile = joinpath(inputdir, "$(str)_injected_tracer_year$year.png")
    @info "Saving injection location as image file:\n  $(outputfile)"
    save(outputfile, fig)
end

x3D = fill(NaN, size(wet3D))
state = src_injections
year = 0
for ksrc in eachindex(keys(src_Ps))
    plotandprint!(x3D, state, ksrc, year)
end


function steponemonth!(u, p, m, y)
    prob = p.stepprobs[m]
    prob.b = u .+ p.δts[m] * src(y) # source only active first 10 years
    u .= solve!(prob).u
    return u
end
function steponeyear!(u, p, y)
    for m in eachindex(p.stepprobs)
        steponemonth!(u, p, m, y)
    end
    return u
end


@info "Time-stepping loop"
# Nyears = 500
Nyears = 11
src_mass = v' * sum(src(y) * sum(δts) for y in 1:Nyears)
u = zeros(N, Nsrc)
# Preallocate what I save? (may be worth it to save to disk instead, especially oif saving full field)
umass = Vector{Float64}(undef, Nyears, Nsrc)
# uprofile = Array{Float64}(undef, Nsrc, length(lev), Nyears)
@showprogress for y in 1:Nyears
    steponeyear!(u, p, y)
    # save mass time series
    umass[y, :] = v' * u
    # x3D[wet3D] .= u
    # uprofile[:, y] .= average(x3D, gridmetrics; dims = (1, 2))
    # if mod(y, 10) == 0
    #     text = "$(round(100umass[y] / src_mass, sigdigits = 2))%"
    #     plotandprint!(x3D, u, y, text)
    # end
end

outputfile = joinpath(inputdir, "timeseries_injected.jld2")
@info "Saving injection time series file:\n  $(outputfile)"
save(outputfile,
    Dict(
        "umass" => umass,
        # "uprofile" => uprofile,
        "src_mass" => src_mass,
        "src_Ps" => src_Ps,
    )
)

# Hovmoller plot of profiles (commented out because not saving profiles anymore)
# fig = Figure(size=(600, 300))
# ax = Axis(fig[1, 1], title="Injected tracer")
# co = contourf!(ax, 1:size(uprofile, 2), Array(lev), uprofile'; colormap = :viridis)
# ylims!(ax, (6000, 0))
# Colorbar(fig[1, 2], co)

# # save plot
# outputfile = joinpath(inputdir, "$(nearesttown)_injected_tracer_year_Hovmoller_profile.png")
# @info "Saving injection location as image file:\n  $(outputfile)"
# save(outputfile, fig)