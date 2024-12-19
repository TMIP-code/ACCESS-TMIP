# qsub -I -P xv83 -q hugemem -l mem=500GB -l storage=scratch/gh0+scratch/xv83 -l walltime=02:00:00 -l ncpus=48

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

# TODO: Update below for adjoint propagator ℊ̃ (\widetilde{\mathcal{g}})
#
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




# script options
@show model = "ACCESS-ESM1-5"
if isempty(ARGS)
    member = "r1i1p1f1"
    # experiment = "historical"
    # time_window = "Jan1850-Dec1859"
    srcname = "Karratha"
    # time_window = "Jan1990-Dec1999"
    experiment = "ssp370"
    time_window = "Jan2030-Dec2039"
    # time_window = "Jan2090-Dec2099"
    finalmonth = "1" # January
    WRITEDATA = "true"
else
    experiment, member, time_window, srcname, finalmonth, WRITEDATA = ARGS
end
finalmonth = parse(Int, finalmonth)
WRITEDATA = parse(Bool, WRITEDATA)
@show experiment
@show member
@show time_window
@show srcname
@show finalmonth




# Load areacello and volcello for grid geometry
fixedvarsinputdir = "/scratch/xv83/TMIP/data/$model"
volcello_ds = open_dataset(joinpath(fixedvarsinputdir, "volcello.nc"))
areacello_ds = open_dataset(joinpath(fixedvarsinputdir, "areacello.nc"))

# Gadi directory for input files
lumpby = "month"
inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)"
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

# Make indices
indices = makeindices(v3D)
(; N, wet3D) = indices

# Make V diagnoal matrix of volumes
V = sparse(Diagonal(v3D[wet3D]))
V⁻¹ = sparse(Diagonal(1 ./ v3D[wet3D]))

issrf = let
    issrf3D = falses(size(wet3D))
    issrf3D[:,:,1] .= true
    issrf3D[wet3D]
end
Ω = sparse(Diagonal(Float64.(issrf)))

months = 1:12

# Build matrices
@time "building M̃s" M̃s = map(months) do m
    inputfile = joinpath(cycloinputdir, "cyclo_matrix_$m.jld2")
    @info "Loading matrix from $inputfile"
    T = load(inputfile, "T")
    V⁻¹ * T' * V + Ω
end
@time "building mean_days_in_months" mean_days_in_months = map(months) do m
    inputfile = joinpath(cycloinputdir, "cyclo_matrix_$m.jld2")
    load(inputfile, "mean_days_in_month")
end
# So the δt that multiplies M̃ₜ is δ(t..t+1)
# which is 0.5 of the mean days in months k and k+1
δts = map(months) do m
    ustrip(s, (mean_days_in_months[mod1(m + 1, 12)] + mean_days_in_months[m]) / 2 * d)
end

function initstepprob(A)
    prob = LinearProblem(A, ones(N))
    return init(prob, MKLPardisoIterate(; nprocs = 48), rtol = 1e-10)
end

p = (;
    δts,
    monthprob = [initstepprob(I + δt * M̃) for (δt, M̃) in zip(δts, M̃s)]
)


# TODO (maybe WONTFIX) update plot below?
# # Plot
# function plotandprint!(x3D, state, year, text="")
#     x3D[wet3D] = state
#     fig = Figure(size=(600, 300))
#     ax = Axis(fig[1, 1], title="Injected tracer (year $year: $text sequestered)")
#     x2D = volumeintegral(x3D, gridmetrics; dim = 3)
#     plt = plotmap!(ax, x2D, gridmetrics; colormap = :viridis)
#     # scatter!(ax, src_P, marker=:star5, markersize=15, color=:red, strokecolor=:black, strokewidth=1)
#     Colorbar(fig[1, 2], plt)

#     # save plot
#     outputfile = joinpath(inputdir, "$(srcname)_injected_tracer_year$year.png")
#     @info "Saving injection location as image file:\n  $(outputfile)"
#     save(outputfile, fig)
# end

# x3D = fill(NaN, size(wet3D))
# year = 0
# # plotandprint!(x3D, src1D, year)

# Get index for M̃ given number of months elapsed t
M̃idx(t) = mod1(finalmonth - t, 12)


function stepbackonemonth!(ℊ̃, p, t)
    m = M̃idx(t) # Get the month (index used for M̃)
    prob = p.monthprob[m]
    prob.b = ℊ̃ # ℊ̃ₘ₋₁ = Ãₘ₋₁⁻¹ ℊ̃ₘ
    ℊ̃ .= solve!(prob).u
    # Write monthly values to 4D Array
    if WRITEDATA
        data4D[:,:,:,t] .= OceanTransportMatrixBuilder.as3D(ℊ̃, wet3D)
    end
    return ℊ̃
end



# Nyears = 1001
Nyears = 22

Nmonths = 12Nyears

# Initial condition
ℊ̃ = V * Ω

# Preallocate what I save? (may be worth it to save to disk instead, especially oif saving full field)
data4D = Array{Float64}(undef, size(wet3D)..., Nmonths)

@showprogress "Time-stepping loop" for t in 1:Nmonths
    steponemonth!(ℊ̃, p, t)
end




# save monthly values

axlist = (dims(volcello_ds["volcello"])..., dims(DimArray(ones(Nmonths), Ti(1:Nmonths)))[1])

cube4D = rebuild(volcello_ds["volcello"];
    data = data4D,
    dims = axlist,
    metadata = Dict(
        "origin" => "cyclo-stationary calg tilde",
        "finalmonth" => finalmonth,
        "model" => model,
        "experiment" => experiment,
        "member" => member,
        "time window" => time_window,
        "units" => "yr",
    )
)

arrays = Dict(:calgtilde => cube4D, :lat => volcello_ds.lat, :lon => volcello_ds.lon)
ds = Dataset(; volcello_ds.properties, arrays...)

# Save Γinyr3D to netCDF file
finalmonthstr = format(finalmonth, width = 2, zeropadding = true)
outputfile = joinpath(inputdir, "calgtilde_$(finalmonthstr).nc")
@info "Saving reemergence time as netCDF file:\n  $(outputfile)"
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)
