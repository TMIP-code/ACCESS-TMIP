# qsub -I -P xv83 -q express -l mem=100GB -l storage=scratch/gh0+scratch/xv83 -l walltime=02:00:00 -l ncpus=24

using Pkg
Pkg.activate(".")
Pkg.instantiate()
const nprocs = 24

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
using StatsBase
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
using GeometryBasics
using LibGEOS
using GeometryOps

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
    # time_window = "Jan1990-Dec1999"
    experiment = "ssp370"
    time_window = "Jan2030-Dec2039"
    # time_window = "Jan2090-Dec2099"
    WRITEDATA = "true"
else
    experiment, member, time_window, WRITEDATA = ARGS
end
WRITEDATA = parse(Bool, WRITEDATA)
@show experiment
@show member
@show time_window

# preferred diffusivities
κVdeep = 3.0e-5 # m^2/s
κVML = 1.0      # m^2/s
κH = 300.0      # m^2/s
@show κVdeep
@show κVML
@show κH
κVdeep_str = "kVdeep" * format(κVdeep, conversion="e")
κVML_str = "kVML" * format(κVML, conversion="e")
κH_str = "kH" * format(κH, conversion="d")

upwind = false
@show upwind
upwind_str = upwind ? "" : "_centered"
upwind_str2 = upwind ? "upwind" : "centered"


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
    inputfile = joinpath(cycloinputdir, "cyclo_matrix$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)_$m.jld2")
    @info "Loading matrix from $inputfile"
    T = load(inputfile, "T")
    V⁻¹ * T' * V + Ω
end
@time "building mean_days_in_months" mean_days_in_months = map(months) do m
    inputfile = joinpath(cycloinputdir, "cyclo_matrix$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)_$m.jld2")
    load(inputfile, "mean_days_in_month")
end
# So the δt that multiplies M̃ₜ is δ(t..t+1)
# which is 0.5 of the mean days in months k and k+1
δts = map(months) do m
    ustrip(s, (mean_days_in_months[mod1(m + 1, 12)] + mean_days_in_months[m]) / 2 * d)
end

const Δt = sum(δts)
const M̃̄ = mean(M̃s, Weights(mean_days_in_months))

matrix_type = Pardiso.REAL_STRUCT_SYM
@show solver = MKLPardisoIterate(; nprocs, matrix_type)

prob = LinearProblem(I + Δt * M̃̄, ones(N))
yearprob = init(prob, solver, rtol = 1e-10)




function stepbackoneyear!(ℊ̃, ∫ℊ̃dt, yearprob, t)
    yearprob.b = ℊ̃ # ℊ̃ₘ₋₁ = Ãₘ₋₁⁻¹ ℊ̃ₘ
    ℊ̃ .= solve!(yearprob).u
    ∫ℊ̃dt .+= ℊ̃ * Δt
    # Write yearly values to 4D Array
    if WRITEDATA
        # data4D[:,:,:,t] .= OceanTransportMatrixBuilder.as3D(ℊ̃, wet3D)
        # Save yearly data at seafloor
        x3D[wet3D] .= ℊ̃
        ℊ̃_seafloor[:,:,t] .= seafloorvalue(x3D, wet3D)
        x3D[wet3D] .= 1 .- ∫ℊ̃dt # sequestration efficiency, ℰ
        ℰ_seafloor[:,:,t] .= seafloorvalue(x3D, wet3D)
    end
    return ℊ̃
end



Nyears = 3000

# Initial condition
ℊ̃ = Ω * ones(N)
∫ℊ̃dt = zeros(N)

# Preallocate what I save? (may be worth it to save to disk instead, especially oif saving full field)
nx, ny, nz = size(wet3D)
ℊ̃_seafloor = fill(NaN, nx, ny, Nyears)
ℰ_seafloor = fill(NaN, nx, ny, Nyears)
x3D = fill(NaN, nx, ny, nz)

@showprogress "Time-stepping loop" for t in 1:Nyears
    stepbackoneyear!(ℊ̃, ∫ℊ̃dt, yearprob, t)
end





# save yearly values of ℊ̃
axlist = (dims(volcello_ds["volcello"])[1:2]..., dims(DimArray(ones(Nyears), Ti(1:Nyears)))[1])
ℊ̃cube = DimensionalData.rebuild(volcello_ds["volcello"];
    data = ℊ̃_seafloor,
    dims = axlist,
    metadata = Dict(
        "description" => "yearly adjoint propagator",
        "model" => model,
        "experiment" => experiment,
        "member" => member,
        "time window" => time_window,
        "upwind" => upwind_str2,
        "units" => "s^-1",
        "Ti unit" => "yr",
    )
)
arrays = Dict(:calgtilde => ℊ̃cube, :lat => volcello_ds.lat, :lon => volcello_ds.lon)
ds = Dataset(; volcello_ds.properties, arrays...)
# Save to netCDF file

outputfile = joinpath(inputdir, "calgtilde$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)_yearly.nc")
@info "Saving adjoint propagrator as netCDF file:\n  $(outputfile)"
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)

# save yearly values of ℰ
axlist = (dims(volcello_ds["volcello"])[1:2]..., dims(DimArray(ones(Nyears), Ti(1:Nyears)))[1])
ℰcube = DimensionalData.rebuild(volcello_ds["volcello"];
    data = ℰ_seafloor,
    dims = axlist,
    metadata = Dict(
        "description" => "yearly sequestration efficiency",
        "model" => model,
        "experiment" => experiment,
        "member" => member,
        "time window" => time_window,
        "units" => "s^-1",
        "Ti unit" => "yr",
    )
)
arrays = Dict(:seqeff => ℰcube, :lat => volcello_ds.lat, :lon => volcello_ds.lon)
ds = Dataset(; volcello_ds.properties, arrays...)
# Save to netCDF file
outputfile = joinpath(inputdir, "seqeff$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)_yearly.nc")

@info "Saving adjoint propagrator as netCDF file:\n  $(outputfile)"
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)

