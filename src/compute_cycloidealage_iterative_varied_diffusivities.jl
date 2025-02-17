# qsub -I -P xv83 -q hugemem -l mem=360GB -l storage=scratch/gh0+scratch/xv83 -l walltime=02:00:00 -l ncpus=48

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
using Format
using Dates
using FileIO
using LinearSolve
import Pardiso # import Pardiso instead of using (to avoid name clash?)
using NonlinearSolve

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
# More succintly, if m represents the months (1..12)
#
#   (I + δt Mₘ₊₁) xₘ₊₁ = xₘ + δt sₘ₊₁          (1)
#
# Here δt = δt(m-1..m) is the time that separates
# the "center time" of climatological months m-1 and m.
# So the δt that multiplies Mₘ is δ(m..m+1).



# script options
@show model = "ACCESS-ESM1-5"
if isempty(ARGS)
    member = "r20i1p1f1"
    experiment = "historical"
    time_window = "Jan1850-Dec1859"
    # time_window = "Jan1990-Dec1999"
    # experiment = "ssp370"
    # time_window = "Jan2030-Dec2039"
    # time_window = "Jan2090-Dec2099"
    κVdeep = 1e-5 # m^2/s
    κVML = 0.1 # m^2/s
    κH = 500 # m^2/s
else
    # experiment, member, time_window, κVdeep_str_in, κH_str_in = ARGS
    # κVdeep = parse(Float64, κVdeep_str_in)
    # κH = parse(Float64, κH_str_in)
    experiment, member, time_window, κVdeep_str_in, κVML_str_in, κH_str_in = ARGS
    κVdeep = parse(Float64, κVdeep_str_in)
    κVML = parse(Float64, κVML_str_in)
    κH = parse(Float64, κH_str_in)
end
@show experiment
@show member
@show time_window
@show κVdeep
@show κVML
@show κH

κVdeep_str = "kVdeep" * format(κVdeep, conversion="e")
κVML_str = "kVML" * format(κVML, conversion="e")
κH_str = "kH" * format(κH, conversion="d")

lumpby = "month"
months = 1:12
Nmonths = length(months)

# Load areacello and volcello for grid geometry
fixedvarsinputdir = "/scratch/xv83/TMIP/data/$model"
volcello_ds = open_dataset(joinpath(fixedvarsinputdir, "volcello.nc"))
areacello_ds = open_dataset(joinpath(fixedvarsinputdir, "areacello.nc"))

# Gadi directory for input files
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


issrf = let
    issrf3D = falses(size(wet3D))
    issrf3D[:,:,1] .= true
    issrf3D[wet3D]
end
Ω = sparse(Diagonal(Float64.(issrf)))


@time "building mean_days_in_months" mean_days_in_months = map(months) do m
    inputfile = joinpath(cycloinputdir, "cyclo_matrix_$(κVML_str)_$m.jld2")
    load(inputfile, "mean_days_in_month")
end
# So the δt that multiplies Mₜ is δ(t-1..t)
# which is 0.5 of the mean days in months t-1 and t
δts = map(eachindex(months)) do m
    ustrip(s, (mean_days_in_months[mod1(m - 1, 12)] + mean_days_in_months[m]) / 2 * d)
end




# Build matrices
@time "building Ms" Ms = map(months) do m
    inputfile = joinpath(cycloinputdir, "cyclo_matrix_$(κVML_str)_$m.jld2")
    @info "Loading matrix from $inputfile"
    T = load(inputfile, "T")
    T + Ω
end

# Preconditioner

# Left Preconditioner needs a new type
struct CycloPreconditioner
    prob
end
Base.eltype(::CycloPreconditioner) = Float64
function LinearAlgebra.ldiv!(Pl::CycloPreconditioner, x::AbstractVector)
    @info "applying Pl"
    Pl.prob.b = x
    solve!(Pl.prob)
    x .= Pl.prob.u .- x # Note the -x (following Bardin et al)
    return x
end
function LinearAlgebra.ldiv!(y::AbstractVector, Pl::CycloPreconditioner, x::AbstractVector)
    Pl.prob.b = x
    solve!(Pl.prob)
    y .= Pl.prob.u .- x # Note the -x (following Bardin et al)
    return y
end
M̄ = mean(Ms) #
Δt = sum(δts)

Plprob = LinearProblem(-Δt * M̄, ones(N))  # following Bardin et al. (M -> -M though)
Plprob = init(Plprob, MKLPardisoIterate(; nprocs), rtol = 1e-10)
Pl = CycloPreconditioner(Plprob)
Pr = I
precs = Returns((Pl, Pr))

@time "initial state solve" u0 = solve(LinearProblem(M̄, ones(N)), MKLPardisoIterate(; nprocs), rtol = 1e-10).u
@show norm(M̄ * u0 - ones(N)) / norm(ones(N))

function initstepprob(A)
    prob = LinearProblem(A, ones(N))
    return init(prob, MKLPardisoIterate(; nprocs), rtol = 1e-10)
end

p = (;
    months,
    δts,
    stepprob = [initstepprob(I + δt * M) for (δt, M) in zip(δts, Ms)]
)
function stepforwardonemonth!(du, u, p, m)
    prob = p.stepprob[m]
    prob.b = u .+ p.δts[m] # xₘ₊₁ = Aₘ₊₁⁻¹ (xₘ + δt 1) # CHECK m index is not off by 1
    du .= solve!(prob).u
    return du
end
function jvpstep!(dv, v, p, m)
    prob = p.stepprob[m]
    prob.b = v # xₘ₊₁ = Aₘ₊₁⁻¹ (xₘ + δt 1) # CHECK m index is not off by 1
    dv .= solve!(prob).u
    return dv
end
function stepforwardoneyear!(du, u, p)
    du .= u
    for m in p.months
        stepforwardonemonth!(du, du, p, m)
    end
    return du
end
function jvponeyear!(dv, v, p)
    dv .= v
    for m in p.months
        jvpstep!(dv, dv, p, m)
    end
    return dv
end
function G!(du, u, p)
    stepforwardoneyear!(du, u, p)
    du .-= u
    return du
end
function jvp!(dv, v, u, p)
    jvponeyear!(dv, v, p)
    dv .-= v
    return dv
end
f! = NonlinearFunction(G!; jvp = jvp!)
nonlinearprob! = NonlinearProblem(f!, u0, p)

@info "solve cyclo-stationary state"
# @time sol = solve(nonlinearprob, NewtonRaphson(linsolve = KrylovJL_GMRES(precs = precs)), verbose = true, reltol=1e-10, abstol=Inf);
@time sol! = solve(nonlinearprob!, NewtonRaphson(linsolve = KrylovJL_GMRES(precs = precs, rtol=1e-12)); show_trace = Val(true), reltol=Inf, abstol=1e-10norm(u0, Inf));


@info "Check the RMS drift, should be order 10⁻¹¹‰ (1e-11 per thousands)"
du = deepcopy(u0)
@show norm(G!(du, sol!.u, p), Inf) / norm(sol!.u, Inf) |> u"permille"


# Save cyclo-stationary age
du = sol!.u # The last month solved for is December (m = 12, implicit in forward time)
data4D = reduce((a, b) -> cat(a, b, dims=4),
    map(months) do m
        stepforwardonemonth!(du, du, p, m)
        Γinyr = ustrip.(yr, du .* s)
        Γinyr3D = OceanTransportMatrixBuilder.as3D(Γinyr, wet3D)
        reshape(Γinyr3D, (size(wet3D)..., 1))
    end
)
axlist = (dims(volcello_ds["volcello"])..., dims(DimArray(ones(Nmonths), Ti(months)))[1])
cube4D = rebuild(volcello_ds["volcello"];
    data = data4D,
    dims = axlist,
    metadata = Dict(
        "origin" => "cyclo-stationary ideal age (by $lumpby) computed from $model $experiment $member $(time_window)",
        "units" => "yr",
    )
)

arrays = Dict(:age => cube4D, :lat => volcello_ds.lat, :lon => volcello_ds.lon)
ds = Dataset(; volcello_ds.properties, arrays...)

# Save Γinyr3D to netCDF file
outputfile = joinpath(cycloinputdir, "ideal_mean_age_$(κVdeep_str)_$(κVML_str)_$(κH_str).nc")
@info "Saving ideal mean age as netCDF file:\n  $(outputfile)"
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)

