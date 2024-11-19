# qsub -I -P xv83 -q hugemem -l mem=720GB -l storage=gdata/gh0+scratch/xv83 -l walltime=02:00:00 -l ncpus=48

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





# Load matrix and grid metrics
model = "ACCESS-ESM1-5"
member = "r1i1p1f1"
experiment = "historical"
time_window = "Jan1990-Dec1999"
lumpby = "season"
steps = ("DJF", "MAM", "JJA", "SON")
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
(; lon_vertices, lat_vertices, v3D) = gridmetrics

# Make indices
indices = makeindices(v3D)
(; N, wet3D) = indices
# The ideal age equation is
#   ∂ₜx + T x = 1 - Ω x
# Applying Backward Euler time step gives


issrf = let
    issrf3D = falses(size(wet3D))
    issrf3D[:,:,1] .= true
    issrf3D[wet3D]
end
Ω = sparse(Diagonal(Float64.(issrf)))


# Build matrices
@time "building Ms" Ms = [
    begin
        inputfile = joinpath(cycloinputdir, "cyclo_matrix_$step.jld2")
        @info "Loading matrices + metrics as $inputfile"
        T = load(inputfile)["T"]
        T + Ω
    end
    for step in steps
]


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
Δt = sum(δt for _ in steps)

Plprob = LinearProblem(-Δt * M̄, ones(N))  # following Bardin et al. (M -> -M though)
Plprob = init(Plprob, MKLPardisoIterate(; nprocs = 48), rtol = 1e-10)
Pl = CycloPreconditioner(Plprob)
Pr = I
precs = Returns((Pl, Pr))

@time "initial state solve" u0 = solve(LinearProblem(M̄, ones(N)), MKLPardisoIterate(; nprocs = 48), rtol = 1e-10).u
@show norm(M̄ * u0 - ones(N)) / norm(ones(N))

function initstepprob(A)
    prob = LinearProblem(A, δt * ones(N))
    return init(prob, MKLPardisoIterate(; nprocs = 48), rtol = 1e-10)
end

p = (;
    δt,
    stepprob = [initstepprob(I + δt * M) for M in Ms]
)
function mystep!(du, u, p, m)
    prob = p.stepprob[m]
    prob.b = u .+ p.δt # xₘ₊₁ = Aₘ₊₁⁻¹ (xₘ + δt 1) # CHECK m index is not off by 1
    du .= solve!(prob).u
    return du
end
function jvpstep!(dv, v, p, m)
    prob = p.stepprob[m]
    prob.b = v # xₘ₊₁ = Aₘ₊₁⁻¹ (xₘ + δt 1) # CHECK m index is not off by 1
    dv .= solve!(prob).u
    return dv
end
function steponeyear!(du, u, p)
    du .= u
    for m in eachindex(p.stepprob)
        mystep!(du, du, p, m)
    end
    return du
end
function jvponeyear!(dv, v, p)
    dv .= v
    for m in eachindex(p.stepprob)
        jvpstep!(dv, dv, p, m)
    end
    return dv
end
function G!(du, u, p)
    steponeyear!(du, u, p)
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

@info "solve seasonal steady state"
# @time sol = solve(nonlinearprob, NewtonRaphson(linsolve = KrylovJL_GMRES(precs = precs)), verbose = true, reltol=1e-10, abstol=Inf);
@time sol! = solve(nonlinearprob!, NewtonRaphson(linsolve = KrylovJL_GMRES(precs = precs, rtol=1e-12)); show_trace = Val(true), reltol=Inf, abstol=1e-10norm(u0, Inf));


@info "Check the RMS drift, should be order 10⁻¹¹‰ (1e-11 per thousands)"
du = deepcopy(u0)
@show norm(G!(du, sol!.u, p), Inf) / norm(sol!.u, Inf) |> u"permille"
