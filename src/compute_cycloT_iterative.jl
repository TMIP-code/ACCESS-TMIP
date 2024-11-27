# qsub -I -P xv83 -q hugemem -l mem=360GB -l storage=gdata/gh0+scratch/xv83 -l walltime=02:00:00 -l ncpus=48

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
# lumpby = "season"
lumpby = "month"
# steps = ["DJF", "MAM", "JJA", "SON"]
steps = 1:12
Nsteps = length(steps)
δt = ustrip(s, 1yr / Nsteps) # TODO maybe use exact mean number of days (more important for monthly because Feb)?


# Gadi directory for input files
inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)"
cycloinputdir = joinpath(inputdir, "cyclo$lumpby")
# Load areacello and volcello for grid geometry
volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))
thetao_ds = open_dataset(joinpath(inputdir, "thetao.nc"))

# Load fixed variables in memory
areacello = readcubedata(areacello_ds.areacello)
volcello = readcubedata(volcello_ds.volcello)
thetao = readcubedata(thetao_ds.thetao)
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
# The equation for conservative θ or S is
#   ∂ₜx + T x = Ω (x_srf - x)
# Applying Backward Euler time step gives
#   (I + Δt M) xₖ₊₁ = xₖ + Δt Ω x_srf
# where M = T + Ω


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

src = Ω * thetao[wet3D]
@time "initial state solve" u0 = solve(LinearProblem(M̄, src), MKLPardisoIterate(; nprocs = 48), rtol = 1e-10).u
@show norm(M̄ * u0 - src) / norm(src)

function initstepprob(A)
    prob = LinearProblem(A, δt * src)
    return init(prob, MKLPardisoIterate(; nprocs = 48), rtol = 1e-10)
end


function mystep!(du, u, p, m)
    prob = p.stepprob[m]
    prob.b = u .+ p.δt * p.src # xₘ₊₁ = Aₘ₊₁⁻¹ (xₘ + δt Ω x_srf)
    du .= solve!(prob).u
    return du
end
function jvpstep!(dv, v, p, m)
    prob = p.stepprob[m]
    prob.b = v # xₘ₊₁ = Aₘ₊₁⁻¹ (xₘ + δt 1)
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

p = (;
    δt,
    stepprob = [initstepprob(I + δt * M) for M in Ms],
    src
)
nonlinearprob! = NonlinearProblem(f!, u0, p)

@info "solve cyclo-stationary state"
# @time sol = solve(nonlinearprob, NewtonRaphson(linsolve = KrylovJL_GMRES(precs = precs)), verbose = true, reltol=1e-10, abstol=Inf);
@time sol! = solve(nonlinearprob!, NewtonRaphson(linsolve = KrylovJL_GMRES(precs = precs, rtol=1e-12)); show_trace = Val(true), reltol=Inf, abstol=1e-10norm(u0, Inf));


@info "Check the RMS drift, should be order 10⁻¹¹‰ (1e-11 per thousands)"
du = deepcopy(u0)
@show norm(G!(du, sol!.u, p), Inf) / norm(sol!.u, Inf) |> u"permille"


# Save mean salinity
du = sol!.u
cube4D = reduce((a, b) -> cat(a, b, dims=Ti),
    (
        begin
            (m > 1) && mystep!(du, du, p, m)
            temperature3D = OceanTransportMatrixBuilder.as3D(du, wet3D)
            temperature4D = reshape(temperature3D, (size(wet3D)..., 1))
            axlist = (dims(volcello_ds["volcello"])..., dims(DimArray(ones(Nsteps), Ti(steps)))[1][m:m])
            temperature_YAXArray = rebuild(volcello_ds["volcello"];
                data = temperature4D,
                dims = axlist,
                metadata = Dict(
                    "origin" => "cyclo-stationary ideal age (by $lumpby) computed from $model $experiment $member $(time_window)",
                )
            )
        end
        for m in eachindex(steps)
    )
)

temperaturemean3D = mean(cube4D, dims=Ti)

arrays = Dict(:thetao => temperaturemean3D, :lat => volcello_ds.lat, :lon => volcello_ds.lon)
ds = Dataset(; volcello_ds.properties, arrays...)

# Save temperature3D to netCDF file
outputfile = joinpath(cycloinputdir, "mean_cyclo_thetao.nc")
@info "Saving ideal mean age as netCDF file:\n  $(outputfile)"
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)
