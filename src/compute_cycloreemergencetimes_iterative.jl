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

# The adjoint tracer equation is
#
#   -∂x̃(t)/∂t + T̃(t) x̃(t) = s(t) - Ω x̃(t)
#
# where Ω "relaxes" x̃ to zero in the top layer.
#
# Applying Forward Euler time step gives
# (which is the way to go since in adjoint equation time goes backwards)
#
#   (I + δt M̃(t-δt)) x̃(t-δt) = x̃(t) + δt s(t-δt)
#
# where M̃(t) = T̃(t) + Ω.
#
# More succintly, if m represents the months (1..12)
#
#   (I + δt Mₘ₋₁) xₘ₋₁ = xₘ + δt sₘ₋₁          (1)
#
# Here δt = δt(m-1..m) is the time that separates
# the "center time" of climatological months m-1 and m.
# So the δt that multiplies Mₘ is δ(m..m+1).



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
else
    experiment, member, time_window = ARGS
end
@show experiment
@show member
@show time_window

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

# Make V diagnoal matrix of volumes
V = sparse(Diagonal(v3D[wet3D]))
V⁻¹ = sparse(Diagonal(1 ./ v3D[wet3D]))


issrf = let
    issrf3D = falses(size(wet3D))
    issrf3D[:,:,1] .= true
    issrf3D[wet3D]
end
Ω = sparse(Diagonal(Float64.(issrf)))


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
# So the δt that multiplies Mₘ is δ(k..k+1)
# which is 0.5 of the mean days in months k and k+1
δts = map(eachindex(months)) do m
    ustrip(s, (mean_days_in_months[mod1(m + 1, 12)] + mean_days_in_months[m]) / 2 * d)
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
M̃̄ = mean(M̃s) #
Δt = sum(δts)

Plprob = LinearProblem(-Δt * M̃̄, ones(N))  # following Bardin et al. (M -> -M though)
Plprob = init(Plprob, MKLPardisoIterate(; nprocs = 48), rtol = 1e-10)
Pl = CycloPreconditioner(Plprob)
Pr = I
precs = Returns((Pl, Pr))

@time "initial state solve" u0 = solve(LinearProblem(M̃̄, ones(N)), MKLPardisoIterate(; nprocs = 48), rtol = 1e-10).u
@show norm(M̃̄ * u0 - ones(N)) / norm(ones(N))

function initstepprob(A)
    prob = LinearProblem(A, ones(N))
    return init(prob, MKLPardisoIterate(; nprocs = 48), rtol = 1e-10)
end

p = (;
    months,
    δts,
    stepprob = [initstepprob(I + δt * M̃) for (δt, M̃) in zip(δts, M̃s)]
)
function stepbackonemonth!(du, u, p, m)
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
function stepbackoneyear!(du, u, p)
    du .= u
    for m in reverse(p.months)
        stepbackonemonth!(du, du, p, m)
    end
    return du
end
function jvponeyear!(dv, v, p)
    dv .= v
    for m in reverse(p.months)
        jvpstep!(dv, dv, p, m)
    end
    return dv
end
function G!(du, u, p)
    stepbackoneyear!(du, u, p)
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


# Save cyclo-stationary reemergence time
du = sol!.u # The last month solved for is January (m = 1, implicit in backward time)
data4D = reduce((a, b) -> cat(b, a, dims=4), # <- note how the order is reversed here
    map(reverse(months)) do m
        stepbackonemonth!(du, du, p, m) # Starting from du = January
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
        "origin" => "cyclo-stationary reemergence time (by $lumpby) computed from $model $experiment $member $(time_window)",
        "units" => "yr",
    )
)

arrays = Dict(:adjointage => cube4D, :lat => volcello_ds.lat, :lon => volcello_ds.lon)
ds = Dataset(; volcello_ds.properties, arrays...)

# Save Γinyr3D to netCDF file
outputfile = joinpath(cycloinputdir, "reemergence_time.nc")
@info "Saving reemergence time as netCDF file:\n  $(outputfile)"
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)
