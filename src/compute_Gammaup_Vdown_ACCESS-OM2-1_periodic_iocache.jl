# qsub -I -P y99 -q normalbw -l mem=256GB -l storage=scratch/gh0+scratch/y99+gdata/xp65 -l walltime=02:00:00 -l ncpus=28
# qsub -I -P y99 -q hugemem -l mem=735GB -l storage=scratch/gh0+scratch/y99+gdata/xp65 -l walltime=02:00:00 -l ncpus=24

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
using ProgressMeter

model = "ACCESS-OM2-1"
experiment = "1deg_jra55_iaf_omip2_cycle6"
time_window = "Jan1960-Dec1979"
@show inputdir = "/scratch/y99/TMIP/data/$model/$experiment/$time_window"

# preferred diffusivities
κVdeep = 3.0e-5 # m^2/s
κVML = 1.0      # m^2/s
κH = 300.0      # m^2/s (grid-scaling by sqrt(area))
@show κVdeep
@show κVML
@show κH
κVdeep_str = "kVdeep" * format(κVdeep, conversion = "e")
κVML_str = "kVML" * format(κVML, conversion = "e")
κH_str = "kH" * format(κH, conversion = "d")

upwind = false
@show upwind
upwind_str = upwind ? "" : "_centered"
upwind_str2 = upwind ? "upwind" : "centered"

# Load areacello and volcello for grid geometry
areacello_ds = open_dataset(joinpath(inputdir, "area_t.nc"))
dht_ds = open_dataset(joinpath(inputdir, "dht.nc")) # <- (new) cell thickness?
lev = dht_ds.st_ocean


# Unfortunately ACCESS-OM2 raw data does not have coordinates of cell vertices
# So instead I go back to the source: the supergrids
include("supergrid.jl")
(; lon, lat, areacello, lon_vertices, lat_vertices) = supergrid(model; dims = dims(dht_ds.dht)[1:2])

# Make indices (from yearly volcello)
dht = readcubedata(dht_ds.dht)
volcello = readcubedata(dht .* areacello)
gridmetrics = makegridmetrics(;
    areacello, volcello, lon, lat, lev,
    lon_vertices, lat_vertices
)
# Make indices
indices = makeindices(gridmetrics.v3D)
(; wet3D) = indices

# surface/interior indices
issrf = let
    issrf3D = falses(size(wet3D))
    issrf3D[:, :, 1] .= true
    issrf3D[wet3D]
end
idx_surface = findall(issrf)
idx_interior = findall(.!issrf)
Nᵢ = length(idx_interior)
Nₛ = length(idx_surface)

months = 1:12

# δts between climatological months
# So the δt that multiplies M̃ₜ is δ(t..t+1)
# which is 0.5 of the mean days in months t and t+1
# Load the monthly dht dataset
dht_periodic_ds = open_dataset(joinpath(inputdir, "dht_periodic.nc")) # <- (new) cell thickness?
mean_days_in_month = dht_periodic_ds.mean_days_in_month |> Array
δts = map(months) do m
    ustrip(s, (mean_days_in_month[mod1(m + 1, 12)] + mean_days_in_month[m]) / 2 * d)
end

# TODO: Write matrices in separate scripts (to allow parallel computations)
# TODO: Also read/write the solver cache to save on expensive factorizations/preconditioners
# But for now, build all the matrices and store everything in memory like I have done before.
# Build matrices
@time "building the monthly TMs" T_periodic = map(months) do m
    inputfile = joinpath(inputdir, "monthly_matrix$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)_$(m).jld2")
    @info "Loading TM from $inputfile"
    load(inputfile, "T")
end
@time "building the monthly volume vectors" v_periodic = map(months) do m
    dht = readcubedata(dht_periodic_ds.dht[month = At(m)])
    volcello = readcubedata(dht .* areacello)
    volcello.data[wet3D]
end


# Solve once the steady-state problem to get initial guess
# TODO: Maybe remove this step if it's too slow at high res?
oceanadjoint(T, v) = sparse(Diagonal(1 ./ v)) * transpose(T) * sparse(Diagonal(v))
v_mean = mean(v_periodic)
V_mean = sparse(Diagonal(v_mean))
V⁻¹_mean = sparse(Diagonal(1 ./ v_mean))
T_mean = mean(T_periodic)
Tᵃ_mean = oceanadjoint(T_mean, v_mean)
Tᵃᵢᵢ_mean = Tᵃ_mean[idx_interior, idx_interior] #
Δt = sum(δts)
matrix_type = Pardiso.REAL_SYM
@show solver = MKLPardisoIterate(; nprocs, matrix_type)
@time "initial state solve" u0 = solve(LinearProblem(Tᵃᵢᵢ_mean, ones(Nᵢ)), solver, rtol = 1.0e-10, verbose = true).u
@show norm(Tᵃᵢᵢ_mean * u0 - ones(Nᵢ)) / norm(ones(Nᵢ))


############## START from periodic solver WIP #####################
# Left Preconditioner needs a new type
struct PeriodicPreconditioner
    prob
end
Base.eltype(::PeriodicPreconditioner) = Float64
function LinearAlgebra.ldiv!(Pl::PeriodicPreconditioner, x::AbstractVector)
    @info "applying Pl"
    Pl.prob.b = x
    solve!(Pl.prob)
    x .= Pl.prob.u .- x # Note the -x (following Bardin et al)
    return x
end
function LinearAlgebra.ldiv!(y::AbstractVector, Pl::PeriodicPreconditioner, x::AbstractVector)
    Pl.prob.b = x
    solve!(Pl.prob)
    y .= Pl.prob.u .- x # Note the -x (following Bardin et al)
    return y
end


function initstepprob(A)
    prob = LinearProblem(A, ones(size(A, 1)))
    return init(prob, solver, rtol = 1.0e-10)
end

p = []
# Instead of collecting all the linear problems into a vector
PBS_JOBID = ENV["PBS_JOBID"]
linearproblemfile(m) = joinpath(inputdir, "stepprob_$(PBS_JOBID)_month$(m).jld2")
for (m, δt, T, v) in zip(months, δts, T_periodic, v_periodic)
    @info "Initializing linear problem for month $m"
    prob = initstepprob(I + δt * oceanadjoint(T, v)[idx_interior, idx_interior])
    save(linearproblemfile(m), Dict("prob" => prob))
end

function stepbackonemonth!(du, u, p, m)
    prob = load(linearproblemfile(m), "prob")
    prob.b = u .+ δts[m] # xₘ₊₁ = Aₘ₊₁⁻¹ (xₘ + δt 1) # CHECK m index is not off by 1
    du .= solve!(prob).u
    save(linearproblemfile(m), Dict("prob" => prob))
    return du
end
function jvpstep!(dv, v, p, m)
    prob = load(linearproblemfile(m), "prob")
    prob.b = v # xₘ₊₁ = Aₘ₊₁⁻¹ (xₘ + δt 1) # CHECK m index is not off by 1
    dv .= solve!(prob).u
    save(linearproblemfile(m), Dict("prob" => prob))
    return dv
end
function stepbackoneyear!(du, u, p)
    du .= u
    for m in reverse(months)
        stepbackonemonth!(du, du, p, m)
    end
    return du
end
function jvponeyear!(dv, v, p)
    dv .= v
    for m in reverse(months)
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

Plprob = LinearProblem(-Δt * Tᵃᵢᵢ_mean, ones(Nᵢ))  # following Bardin et al. (M -> -M though)
Plprob = init(Plprob, solver, rtol = 1.0e-10)
Pl = PeriodicPreconditioner(Plprob)
Pr = I
precs = Returns((Pl, Pr))

@info "solve periodic state"
# @time sol = solve(nonlinearprob, NewtonRaphson(linsolve = KrylovJL_GMRES(precs = precs)), verbose = true, reltol=1e-10, abstol=Inf);
@time sol! = solve(nonlinearprob!, NewtonRaphson(linsolve = KrylovJL_GMRES(precs = precs, rtol = 1.0e-12)); show_trace = Val(true), reltol = Inf, abstol = 1.0e-10norm(u0, Inf));


@info "Check the RMS drift, should be order 10⁻⁹‰ (1e-9 per thousands)"
du = deepcopy(u0)
@show norm(G!(du, sol!.u, p), Inf) / norm(sol!.u, Inf) |> u"permille"

# Save periodic reemergence time
du = sol!.u # The last month solved for is January (m = 1, implicit in backward time)
Γꜛ4D = reduce(
    (a, b) -> cat(b, a, dims = 4), # <- note how the order is reversed here
    map(reverse(months)) do m
        stepbackonemonth!(du, du, p, m) # Starting from du = January
        Γꜛ3D = OceanTransportMatrixBuilder.as3D([zeros(Nₛ); du], wet3D)
        reshape(Γꜛ3D, (size(wet3D)..., 1))
    end
)
Γꜛyax = YAXArray(
    dims(dht_periodic_ds.dht),
    ustrip.(yr, Γꜛ4D * s),
    Dict(
        "description" => "periodic reemergence time (time until next surface contact)",
        "solver" => "MKLPardisoIterate",
        "model" => model,
        "experiment" => experiment,
        "time window" => time_window,
        "upwind" => upwind_str2,
        "units" => "yr",
    )
)

arrays = Dict(:Gammaup => Γꜛyax, :lat => lat, :lon => lon)
ds = Dataset(; properties = Dict(), arrays...)
# Save to netCDF file
outputfile = joinpath(inputdir, "periodic_Gup_$(κVdeep_str)_$(κH_str)_$(κVML_str)_MKLPardisoIterate_iocache.nc")
@info "Saving age as netCDF file:\n  $(outputfile)"
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)

# Following Holzer et al. (2020) or Pasquier et al. (2024) the volume 𝒱↓ is given by
#     𝒱ꜜ = −Aₛ⁻¹ Vₛ Tᵃₛᵢ Tᵃᵢᵢ⁻¹ 1ᵢ
# But this is the same as
#     𝒱ꜜ = −Aₛ⁻¹ Vₛ Tᵃₛᵢ Γꜛᵢ
# So I might as well compute 𝒱ꜜ now since I just computed Γꜛᵢ
# Unit is m⁻² m³ s⁻¹ s = interior volume (m³) / surface area (m²)
# Note: In Pasquier et al. (2024) I plot this as %(interior volume) / 10,000km²


wet2D = wet3D[:, :, 1]
isurface2D = findall(wet2D)
Aₛ⁻¹ = sparse(Diagonal(1 ./ areacello.data[isurface2D]))
function as2D(xₛ)
    x2D = fill(NaN, size(wet2D))
    x2D[isurface2D] .= xₛ
    return x2D
end


𝒱ꜜ3D = reduce(
    (a, b) -> cat(a, b, dims = 3), # no need to reverse order here
    map(months) do m
        Γꜛ3D = Γꜛ4D[:, :, :, m]
        Γꜛᵢ = Γꜛ3D[wet3D][idx_interior]
        T = T_periodic[m]
        v = v_periodic[m]
        Tᵃ = oceanadjoint(T, v)
        Tᵃₛᵢ = Tᵃ[idx_surface, idx_interior]
        Vₛ = sparse(Diagonal(v[idx_surface]))
        𝒱ꜜ = -Aₛ⁻¹ * Vₛ * Tᵃₛᵢ * Γꜛᵢ
        reshape(as2D(𝒱ꜜ), (size(wet2D)..., 1))
    end
)

# Save 𝒱↑ as netCDF file
𝒱ꜜyax = YAXArray(
    dims(dht_periodic_ds.dht)[[1, 2, 4]],
    𝒱ꜜ3D,
    Dict(
        "description" => "periodic ocean volume ventilated down by unit area",
        "model" => model,
        "experiment" => experiment,
        "time window" => time_window,
        "units" => "m^3/m^2",
        "solver" => "MKLPardisoIterate",
        "upwind" => upwind_str2,
    )
)
arrays = Dict(:Vdown => 𝒱ꜜyax, :lat => lat, :lon => lon)
ds = Dataset(; properties = Dict(), arrays...)
# Save to netCDF file
outputfile = joinpath(inputdir, "periodic_Vdown_$(κVdeep_str)_$(κH_str)_$(κVML_str)_MKLPardisoIterate_iocache.nc")
@info "Saving Vdown as netCDF file:\n  $(outputfile)"
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)

# # clean up
# for m in months
#     rm(linearproblemfile(m); force = true)
# end