
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
using CUDSS
using CUDA
using LinearOperators
CUDA.set_runtime_version!(v"12.9.1")
@show CUDA.versioninfo()
using CUDA.CUSPARSE, CUDA.CUSOLVER
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
(; wet3D, N) = indices

# surface/interior indices
issrf = let
    issrf3D = falses(size(wet3D))
    issrf3D[:, :, 1] .= true
    issrf3D[wet3D]
end
Nₛ = findfirst(.!issrf) - 1
Nᵢ = N - Nₛ
idx_surface = 1:Nₛ
idx_interior = (Nₛ + 1):N

# Figure out the k-first ordering
# By that I mean loop over z (depth) first
# The goal is for vertical diffusion
# to be represented by a tridiagonal matrix
# Note that here I only index the interior points!
L = findall(wet3D[:, :, 2:end])
wet3D_kfirst = permutedims(wet3D[:, :, 2:end], (3, 1, 2))
L_kfirst = findall(wet3D_kfirst)
# assign indices to L3D and get k-first ordering
L3D = zeros(Int, size(wet3D[:, :, 2:end]))
L3D[wet3D[:, :, 2:end]] .= 1:length(L)
idx = permutedims(L3D, (3, 1, 2))[L_kfirst]
# assign k-first indices to L3D_kfirst and get reverse ordering
L3D_kfirst = zeros(Int, size(wet3D_kfirst))
L3D_kfirst[L_kfirst] .= 1:length(L_kfirst)
idx_kfirst = permutedims(L3D_kfirst, (2, 3, 1))[L]
# Figure out the block indices ()
# for vertical diffusion tridiagonal
# The size of each block is the same as
# the number of vertical levels with in each water column, which is
K2D = sum(wet3D(:, :, 2:end), dims = 3)
# cumulative sum over wet points should give last index of each block
blk_last_idx = cumsum(K2D[K2D .≠ 0])
blk_first_idx = [1; blk_last_idx[1:end-1] .+ 1]

foo

@assert Nₛ == length(idx_surface)
@assert Nᵢ == length(idx_interior)
@assert collect(idx_surface) == findall(issrf)
@assert collect(idx_interior) == findall(.!issrf)

months = 1:12

# δts between climatological months
# So the δt that multiplies M̃ₜ is δ(t..t+1)
# which is 0.5 of the mean days in months t and t+1
# Load the monthly dht dataset
dht_periodic_ds = open_dataset(joinpath(inputdir, "dht_periodic.nc")) # <- (new) cell thickness?
mean_days_in_month = dht_periodic_ds.mean_days_in_month |> Array


# Follow Bardin et al. and do
# implicit time stepping for vertical diffusion
# and explicit time stepping for everything else
# so with T = A + D
#
# Age equation is
# Explicit equation is
#     xₜ₊₁ - xₜ + δt T xₜ = δt
# so rearranged is
#     xₜ₊₁ = xₜ - δt T xₜ + δt
#
# 3rd order Adams-Bashforth is
#     xₜ₊₁ - xₜ + δt/12 T [23 xₜ - 16 xₜ₋₁ + 5 xₜ₋₂] = δt
# so rearranged is
#     xₜ₊₁ = xₜ - δt/12 T [23 xₜ - 16 xₜ₋₁ + 5 xₜ₋₂] + δt
#
# Implicit is
#     xₜ₊₁ - xₜ + δt T xₜ₊₁ = δt
# so rearranged is
#     (I + δt T) xₜ₊₁ = xₜ + δt
#
# Implicit-explicit is
#     xₜ₊₁ - xₜ + δt D xₜ₊₁ + δt/12 (A + H) (23 xₜ - 16 xₜ₋₁ + 5 xₜ₋₂) = δt
# so rearranged is
#     D⁺ xₜ₊₁ = xₜ - AH (23 xₜ - 16 xₜ₋₁ + 5 xₜ₋₂) + δt
# with D⁺ = I + δt D and AH = δt/12 (A + H)

# "The tests used baroclinic time steps of 5400, 1800, and 400 s at 1, 0.25, and 0.1∘ (respectively)."
# Kiss et al. (2020)
@show δt = (model == "ACCESS-OM2-1") ? 5400.0 :
           (model == "ACCESS-OM2-025") ? 1800.0 :
           (model == "ACCESS-OM2-01") ? 400.0 :
           error("What's the time step?") # s


# TODO: Write matrices in separate scripts (to allow parallel computations)
# TODO: Also read/write the solver cache to save on expensive factorizations/preconditioners
# But for now, build all the matrices and store everything in memory like I have done before.
# Build matrices
@time "building the monthly TMs" TMs = map(months) do m
    inputfile = joinpath(inputdir, "OM2-1_TM_split_$(upwind_str)_$(κVdeep_str)_$(κH_str)_$(κVML_str)_$(m).jld2")
    @info "Loading TM from $inputfile"
    A = load(inputfile, "A")[idx_interior[idx], idx_interior[idx]]
    H = load(inputfile, "H")[idx_interior[idx], idx_interior[idx]]
    D = load(inputfile, "D")[idx_interior[idx], idx_interior[idx]]
    AH = δt / 12 * (A + H)
    D⁺ = I + δt * D
    # D⁺ = Tridiagonal(I + δt * D)
    # @assert Tridiagonal(I + δt * D) == I + δt * D
    (; AH, D⁺)
end
# @time "building the monthly volume vectors" v_periodic = map(months) do m
#     dht = readcubedata(dht_periodic_ds.dht[month = At(m)])
#     volcello = readcubedata(dht .* areacello)
#     volcello.data[wet3D]
# end



# 4 vectors on the GPU for 3rd-order Adams-Bashforth
xₜ₊₁ = CuVector(ones(Nᵢ))
xₜ = CuVector(ones(Nᵢ))
xₜ₋₁ = CuVector(ones(Nᵢ))
xₜ₋₂ = CuVector(ones(Nᵢ))

# 2 matrices on the GPU for D and A+H
AH = CuSparseMatrixCSR(TMs[1].AH)
D⁺ = CuSparseMatrixCSR(TMs[1].D⁺)
D⁺fac = lu(D⁺)




function timestep!(xₜ₊₁, xₜ, xₜ₋₁, xₜ₋₂, D⁺fac, AH, δt)
    # First do xₜ₊₁ <- AH (23 xₜ - 16 xₜ₋₁ + 5 xₜ₋₂)
    @time "mul!" mul!(xₜ₊₁, AH, @. 23.0 * xₜ - 16.0 * xₜ₋₁ + 5.0 * xₜ₋₂)
    # Then xₜ <- xₜ + xₜ₊₁ + δt
    # which is   xₜ + AH (23 xₜ - 16 xₜ₋₁ + 5 xₜ₋₂) + δt
    @time "add dt" xₜ .+= xₜ₊₁ .+ δt
    # Finally solve D⁺ xₜ₊₁ = ...
    @time "ldiv!" ldiv!(xₜ₊₁, D⁺fac, xₜ)
    # Update previous time steps
    @time "update t-2" xₜ₋₂ .= xₜ₋₁
    @time "update t-1" xₜ₋₁ .= xₜ
    @time "update t" xₜ .= xₜ₊₁
    return nothing
end

@time "warming up the GPU" timestep!(xₜ₊₁, xₜ, xₜ₋₁, xₜ₋₂, D⁺fac, AH, δt)
@time "time-step" timestep!(xₜ₊₁, xₜ, xₜ₋₁, xₜ₋₂, D⁺fac, AH, δt)

foo


function stepbackonemonth!(du, u, p, m)
    prob = stepprob[m]
    prob.b = u .+ δts[m] # xₘ₊₁ = Aₘ₊₁⁻¹ (xₘ + δt 1) # CHECK m index is not off by 1
    du .= solve!(prob).u
    return du
end
function jvpstep!(dv, v, p, m)
    prob = stepprob[m]
    prob.b = v # xₘ₊₁ = Aₘ₊₁⁻¹ (xₘ + δt 1) # CHECK m index is not off by 1
    dv .= solve!(prob).u
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
outputfile = joinpath(inputdir, "periodic_Gup_$(κVdeep_str)_$(κH_str)_$(κVML_str)_MKLPardisoIterate.nc")
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
outputfile = joinpath(inputdir, "periodic_Vdown_$(κVdeep_str)_$(κH_str)_$(κVML_str)_MKLPardisoIterate.nc")
@info "Saving Vdown as netCDF file:\n  $(outputfile)"
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)
