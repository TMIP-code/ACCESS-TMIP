# qsub -I -P y99 -q normal -l mem=190GB -l storage=scratch/gh0+scratch/y99+gdata/xp65 -l walltime=02:00:00 -l ncpus=48

using Pkg
Pkg.activate(".")
Pkg.instantiate()
const nprocs = 48

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
# TODO: caputre correlations between transport and dht
# z* coordinate varies with time in ACCESS-OM2
# volcello_ds = open_dataset(joinpath(fixedvarsinputdir, "volcello.nc")) # <- not in ACCESS-OM2; must be built from dht * area


# Load fixed variables in memory
areacello_OM2 = replace(readcubedata(areacello_ds.area_t), missing => NaN) # This is required for later multiplication
dht = readcubedata(dht_ds.dht)
lon_OM2 = readcubedata(areacello_ds.geolon_t)
lat_OM2 = readcubedata(areacello_ds.geolat_t)
lev = dht_ds.st_ocean

# Unfortunately ACCESS-OM2 raw data does not have coordinates of cell vertices
# So instead I go back to the source: the supergrids
include("supergrid.jl")
(; lon, lat, areacello, lon_vertices, lat_vertices) = supergrid(model; dims = dims(dht_ds.dht)[1:2])

@show size(lon_vertices)

volcello = readcubedata(dht .* areacello)

# Make makegridmetrics
gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
(; lon_vertices, lat_vertices, v3D) = gridmetrics

# Make indices
indices = makeindices(v3D)
(; N, wet3D) = indices

# Make V diagnoal matrix of volumes
V = sparse(Diagonal(v3D[wet3D]))

issrf = let
    issrf3D = falses(size(wet3D))
    issrf3D[:, :, 1] .= true
    issrf3D[wet3D]
end
Ω = sparse(Diagonal(Float64.(issrf)))

TMfile = joinpath(inputdir, "yearly_matrix_$(κVdeep_str)_$(κH_str)_$(κVML_str).jld2")
@info "Loading matrix from $TMfile"
T = load(TMfile, "T")
@info "Matrix size: $(size(T)), nnz = $(nnz(T))"

M = T + Ω

matrix_type = Pardiso.REAL_SYM
@show solver = MKLPardisoFactorize(; nprocs, matrix_type)

@info "Solve full problem but with GMRES + LU-by-layer preconditioner"
maxindexinlayer = [sum(wet3D[:, :, 1:k]) for k in axes(wet3D, 3)]
minindexinlayer = [1; maxindexinlayer[1:(end - 1)] .+ 1]
indices = [min:max for (min, max) in zip(minindexinlayer, maxindexinlayer)]
# Define Preconditioner
struct LayerPreconditioner
    indices
    prob_layers
end
Base.eltype(::LayerPreconditioner) = Float64
@views function LinearAlgebra.ldiv!(Pl::LayerPreconditioner, x::AbstractVector)
    for (idx, prob) in zip(Pl.indices, Pl.prob_layers)
        prob.b .= x[idx]
        solve!(prob)
        x[idx] .= prob.u
    end
    return x
end
@views function LinearAlgebra.ldiv!(y::AbstractVector, Pl::LayerPreconditioner, x::AbstractVector)
    for (idx, prob) in zip(Pl.indices, Pl.prob_layers)
        prob.b .= x[idx]
        solve!(prob)
        y[idx] .= prob.u
    end
    return y
end
@time "MKLPardisoFactorize init" Pl_layer = LayerPreconditioner(
    indices,
    [init(LinearProblem(M[idx, idx], ones(length(idx))), solver, rtol = 1e-10) for idx in indices],
)

prob = LinearProblem(M, ones(N), u0 = ustrip(s, 1000yr) * ones(N))
@time "GMRES + MKLPardisoFactorize(k) solve" sol = solve(prob, KrylovJL_GMRES(); Pl = Pl_layer, maxiters = 500, restarts = 40, verbose = true, reltol = 1.0e-12).u

@show sol_error = norm(M * sol - ones(N)) / norm(ones(N))

# turn the age solution vector back into a 3D cube
agecube = YAXArray(
    dims(volcello),
    ustrip.(yr, OceanTransportMatrixBuilder.as3D(sol, wet3D) * s),
    Dict(
        "description" => "steady-state mean age",
        "solver" => "GMRES + MKLPardisoFactorize-by-layer",
        "error" => string(sol_error),
        # "model" => model,
        # "experiment" => experiment,
        # "member" => member,
        # "time window" => time_window,
        "upwind" => upwind_str2,
        "units" => "yr",
    )
)
arrays = Dict(:age => agecube, :lat => lat, :lon => lon)
ds = Dataset(; properties = Dict(), arrays...)
# Save to netCDF file
outputfile = joinpath(inputdir, "steady_age_$(κVdeep_str)_$(κH_str)_$(κVML_str)_GMRES_PFk.nc")
@info "Saving age as netCDF file:\n  $(outputfile)"
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)

