
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
# using LinearSolve
# using IncompleteLU
# using NonlinearSolve
using ProgressMeter
using Krylov
using KrylovPreconditioners
using SparseMatricesCSR
using ThreadedSparseCSR
using InteractiveUtils: @which

ThreadedSparseCSR.multithread_matmul(BaseThreads())

model = "ACCESS-OM2-01"
experiment = "01deg_jra55v140_iaf_cycle4"
time_window = "Jan1960-Dec1979"
@show inputdir = "/scratch/y99/TMIP/data/$model/$experiment/$time_window"

# preferred diffusivities
κVdeep = 3.0e-5 # m^2/s
κVML = 1.0      # m^2/s
κH = 300.0 / 4  # m^2/s (grid-scaling by sqrt(area))
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
dzt_ds = open_dataset(joinpath(inputdir, "dzt.nc")) # <- (new) cell thickness?

# Load fixed variables in memory
dzt = readcubedata(dzt_ds.dzt)
lev = dzt_ds.st_ocean

# Unfortunately ACCESS-OM2 raw data does not have coordinates of cell vertices
# So instead I go back to the source: the supergrids
include("supergrid.jl")
(; lon, lat, areacello, lon_vertices, lat_vertices) = supergrid(model; dims = dims(dzt_ds.dzt)[1:2])

@show size(lon_vertices)

volcello = readcubedata(dzt .* areacello)

# Make makegridmetrics
gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
(; lon_vertices, lat_vertices, v3D) = gridmetrics

# Make indices
indices = makeindices(v3D)
(; wet3D) = indices

# Make V diagnoal matrix of volumes
V = sparse(Diagonal(v3D[wet3D]))

issrf = let
    issrf3D = falses(size(wet3D))
    issrf3D[:, :, 1] .= true
    issrf3D[wet3D]
end
idx_surface = findall(issrf)
idx_interior = findall(.!issrf)
Nᵢ = length(idx_interior)
Nₛ = length(idx_surface)

TMfile = joinpath(inputdir, "yearly_matrix_$(κVdeep_str)_$(κH_str)_$(κVML_str).jld2")
@info "Loading matrix from $TMfile"
T = load(TMfile, "T")
function showsparsitypattern(str, A)
    println(str)
    display(A)
    println()
end
showsparsitypattern("Transport matrix T =", T)

Tᵢᵢ = T[idx_interior, idx_interior]
showsparsitypattern("Tᵢᵢ =", Tᵢᵢ)

println("queue: ", ENV["PBS_O_QUEUE"])
println("VMEM: ", parse(Int, ENV["PBS_VMEM"]) ÷ 10^9, "GB")
println("NCPUS: ", parse(Int, ENV["PBS_NCPUS"]))

@show τ = parse(Float64, ENV["ILUTAU"]) # drop tolerance for ILU
@time "ILU preconditioner" Pl = KrylovPreconditioners.ilu(Tᵢᵢ; τ)  # Block-Jacobi preconditioner
# @time "GMRES + ILU" Γꜜᵢ, stats = Krylov.gmres(Tᵢᵢ, ones(Nᵢ);
#     M = Pl, # linear operator that models a nonsingular matrix of size n used for left preconditioning;
#     memory = 40, # if restart = true, the restarted version GMRES(k) is used with k = memory. If restart = false, the parameter memory should be used as a hint of the number of iterations to limit dynamic memory allocations. Additional storage will be allocated if the number of iterations exceeds memory;
#     # N = # linear operator that models a nonsingular matrix of size n used for right preconditioning;
#     ldiv = true, # define whether the preconditioners use ldiv! or mul!;
#     restart = true, # restart the method after memory iterations;
#     reorthogonalization = true, # reorthogonalize the new vectors of the Krylov basis against all previous vectors;
#     # atol = # absolute stopping tolerance based on the residual norm;
#     rtol = 1e-10, # relative stopping tolerance based on the residual norm;
#     itmax = 500, # the maximum number of iterations. If itmax=0, the default number of iterations is set to 2n;
#     timemax = 360.0, # the time limit in seconds;
#     verbose = 10, # additional details can be displayed if verbose mode is enabled (verbose > 0). Information will be displayed every verbose iterations;
#     # history = # collect additional statistics on the run such as residual norms, or Aᴴ-residual norms;
#     # callback = # function or functor called as callback(workspace) that returns true if the Krylov method should terminate, and false otherwise;
#     # iostream = # stream to which output is logged.
# )
showsparsitypattern("Incomplete U =", Pl.U)
showsparsitypattern("Incomplete L =", Pl.L)

Tᵢᵢ = sparsecsr(findnz(Tᵢᵢ)..., size(Tᵢᵢ)...)

@time "GMRES + ILU" Γꜜᵢ, stats = Krylov.gmres(Tᵢᵢ, ones(Nᵢ);
    M = Pl,
    memory = 40,
    ldiv = true,
    restart = true,
    reorthogonalization = true,
    rtol = 1e-10,
    itmax = 2000,
    # timemax = 360.0,
    verbose = 10,
)

@show stats

# Check error magnitude
@show sol_error = norm(Tᵢᵢ * Γꜜᵢ - ones(Nᵢ)) / norm(ones(Nᵢ))

# turn the age solution vector back into a 3D yax
Γꜜyax = YAXArray(
    dims(volcello),
    ustrip.(yr, OceanTransportMatrixBuilder.as3D([zeros(Nₛ); Γꜜᵢ], wet3D) * s),
    Dict(
        "description" => "steady-state mean age",
        "solver" => "GMRES + ILU",
        "model" => model,
        "experiment" => experiment,
        "time window" => time_window,
        "upwind" => upwind_str2,
        "units" => "yr",
    )
)
arrays = Dict(:Gammadown => Γꜜyax, :lat => lat, :lon => lon)
ds = Dataset(; properties = Dict(), arrays...)
# Save to netCDF file
outputfile = joinpath(inputdir, "steady_age_$(κVdeep_str)_$(κH_str)_$(κVML_str)_KrylovGMRES_ILU.nc")
@info "Saving age as netCDF file:\n  $(outputfile)"
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)

