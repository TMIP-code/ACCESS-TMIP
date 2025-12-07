# # qsub -I -P y99 -q normal -l mem=190GB -l storage=scratch/y99+gdata/xp65 -l walltime=02:00:00 -l ncpus=48

# using Pkg
# Pkg.activate(".")
# Pkg.instantiate()

# ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
# using OceanTransportMatrixBuilder
# using NetCDF
# using YAXArrays
# using DataFrames
# using DimensionalData
# using SparseArrays
# using LinearAlgebra
# using Unitful
# using Unitful: s, yr, d
# using Statistics
# using Format
# using Dates
# using FileIO
# # using LinearSolve
# # using IncompleteLU
# # using NonlinearSolve
# using ProgressMeter
# using Krylov
# using KrylovPreconditioners
# using SparseMatricesCSR
# using ThreadedSparseCSR
# using InteractiveUtils: @which

# ThreadedSparseCSR.multithread_matmul(BaseThreads())

# model = "ACCESS-OM2-1"
# experiment = "1deg_jra55_iaf_omip2_cycle6"
# time_window = "Jan1960-Dec1979"
# @show inputdir = "/scratch/y99/TMIP/data/$model/$experiment/$time_window"

# # preferred diffusivities
# κVdeep = 3.0e-5 # m^2/s
# κVML = 1.0      # m^2/s
# κH = 300.0      # m^2/s (grid-scaling by sqrt(area))
# @show κVdeep
# @show κVML
# @show κH
# κVdeep_str = "kVdeep" * format(κVdeep, conversion = "e")
# κVML_str = "kVML" * format(κVML, conversion = "e")
# κH_str = "kH" * format(κH, conversion = "d")

# upwind = false
# @show upwind
# upwind_str = upwind ? "" : "_centered"
# upwind_str2 = upwind ? "upwind" : "centered"

# # Load areacello and volcello for grid geometry
# areacello_ds = open_dataset(joinpath(inputdir, "area_t.nc"))
# dht_ds = open_dataset(joinpath(inputdir, "dht.nc")) # <- (new) cell thickness?
# # TODO: caputre correlations between transport and dht
# # z* coordinate varies with time in ACCESS-OM2
# # volcello_ds = open_dataset(joinpath(fixedvarsinputdir, "volcello.nc")) # <- not in ACCESS-OM2; must be built from dht * area


# # Load fixed variables in memory
# areacello_OM2 = replace(readcubedata(areacello_ds.area_t), missing => NaN) # This is required for later multiplication
# dht = readcubedata(dht_ds.dht)
# lon_OM2 = readcubedata(areacello_ds.geolon_t)
# lat_OM2 = readcubedata(areacello_ds.geolat_t)
# lev = dht_ds.st_ocean

# # Unfortunately ACCESS-OM2 raw data does not have coordinates of cell vertices
# # So instead I go back to the source: the supergrids
# include("supergrid.jl")
# (; lon, lat, areacello, lon_vertices, lat_vertices) = supergrid(model; dims = dims(dht_ds.dht)[1:2])

# @show size(lon_vertices)

# volcello = readcubedata(dht .* areacello)

# # Make makegridmetrics
# gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
# (; lon_vertices, lat_vertices, v3D) = gridmetrics

# # Make indices
# indices = makeindices(v3D)
# (; wet3D) = indices

# # Make V diagnoal matrix of volumes
# V = sparse(Diagonal(v3D[wet3D]))

# issrf = let
#     issrf3D = falses(size(wet3D))
#     issrf3D[:, :, 1] .= true
#     issrf3D[wet3D]
# end
# idx_surface = findall(issrf)
# idx_interior = findall(.!issrf)
# Nᵢ = length(idx_interior)
# Nₛ = length(idx_surface)

# TMfile = joinpath(inputdir, "yearly_matrix_$(κVdeep_str)_$(κH_str)_$(κVML_str).jld2")
# @info "Loading matrix from $TMfile"
# T = load(TMfile, "T")
# function showsparsitypattern(str, A)
#     println(str)
#     display(A)
#     println()
# end
# showsparsitypattern("Transport matrix T =", T)

Tᵢᵢ = T[idx_interior, idx_interior]
showsparsitypattern("Tᵢᵢ =", Tᵢᵢ)

println("queue: ", ENV["PBS_O_QUEUE"])
println("VMEM: ", parse(Int, ENV["PBS_VMEM"]) ÷ 10^9, "GB")
println("NCPUS: ", parse(Int, ENV["PBS_NCPUS"]))

@show τ = parse(Float64, get(ENV, "ILUTAU", "1e-8")) # drop tolerance for ILU
@time "ILU preconditioner" Pr = KrylovPreconditioners.ilu(Tᵢᵢ; τ)

showsparsitypattern("Incomplete U =", Pr.U)
showsparsitypattern("Incomplete L =", Pr.L)

Tᵢᵢ = sparsecsr(findnz(Tᵢᵢ)..., size(Tᵢᵢ)...)

@time "BICGSTAB + ILU" Γꜜᵢ, stats = Krylov.bicgstab(Tᵢᵢ, ones(Nᵢ);
    N = Pr,
    ldiv = true,
    rtol = 1e-10,
    itmax = 2000,
    timemax = 360.0,
    verbose = 1,
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
        "solver" => "BICGSTAB + ILU",
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
outputfile = joinpath(inputdir, "steady_age_$(κVdeep_str)_$(κH_str)_$(κVML_str)_BICGSTAB_ILU.nc")
@info "Saving age as netCDF file:\n  $(outputfile)"
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)

