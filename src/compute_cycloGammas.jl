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
using NaNStatistics
using Format
using Dates
using FileIO
using BlockArrays
using LinearSolve
import Pardiso # import Pardiso instead of using (to avoid name clash?)



# Load matrix and grid metrics
model = "ACCESS-ESM1-5"
member = "r1i1p1f1"
experiment = "historical"
time_window = "Jan1990-Dec1999"
lumpby = "season"
seasons = ("DJF", "MAM", "JJA", "SON")

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
#   ∂ₜx + T x = 1 - M x
# Applying Backward Euler time step gives


issrf = let
    issrf3D = zeros(size(wet3D))
    issrf3D[:,:,1] .= 1
    issrf3D[wet3D]
end
M = sparse(Diagonal(issrf))

# Preconditioner for cyclostationary solver
struct CycloPreconditioner
    DJF # Pardiso factors for DJF
    MAM # Pardiso factors for MAM
    JJA # Pardiso factors for JJA
    SON # Pardiso factors for SON
end
Nseasons = length(seasons)
δt = ustrip(s, 365d / Nseasons) # TODO maybe use exact mean number of days (more important for monthly because Feb)?


# Build matrix to multiply
A = BlockArray(spzeros(Nseasons * N, Nseasons * N), fill(N, Nseasons), fill(N, Nseasons))
@time "building blocks" for (i, season) in enumerate(seasons)
    inputfile = joinpath(cycloinputdir, "cyclo_matrix_$season.jld2")
    @info "Loading matrices + metrics as $inputfile"
    T = load(inputfile)["T"]
    A[Block(i, i)] = I + δt * (T + M)
    A[Block(mod1(i+1, Nseasons), i)] = -I(N)
end
# @time "converting BlockArray to standard sparse" A = sparse(A)
function initseasonlinprob(A, i)
    prob = LinearProblem(A[Block(i,i)], δt * ones(N))
    return init(prob, MKLPardisoFactorize(; nprocs = 48))
end
Pl = CycloPreconditioner(
    initseasonlinprob(A, 1),
    initseasonlinprob(A, 2),
    initseasonlinprob(A, 3),
    initseasonlinprob(A, 4),
)


Base.eltype(::CycloPreconditioner) = Float64
function LinearAlgebra.ldiv!(Pl::CycloPreconditioner, x::AbstractVector)
    # Making a blocked view into x to help with blocks
    xᵇ = BlockedVector(x, fill(N, Nseasons))
    for (i, block) in enumerate(propertynames(Pl))
        # Grab the $block seasonal linear problem
        linprob = getproperty(Pl, block)
        # update the RHS
        if i == 1 # First season is independent of other seasons
            linprob.b = xᵇ[Block(i)]
        else # Other seasons depend on the previous one
            linprob.b = xᵇ[Block(i)] + xᵇ[Block(i - 1)]
        end
        # @time "solve $block" solve!(linprob)
        solve!(linprob)
        # Update input vector
        xᵇ[Block(i)] .= linprob.u
    end
    return x
end
function LinearAlgebra.ldiv!(y::AbstractVector, Pl::CycloPreconditioner, x::AbstractVector)
    # Making blocked views into x and y to help with blocks
    xᵇ = BlockedVector(x, fill(N, Nseasons)) # a blocked view into x to help with blocks
    yᵇ = BlockedVector(y, fill(N, Nseasons)) # a blocked view into x to help with blocks
    for (i, block) in enumerate(propertynames(Pl))
        # Grab the $block seasonal linear problem
        linprob = getproperty(Pl, block)
        # update the RHS
        if i == 1 # First season is independent of other seasons
            linprob.b = xᵇ[Block(i)]
        else # Other seasons depend on the previous one
            linprob.b = xᵇ[Block(i)] + xᵇ[Block(i - 1)]
        end
        # @time "solve $block" solve!(linprob)
        solve!(linprob)
        # Update input vector
        yᵇ[Block(i)] .= linprob.u
    end
    return y
end

@time "converting BlockArray to standard sparse" A = reduce(hcat, reduce(vcat, A[Block(i, j)] for i in 1:blocksize(A, 1)) for j in 1:blocksize(A, 2))
b = δt * ones(N * Nseasons)

@info "Setting up full seasonal problem"
prob = LinearProblem(A, b)

@time "initialize full problem" linsolve = init(prob, KrylovJL_GMRES(gmres_restart = 40), Pl = Pl)

@info "Now attempting seasonal solve"

@time "solve" solve!(linsolve)





FOOFOFOF




# Load volcello to use as the basis for building YAXArrays
# TODO: Learn how to build them from scratch to fix that
intake = pyimport("intake")
catalog = intake.cat.access_nri["cmip6_fs38"]
searched_cat = catalog.search(
    source_id = model,
    member_id = member,
    experiment_id = experiment,
    file_type = "l"
)
volcellodf = DataFrame(PyTable(searched_cat.search(variable_id = "volcello", table_id = "Ofx").df)) # Note I am using the "fixed" variable here.
volcello_ds = open_dataset(only(volcellodf.path))
volcello = readcubedata(volcello_ds.volcello)

# solve for ideal age / reemergence time

# unpack model grid
(; lon, lat, zt, v3D,) = gridmetrics
lev = zt
# unpack indices
(; wet3D, N) = indices

v = v3D[wet3D]

# surface mask
issrf3D = copy(wet3D)
issrf3D[:,:,2:end] .= false
issrf = issrf3D[wet3D]
# Ideal mean age Γ↓ is governed by
# 	∂Γ↓/∂t + T Γ↓ = 1 - M Γ↓
# so
# 	Γ↓ = A⁻¹ 1
# where A = T + M and M is matrix mask of surface with short timescale (1s)
sΓin = ones(size(v))
M = sparse(Diagonal(issrf))
@info "Factorizing A = T + M"
A = factorize(T + M)
@info "Solving ideal mean age"
Γin = A \ sΓin
Γinyr = ustrip.(yr, Γin .* s)
Γinyr3D = OceanTransportMatrixBuilder.as3D(Γinyr, wet3D)

# global mean of the mean age should be order 1000yr
(v' * Γinyr) / sum(v)

# Mean reemergence time is governed by
# 	-∂Γ↑/∂t + (T̃ + M) Γ↑ = 1    (note the tilde on T̃ = V⁻¹ Tᵀ V, the adjoint of T)
# For computational efficiency, reuse factors of A = T + M:
#   Γ↑ = V⁻¹ A⁻ᵀ V 1
V = sparse(Diagonal(v))
V⁻¹ = sparse(Diagonal(1 ./ v))
Γout = V⁻¹ * (transpose(A) \ (V * sΓin))
Γoutyr = ustrip.(yr, Γout .* s)
Γoutyr3D = OceanTransportMatrixBuilder.as3D(Γoutyr, wet3D)

# global mean of the mean reemergence time should also be order 1000yr
(v' * Γout) / sum(v)

# Turn Γin into a YAXArray by rebuilding from volcello
Γinyr_YAXArray = rebuild(volcello_ds["volcello"];
    data = Γinyr3D,
    dims = dims(volcello_ds["volcello"]),
    metadata = Dict(
        "origin" => "ideal_mean_age computed from $model $experiment $member $(time_window)",
        "units" => "yr",
    )
)
arrays = Dict(:age => Γinyr_YAXArray, :lat => volcello_ds.latitude, :lon => volcello_ds.longitude)
Γin_ds = Dataset(; volcello_ds.properties, arrays...)

# Turn Γout into a YAXArray by rebuilding from volcello
Γoutyr_YAXArray = rebuild(volcello_ds["volcello"];
    data = Γoutyr3D,
    dims = dims(volcello_ds["volcello"]),
    metadata = Dict(
        "origin" => "mean_reemergence_time computed from $model $experiment $member $(time_window)",
        "units" => "yr",
    )
)
arrays = Dict(:age => Γoutyr_YAXArray, :lat => volcello_ds.latitude, :lon => volcello_ds.longitude)
Γout_ds = Dataset(; volcello_ds.properties, arrays...)

str = "monthlymatrices"

# Save Γinyr3D to netCDF file
outputfile = joinpath(CMIP6outputdir, "ideal_mean_age_frommonthlymatrices.nc")
@info "Saving ideal mean age as netCDF file:\n  $(outputfile)"
savedataset(Γin_ds, path = outputfile, driver = :netcdf, overwrite = true)

# Save Γoutyr3D to netCDF file
outputfile = joinpath(CMIP6outputdir, "mean_reemergence_time_frommonthlymatrices.nc")
@info "Saving mean reemergence time as netCDF file:\n  $(outputfile)"
savedataset(Γout_ds, path = outputfile, driver = :netcdf, overwrite = true)



