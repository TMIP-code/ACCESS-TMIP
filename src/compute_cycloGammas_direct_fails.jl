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

v = v3D[wet3D];

issrf = let
    issrf3D = zeros(size(wet3D))
    issrf3D[:, :, 1] .= 1
    issrf3D[wet3D]
end
M = sparse(Diagonal(issrf))

# Now let's imagine we have 4 seasons of OCCA where T gets more or less intense
Nseasons = length(seasons)
δt = ustrip(s, 365d / Nseasons) # TODO maybe use exact mean number of days (more important for monthly because Feb)?
A = BlockArray(spzeros(Nseasons * N, Nseasons * N), fill(N, Nseasons), fill(N, Nseasons))
@time "building blocks" for (i, season) in enumerate(seasons)
    inputfile = joinpath(cycloinputdir, "cyclo_matrix_$season.jld2")
    @info "Loading matrices + metrics as $inputfile"
    T = load(inputfile)["T"]
    A[Block(i, i)] = I + δt * (T + M)
    A[Block(mod1(i + 1, Nseasons), i)] = -I(N)
end
# @time "converting BlockArray to standard sparse" A = sparse(A)
A = @time "converting BlockArray to standard sparse" reduce(hcat, reduce(vcat, A[Block(i, j)] for i in 1:blocksize(A, 1)) for j in 1:blocksize(A, 2))

b = δt * ones(size(A, 1))

# Solve the system using the LinearSolve.jl package
prob = LinearProblem(A, b)

# extra test
nprocs = 48
@time "init" linsolve = init(prob, MKLPardisoFactorize(; nprocs))
@info "Now attempting seasonal solve"

# Line below fails for 180GB memory
# @time "solve" solve!(linsolve)

foo
