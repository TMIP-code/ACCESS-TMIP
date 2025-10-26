# qsub -I -P y99 -q nornal -l mem=190GB -l storage=scratch/gh0+scratch/y99+scratch/p66 -l walltime=10:00:00 -l ncpus=48

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
using ProgressMeter
using MPI
using PETSc

MPI.Initialized() || MPI.Init()
PETSc.initialize()

# script options
@show model = "ACCESS-OM2-025"
if isempty(ARGS)
    experiment = "omip2"
    member = "r1i1p1f1"
    time_window = "Jan0200-Dec0209"
else
    experiment, member, time_window = ARGS
end

@show experiment
@show member
@show time_window

# preferred diffusivities
κVdeep = 3.0e-5 # m^2/s
κVML = 1.0      # m^2/s
κH = 300.0 / 4  # m^2/s (grid-scaling by sqrt(area))
@show κVdeep
@show κVML
@show κH
κVdeep_str = "kVdeep" * format(κVdeep, conversion="e")
κVML_str = "kVML" * format(κVML, conversion="e")
κH_str = "kH" * format(κH, conversion="d")


upwind = false
@show upwind
upwind_str = upwind ? "" : "_centered"
upwind_str2 = upwind ? "upwind" : "centered"

# Load areacello and volcello for grid geometry
inputdir = "/scratch/y99/TMIP/data/ACCESS-OM2-025/omip2/r1i1p1f1/Jan0200-Dec0209/"
areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))
volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))

# Load fixed variables in memory
areacello = readcubedata(areacello_ds.areacello)
volcello = readcubedata(volcello_ds.volcello)
lon = readcubedata(volcello_ds.lon)
lat = readcubedata(volcello_ds.lat)
lev = volcello_ds.lev
# Identify the vertices keys (vary across CMIPs / models)
# FIXME using test CMORized data for vertices. Hopefully they match!
CMORtestfile = "/scratch/p66/yz9299/OM2_CMORised/umo_Omon_ACCESS-OM2-025_omip1_r1i1p1f1_gn_190001-190112.nc"
CMORtest_ds = open_dataset(CMORtestfile)
volcello_keys = propertynames(CMORtest_ds)
lon_vertices_key = volcello_keys[findfirst(x -> occursin("lon", x) & occursin("vert", x), string.(volcello_keys))]
lat_vertices_key = volcello_keys[findfirst(x -> occursin("lat", x) & occursin("vert", x), string.(volcello_keys))]
lon_vertices = readcubedata(getproperty(CMORtest_ds, lon_vertices_key))
lat_vertices = readcubedata(getproperty(CMORtest_ds, lat_vertices_key))


# Make makegridmetrics
gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
(; lon_vertices, lat_vertices, v3D, lat, lon, zt) = gridmetrics

# Make indices
indices = makeindices(v3D)
(; N, wet3D) = indices

# Make V diagnoal matrix of volumes
v = v3D[wet3D]
V = sparse(Diagonal(v))

issrf = let
    issrf3D = falses(size(wet3D))
    issrf3D[:,:,1] .= true
    issrf3D[wet3D]
end
Ω = sparse(Diagonal(Float64.(issrf)))

TMfile = joinpath(inputdir, "yearly_matrix_$(κVdeep_str)_$(κH_str)_$(κVML_str).jld2")
@info "Loading matrix from $TMfile"
T = load(TMfile, "T")
@info "Matrix size: $(size(T)), nnz = $(nnz(T))"

M = T + Ω
b = ones(N)

# Solve the linear system using PETSc + ILU
# don't print these below
M_PETSc = PETSc.MatSeqAIJ(M);
b_PETSc = PETSc.VecSeq(b);
x_PETSc = PETSc.VecSeq(zeros(size(b)));

ksp = PETSc.KSP(
    M_PETSc;
    ksp_monitor = false,    # set to true for output
    ksp_monitor_true_residual = false,
    ksp_view = false,
    ksp_type = "gmres",
    ksp_atol = 1e-10,
    pc_type = "ilu",
);
PETSc.solve!(x_PETSc, ksp, b_PETSc);
sol = x_PETSc.array

# turn the age solution vector back into a 3D cube
agecube = DimensionalData.rebuild(volcello_ds["volcello"];
    data = ustrip.(yr, OceanTransportMatrixBuilder.as3D(sol, wet3D) * s),
    dims = dims(volcello),
    metadata = Dict(
        "description" => "steady-state mean age with PETSc",
        "model" => model,
        "experiment" => experiment,
        "member" => member,
        "time window" => time_window,
        "upwind" => upwind_str2,
        "units" => "yr",
    )
)
arrays = Dict(:age => agecube, :lat => volcello_ds.lat, :lon => volcello_ds.lon)
ds = Dataset(; volcello_ds.properties, arrays...)
# Save to netCDF file
outputfile = joinpath(inputdir, "steady_age_$(κVdeep_str)_$(κH_str)_$(κVML_str)_OM2-025_PETSc.nc")
@info "Saving age as netCDF file:\n  $(outputfile)"
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)
