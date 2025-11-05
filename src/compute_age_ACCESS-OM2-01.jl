# qsub -I -P y99 -q megamem -l mem=2990GB -l storage=scratch/gh0+scratch/y99+scratch/p66 -l walltime=1:00:00 -l ncpus=48

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

model = "ACCESS-OM2-01"


# preferred diffusivities
κVdeep = 3.0e-5 # m^2/s
κVML = 1.0      # m^2/s
κH = 300.0 / 10 # m^2/s (grid-scaling by sqrt(area))
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
# But not available in raw output, so must rebuild from area_t and dzt
inputdir = "/scratch/y99/TMIP/data/$model/01deg_jra55v13_ryf9091_qian_wthmp/Jan2150-Dec2159"
area_t_ds = open_dataset(joinpath(inputdir, "area_t.nc"))
dzt_ds = open_dataset(joinpath(inputdir, "dzt.nc"))
areacello = readcubedata(area_t_ds.area_t)
dzt = readcubedata(dzt_ds.dzt)
volcello = areacello .* dzt

# Load fixed variables in memory
# Identify the vertices keys (vary across CMIPs / models)
# FIXME using test CMORized data for vertices. Hopefully they match!
CMORtestfile = "/scratch/p66/yz9299/OM2_CMORised/umo_Omon_$(model)_omip1_r1i1p1f1_gn_214501-214503.nc"
CMORtest_ds = open_dataset(CMORtestfile)
lon = readcubedata(CMORtest_ds.longitude)
lat = readcubedata(CMORtest_ds.latitude)
lev = CMORtest_ds.lev
volcello_keys = propertynames(CMORtest_ds)
lon_vertices_key = volcello_keys[findfirst(x -> occursin("lon", x) & occursin("vert", x), string.(volcello_keys))]
lat_vertices_key = volcello_keys[findfirst(x -> occursin("lat", x) & occursin("vert", x), string.(volcello_keys))]
lon_vertices = readcubedata(getproperty(CMORtest_ds, lon_vertices_key))
lat_vertices = readcubedata(getproperty(CMORtest_ds, lat_vertices_key))


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
@show solver = MKLPardisoIterate(; nprocs, matrix_type)

prob = init(LinearProblem(M, ones(N)), solver, rtol = 1.0e-10)
sol = solve!(prob).u


# turn the age solution vector back into a 3D cube
agecube = DimensionalData.rebuild(
    volcello_ds["volcello"];
    data = ustrip.(yr, OceanTransportMatrixBuilder.as3D(sol, wet3D) * s),
    dims = dims(volcello),
    metadata = Dict(
        "description" => "steady-state mean age",
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
outputfile = joinpath(inputdir, "steady_age_$(κVdeep_str)_$(κH_str)_$(κVML_str).nc")
@info "Saving age as netCDF file:\n  $(outputfile)"
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)
