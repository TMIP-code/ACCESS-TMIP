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
time_window = "Jan2150-Dec2159"
experiment = "01deg_jra55v13_ryf9091_qian_wthmp"


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
inputdir = "/scratch/y99/TMIP/data/$model/$experiment/$time_window"
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
v = v3D[wet3D]
V = sparse(Diagonal(v))

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
b = ones(N)

# Coarsen 2x2 North of 45°S (so effectively 0.1° -> 0.2°)
@info "coarsening grid 2x2 north of 45°S"
SOmask = lat.data .< -45
mymask = .!SOmask .& trues(size(wet3D))
LUMP, SPRAY, v_c = OceanTransportMatrixBuilder.lump_and_spray(wet3D, v, T, mymask; di = 2, dj = 2, dk = 1)
M_c = LUMP * M * SPRAY
b_c = LUMP * b

matrix_type = Pardiso.REAL_SYM
@show solver = MKLPardisoIterate(; nprocs, matrix_type)

prob = init(LinearProblem(M_c, b_c), solver, rtol = 1.0e-10)
sol_c = solve!(prob).u
sol = SPRAY * sol_c


# turn the age solution vector back into a 3D cube
agecube = DimensionalData.rebuild(
    dzt; # FIXME should be volcello
    data = ustrip.(yr, OceanTransportMatrixBuilder.as3D(sol, wet3D) * s),
    dims = dims(dzt), # FIXME should be volcello
    metadata = Dict(
        "description" => "steady-state mean age",
        "model" => model,
        "experiment" => experiment,
        "time window" => time_window,
        "upwind" => upwind_str2,
        "units" => "yr",
    )
)
arrays = Dict(:age => agecube, :lat => lat, :lon => lon)
ds = Dataset(; dzt_ds.properties, arrays...) # FIXME should be volcello
# Save to netCDF file
outputfile = joinpath(inputdir, "steady_age_$(κVdeep_str)_$(κH_str)_$(κVML_str)_2x245S.nc")
@info "Saving age as netCDF file:\n  $(outputfile)"
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)
