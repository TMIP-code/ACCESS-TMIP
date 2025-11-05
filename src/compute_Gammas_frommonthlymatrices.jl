using Pkg
Pkg.activate(".")
Pkg.instantiate()


ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
using PythonCall
using OceanTransportMatrixBuilder
using NetCDF
using YAXArrays
using DataFrames
using DimensionalData
using SparseArrays
using LinearAlgebra
using Unitful
using Unitful: s, yr
using NaNStatistics
using Format
using Dates
using FileIO

# Load matrix and grid metrics
model = "ACCESS-ESM1-5"
member = "r1i1p1f1"
experiment = "historical"
time_window = "Jan1990-Dec1999"
CMIP6outputdir = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)"
inputfile = joinpath(CMIP6outputdir, "transportmatrix_frommonthlymatrices.jld2")
@info "Loading matrices + metrics as $inputfile"
input = load(inputfile)
T = input["T"]
gridmetrics = input["gridmetrics"]
indices = input["indices"]
note = input["note"]
OceanTransportMatrixBuilder_version = input["OceanTransportMatrixBuilder"]

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
(; lon, lat, zt, v3D) = gridmetrics
lev = zt
# unpack indices
(; wet3D, N) = indices

v = v3D[wet3D]

# surface mask
issrf3D = copy(wet3D)
issrf3D[:, :, 2:end] .= false
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
Γinyr_YAXArray = rebuild(
    volcello_ds["volcello"];
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
Γoutyr_YAXArray = rebuild(
    volcello_ds["volcello"];
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
