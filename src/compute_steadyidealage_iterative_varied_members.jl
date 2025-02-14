# qsub -I -P xv83 -q hugemem -l mem=360GB -l storage=scratch/gh0+scratch/xv83 -l walltime=02:00:00 -l ncpus=48

using Pkg
Pkg.activate(".")
Pkg.instantiate()
const nprocs = 24

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

# The tracer equation is
#
#   ∂x(t)/∂t + T(t) x(t) = s(t) - Ω x(t)
#
# where Ω "relaxes" x to zero in the top layer.
#
# Applying Backward Euler time step gives
#
#   (I + δt M(t+δt)) x(t+δt) = x(t) + δt s(t+δt)
#
# where M(t) = T(t) + Ω.
#
# More succintly, if m represents the months (1..12)
#
#   (I + δt Mₘ₊₁) xₘ₊₁ = xₘ + δt sₘ₊₁          (1)
#
# Here δt = δt(m-1..m) is the time that separates
# the "center time" of climatological months m-1 and m.
# So the δt that multiplies Mₘ is δ(m..m+1).



# script options
@show model = "ACCESS-ESM1-5"
if isempty(ARGS)
    member = "r20i1p1f1"
    experiment = "historical"
    time_window = "Jan1850-Dec1859"
    # time_window = "Jan1990-Dec1999"
    # experiment = "ssp370"
    # time_window = "Jan2030-Dec2039"
    # time_window = "Jan2090-Dec2099"
    average_over = "monthlymatrices"
else
    experiment, member, time_window, average_over = ARGS
end
κVML = 1e-7   # m^2/s
κVdeep = 1e-7 # m^2/s
κH = 5        # m^2/s
@show experiment
@show member
@show time_window
@show κVdeep
@show κH
@show κVML
@show average_over

κVdeep_str = "kVdeep" * format(κVdeep, conversion="e")
κH_str = "kH" * format(κH, conversion="d")
κVML_str = "kVML" * format(κVML, conversion="e")

lumpby = "month"
months = 1:12
Nmonths = length(months)

# Load areacello and volcello for grid geometry
fixedvarsinputdir = "/scratch/xv83/TMIP/data/$model"
volcello_ds = open_dataset(joinpath(fixedvarsinputdir, "volcello.nc"))
areacello_ds = open_dataset(joinpath(fixedvarsinputdir, "areacello.nc"))

# Gadi directory for input files
inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)"
cycloinputdir = joinpath(inputdir, "cyclo$lumpby")

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

issrf = let
    issrf3D = falses(size(wet3D))
    issrf3D[:,:,1] .= true
    issrf3D[wet3D]
end
Ω = sparse(Diagonal(Float64.(issrf)))


M̄ = if average_over == "monthlymatrices"
    # Load monthly matrices
    @time "Loading monthly matrices" Ms = map(months) do m
        inputfile = joinpath(cycloinputdir, "cyclo_matrix_$(κVdeep_str)_$(κH_str)_$(κVML_str)_$m.jld2")
        @info "  Loading T from $inputfile"
        T = load(inputfile, "T")
        T + Ω
    end
    mean(Ms)
elseif average_over == "flow"
    # Or load the mean flow matrix
    inputfile = joinpath(cycloinputdir, "cyclo_matrix_$(κVdeep_str)_$(κH_str)_$(κVML_str)_meanflow.jld2")
    @info "Loading T from $inputfile"
    T = load(inputfile, "T")
    T + Ω
end

@time "steady state solve" u0 = solve(LinearProblem(M̄, ones(N)), MKLPardisoIterate(; nprocs), rtol = 1e-10).u

# Save cyclo-stationary age
Γinyr = ustrip.(yr, u0 .* s)
Γinyr3D = OceanTransportMatrixBuilder.as3D(Γinyr, wet3D)

cube3D = rebuild(volcello_ds["volcello"];
    data = Γinyr3D,
    dims = dims(volcello_ds["volcello"]),
    metadata = Dict(
        "origin" => "steady-state ideal age from mean $average_over computed from $model $experiment $member $(time_window)",
        "units" => "yr",
    )
)

arrays = Dict(:age => cube3D, :lat => volcello_ds.lat, :lon => volcello_ds.lon)
ds = Dataset(; volcello_ds.properties, arrays...)

# Save Γinyr3D to netCDF file
outputfile = joinpath(cycloinputdir, "steady_state_ideal_mean_age_$(κVdeep_str)_$(κH_str)_$(κVML_str)_mean$(average_over).nc")
@info "Saving ideal mean age as netCDF file:\n  $(outputfile)"
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)

