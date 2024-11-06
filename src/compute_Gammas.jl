using Pkg
Pkg.activate(".")
Pkg.instantiate()

using OceanTransportMatrixBuilder
using NetCDF
using YAXArrays
using DataFrames
using DimensionalData
using SparseArrays
using LinearAlgebra
using Unitful
using Unitful: s, yr

# @show model = ARGS[1]
# model = "ACCESS1-0"
# model = "ACCESS1-3"
# model = "ACCESS-CM2"
model = "ACCESS-ESM1-5"
# model = "CESM2"
# model = "CESM2-FV2"
# model = "CESM2-WACCM-FV2"
# model = "CMCC-CM2-HR4"
# model = "CMCC-CM2-SR5"
# model = "CMCC-ESM2"
# model = "FGOALS-f3-L"
# model = "FGOALS-g3"
# model = "MPI-ESM-1-2-HAM"
# model = "MPI-ESM1-2-HR"
# model = "MPI-ESM1-2-LR"
# model = "NorCPM1"
# model = "NorESM2-LM"
# model = "NorESM2-MM"

# CMIP_version = "CMIP5"
CMIP_version = "CMIP6"

experiment = "historical"
# experiment = "piControl"

time_window = "Jan1990-Dec1999"
# time_window = "Jan1071-Dec1100" # <- last 30 years of ACCESS-ESM1-5 piControl
# time_window = "Jan1420-Dec1449" # <- last 30 years of ACCESS-CM2 piControl


# Gadi directory for input files
inputdirfun(member) = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)"

# find all members for which the inputdir contains umo.nc, vmo.nc, mlotst.nc, volcello.nc, and areacello.nc
requiredvariables = ["umo", "vmo", "mlotst", "volcello", "areacello"]
hasrequireddata(member, variable_name) = isfile(joinpath(inputdirfun(member), "$variable_name.nc"))
hasrequireddata(member) = all(variable_name -> hasrequireddata(member, variable_name), requiredvariables)
members = readdir("/scratch/xv83/TMIP/data/$model/$experiment")

# sort members by r, i, p[, f]

member_regex = CMIP_version == "CMIP6" ? r"r(\d+)i(\d+)p(\d+)f(\d+)" : r"r(\d+)i(\d+)p(\d+)"
parse_member(member) = parse.(Int, match(member_regex, member).captures)
members = sort(members, by = x -> parse_member(x))
dataavailability = DataFrame(
    :member => members,
    :has_it_all => hasrequireddata.(members),
    [Symbol(var) => [hasrequireddata(member, var) for member in members] for var in requiredvariables]...,
)
show(dataavailability; allrows = true)

all(.!dataavailability.has_it_all) && @warn "Nothing to do: missing something for all members"

for member in members[dataavailability.has_it_all]
# for member in [first(members)]

    inputdir = inputdirfun(member)

    # Load umo, vmo, mlotst, volcello, and areacello
    umo_ds = open_dataset(joinpath(inputdir, "umo.nc"))
    vmo_ds = open_dataset(joinpath(inputdir, "vmo.nc"))
    mlotst_ds = open_dataset(joinpath(inputdir, "mlotst.nc"))
    volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
    areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))

    # Load variables in memory
    umo = readcubedata(umo_ds.umo)
    vmo = readcubedata(vmo_ds.vmo)
    mlotst = readcubedata(mlotst_ds.mlotst)
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

    # Some parameter values
    ρ = 1035.0    # kg/m^3
    κH = 500.0    # m^2/s
    κVML = 0.1    # m^2/s
    κVdeep = 1e-5 # m^2/s

    # Make makegridmetrics
    gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
    (; lon_vertices, lat_vertices) = gridmetrics

    # Make fuxes from all directions
    ϕ = facefluxesfrommasstransport(; umo, vmo, gridmetrics, indices)

    # Make indices
    indices = makeindices(gridmetrics.v3D)

    # Make transport matrix
    (; T, Tadv, TκH, TκVML, TκVdeep) = transportmatrix(; ϕ, mlotst, gridmetrics, indices, ρ, κH, κVML, κVdeep)

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
    arrays = Dict(:age => Γinyr_YAXArray, :lat => volcello_ds.lat, :lon => volcello_ds.lon)
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
    arrays = Dict(:age => Γoutyr_YAXArray, :lat => volcello_ds.lat, :lon => volcello_ds.lon)
    Γout_ds = Dataset(; volcello_ds.properties, arrays...)

    # Save Γinyr3D to netCDF file
    outputfile = joinpath(inputdir, "ideal_mean_age.nc")
    @info "Saving ideal mean age as netCDF file:\n  $(outputfile)"
    savedataset(Γin_ds, path = outputfile, driver = :netcdf, overwrite = true)

    # Save Γoutyr3D to netCDF file
    outputfile = joinpath(inputdir, "mean_reemergence_time.nc")
    @info "Saving mean reemergence time as netCDF file:\n  $(outputfile)"
    savedataset(Γout_ds, path = outputfile, driver = :netcdf, overwrite = true)

end


