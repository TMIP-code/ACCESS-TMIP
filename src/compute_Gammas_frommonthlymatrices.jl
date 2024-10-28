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
using NaNStatistics
using Format
using Dates

ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
using PythonCall

# Making monthly matrices for every month from 1990s and averaging them into a single matrix, and computing the age

model = "ACCESS-ESM1-5"
member = "r1i1p1f1"
CMIP_version = "CMIP5"
experiment = "historical"
time_window = "Jan1990-Dec1999"

# 1. get monthly data directly from CMIP archive
intake = pyimport("intake")

catalog = intake.cat.access_nri["cmip6_fs38"]

searched_cat = catalog.search(
    source_id = model,
    member_id = member,
    experiment_id = experiment,
    file_type = "l"
)

# Get umo and vmo for the 1990s
umodf = DataFrame(PyTable(searched_cat.search(variable_id = "umo").df))
umopath = only([x for x in umodf.path if contains(x, "199")]) # Lucky there a single path since saved every decade.
umo_ds = open_dataset(umopath)

vmodf = DataFrame(PyTable(searched_cat.search(variable_id = "vmo").df))
vmopath = only([x for x in vmodf.path if contains(x, "199")]) # Lucky there a single path since saved every decade.
vmo_ds = open_dataset(vmopath)

mlotstdf = DataFrame(PyTable(searched_cat.search(variable_id = "mlotst").df))
mlotst_ds = open_dataset(only(mlotstdf.path))

volcellodf = DataFrame(PyTable(searched_cat.search(variable_id = "volcello", table_id = "Ofx").df)) # Note I am using the "fixed" variable here.
volcello_ds = open_dataset(only(volcellodf.path))

areacellodf = DataFrame(PyTable(searched_cat.search(variable_id = "areacello", table_id = "Ofx").df)) # Note I am using the "fixed" variable here.
areacello_ds = open_dataset(only(areacellodf.path))

# Gadi directory for output files
CMIP6outputdir = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)"

# 2. Load Tilo data (GM + submeso transport)
CSIRO_member = CMIP6member2CSIROmember(member)
Tilodatainputdir = "/g/data/xv83/TMIP/data/$model/$CSIRO_member"
ϕᵢGM_ds = open_dataset(joinpath(Tilodatainputdir, "month_tx_trans_gm_1990s.nc"))
ϕⱼGM_ds = open_dataset(joinpath(Tilodatainputdir, "month_ty_trans_gm_1990s.nc"))
ϕᵢsubmeso_ds = open_dataset(joinpath(Tilodatainputdir, "month_tx_trans_submeso_1990s.nc"))
ϕⱼsubmeso_ds = open_dataset(joinpath(Tilodatainputdir, "month_ty_trans_submeso_1990s.nc"))

# Load fixed variables in memory
areacello = readcubedata(areacello_ds.areacello)
volcello = readcubedata(volcello_ds.volcello)
lon = readcubedata(volcello_ds.longitude)
lat = readcubedata(volcello_ds.latitude)
lev = volcello_ds.lev
# Identify the vertices keys (vary across CMIPs / models)
volcello_keys = propertynames(volcello_ds)
lon_vertices_key = volcello_keys[findfirst(x -> occursin("lon", x) & occursin("vert", x), string.(volcello_keys))]
lat_vertices_key = volcello_keys[findfirst(x -> occursin("lat", x) & occursin("vert", x), string.(volcello_keys))]
lon_vertices = readcubedata(getproperty(volcello_ds, lon_vertices_key))
lat_vertices = readcubedata(getproperty(volcello_ds, lat_vertices_key))

# Make makegridmetrics
gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
(; lon_vertices, lat_vertices) = gridmetrics

# Make indices
indices = makeindices(gridmetrics.v3D)

# Some parameter values
ρ = 1035.0    # kg/m^3
κH = 500.0    # m^2/s
κVML = 0.1    # m^2/s
κVdeep = 1e-5 # m^2/s

monthly_T = [];
# for year in 1990:1999
YEAR = 1990
    # for month in 1:12
    MONTH = 1

        # Load variables in memory
        mlotst = readcubedata(mlotst_ds.mlotst[Ti=Near(Date(YEAR, MONTH, 15))])
        umo = readcubedata(umo_ds.umo[Ti=Near(Date(YEAR, MONTH, 15))])
        vmo = readcubedata(vmo_ds.vmo[Ti=Near(Date(YEAR, MONTH, 15))])

        ϕᵢGM = readcubedata(ϕᵢGM_ds.tx_trans_gm[Ti=Near(Date(YEAR, MONTH, 15))])
        ϕⱼGM = readcubedata(ϕⱼGM_ds.ty_trans_gm[Ti=Near(Date(YEAR, MONTH, 15))])
        ϕᵢsubmeso = readcubedata(ϕᵢsubmeso_ds.tx_trans_submeso[Ti=Near(Date(YEAR, MONTH, 15))])
        ϕⱼsubmeso = readcubedata(ϕⱼsubmeso_ds.ty_trans_submeso[Ti=Near(Date(YEAR, MONTH, 15))])

        # TODO take the vertical diff first!
        WIP WIP

        # TODO fix incompatible dimensions betwewen umo and ϕᵢGM/ϕᵢsubmeso Dim{:i} and Dim{:xu_ocean}
        ϕ = let umo = umo + ϕᵢGM + ϕᵢsubmeso, vmo = vmo + ϕⱼGM + ϕⱼsubmeso
            facefluxesfrommasstransport(; umo, vmo)
        end

        weight = daysinmonth(Date(YEAR, MONTH, 15))

        (; T, Tadv, TκH, TκVML, TκVdeep) = transportmatrix(; ϕ, mlotst, gridmetrics, indices, ρ, κH, κVML, κVdeep)

    end





# 3. Then solve for ideal age / reemergence time

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
    outputfile = joinpath(CMIP6outputdir, "ideal_mean_age_$str.nc")
    @info "Saving ideal mean age as netCDF file:\n  $(outputfile)"
    savedataset(Γin_ds, path = outputfile, driver = :netcdf, overwrite = true)

    # Save Γoutyr3D to netCDF file
    outputfile = joinpath(CMIP6outputdir, "mean_reemergence_time_$str.nc")
    @info "Saving mean reemergence time as netCDF file:\n  $(outputfile)"
    savedataset(Γout_ds, path = outputfile, driver = :netcdf, overwrite = true)


end
