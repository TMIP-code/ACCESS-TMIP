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
using FileIO
using Dates

# Load functions for GM terms
include("GentMcWilliams.jl")

model = "ACCESS-ESM1-5"
CMIP_version = "CMIP5"
experiment = "historical"
time_window = "Jan1990-Dec1999"

# members = [1, 3, 4, 5, 6, 7, 8]
members = [1, 3, 4]
members = [5, 7]
# member = "r1i1p1f1"

for member in members

    CMIP6_member = "r$(member)i1p1f1"
    @info "member $CMIP6_member"

    # Gadi directory for input files
    inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/$CMIP6_member/$(time_window)"

    @info "Loading files lazily"
    umo_ds = open_dataset(joinpath(inputdir, "umo.nc"))
    vmo_ds = open_dataset(joinpath(inputdir, "vmo.nc"))
    mlotst_ds = open_dataset(joinpath(inputdir, "mlotst.nc"))
    volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
    areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))
    tx_trans_gm_ds = open_dataset(joinpath(inputdir, "tx_trans_gm.nc"))
    ty_trans_gm_ds = open_dataset(joinpath(inputdir, "ty_trans_gm.nc"))
    tx_trans_submeso_ds = open_dataset(joinpath(inputdir, "tx_trans_submeso.nc"))
    ty_trans_submeso_ds = open_dataset(joinpath(inputdir, "ty_trans_submeso.nc"))

    @info "Loading data in memory"
    umo = readcubedata(umo_ds.umo)
    vmo = readcubedata(vmo_ds.vmo)
    tx_trans_gm = readcubedata(tx_trans_gm_ds.tx_trans_gm)
    ty_trans_gm = readcubedata(ty_trans_gm_ds.ty_trans_gm)
    tx_trans_submeso = readcubedata(tx_trans_submeso_ds.tx_trans_submeso)
    ty_trans_submeso = readcubedata(ty_trans_submeso_ds.ty_trans_submeso)
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

    @info "Taking vertical diff of GM and submeso terms"
    (nx, ny, _) = size(tx_trans_gm)
    ϕᵢ_gm = diff([fill(0.0, nx, ny, 1);;; tx_trans_gm], dims = 3)
    ϕⱼ_gm = diff([fill(0.0, nx, ny, 1);;; ty_trans_gm], dims = 3)
    ϕᵢ_submeso = diff([fill(0.0, nx, ny, 1);;; tx_trans_submeso], dims = 3)
    ϕⱼ_submeso = diff([fill(0.0, nx, ny, 1);;; ty_trans_submeso], dims = 3)

    # Some parameter values
    ρ = 1035.0    # kg/m^3
    κH = 500.0    # m^2/s
    κVML = 0.1    # m^2/s
    κVdeep = 1e-5 # m^2/s

    @info "Making grid metrics"
    gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
    (; lon_vertices, lat_vertices) = gridmetrics

    @info "Making indices"
    indices = makeindices(gridmetrics.v3D)

    @info "Adding GM + submeso terms to umo and vmo"
    umo = umo + ϕᵢ_gm + ϕᵢ_submeso
    vmo = vmo + ϕⱼ_gm + ϕⱼ_submeso
    str = "resolved_GM_submeso"

    @info "Build resolved + parameterized fluxes across all faces"
    ϕ = facefluxesfrommasstransport(; umo, vmo, gridmetrics, indices)

    @info "Build transport matrix"
    (; T, Tadv, TκH, TκVML, TκVdeep) = transportmatrix(; ϕ, mlotst, gridmetrics, indices, ρ, κH, κVML, κVdeep)

    @info "Save matrices"
    outputfile = joinpath(inputdir, "transportmatrix_$str.jld2")
    @info "Saving matrices + metrics as $outputfile"
    save(outputfile,
        Dict(
            "T" => T,
            "Tadv" => Tadv,
            "TκH" => TκH,
            "TκVML" => TκVML,
            "TκVdeep" => TκVdeep,
            "gridmetrics" => gridmetrics,
            "indices" => indices,
            "κH" => κH,
            "κVML" => κVML,
            "κVdeep" => κVdeep,
            "model" => model,
            "experiment" => experiment,
            "member" => member,
            "time_window" => time_window,
            "note" => "built from averaging mass transport, $str",
            "OceanTransportMatrixBuilder" => "v$(pkgversion(OceanTransportMatrixBuilder))",
            "Date" => Dates.now()
        )
    )


    @info "Computing timescales (ideal mean age + reemergence time)"

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

    @info "Solving for ideal mean age"
    Γin = A \ sΓin
    Γinyr = ustrip.(yr, Γin .* s)
    Γinyr3D = OceanTransportMatrixBuilder.as3D(Γinyr, wet3D)

    # global mean of the mean age should be order 1000yr
    @info "global average ideal mean age is $((v' * Γinyr) / sum(v)) (should be order 1000yr)"

    # Mean reemergence time is governed by
    # 	-∂Γ↑/∂t + (T̃ + M) Γ↑ = 1    (note the tilde on T̃ = V⁻¹ Tᵀ V, the adjoint of T)
    # For computational efficiency, reuse factors of A = T + M:
    #   Γ↑ = V⁻¹ A⁻ᵀ V 1
    V = sparse(Diagonal(v))
    V⁻¹ = sparse(Diagonal(1 ./ v))
    @info "Solving for reemergence time"
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
    # Turn Γout into a YAXArray by rebuilding from volcello
    Γoutyr_YAXArray = rebuild(volcello_ds["volcello"];
        data = Γoutyr3D,
        dims = dims(volcello_ds["volcello"]),
        metadata = Dict(
            "origin" => "mean_reemergence_time computed from $model $experiment $member $(time_window)",
            "units" => "yr",
        )
    )

    arrays = Dict(:ideal_age => Γinyr_YAXArray, :reemergence_time => Γoutyr_YAXArray, :lat => volcello_ds.lat, :lon => volcello_ds.lon)
    Γ_ds = Dataset(; volcello_ds.properties, arrays...)

    # Save to netCDF file
    outputfile = joinpath(inputdir, "circulation_timescales_$str.nc")
    @info "Saving as netCDF file:\n  $(outputfile)"
    savedataset(Γ_ds, path = outputfile, driver = :netcdf, overwrite = true)


end