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
include("GentMcWilliams.jl")
CSIRO_member = CMIP6member2CSIROmember(member)
Tilodatainputdir = "/g/data/xv83/TMIP/data/$model/$CSIRO_member"
ψᵢGM_ds = open_dataset(joinpath(Tilodatainputdir, "month_tx_trans_gm_1990s.nc"))
ψⱼGM_ds = open_dataset(joinpath(Tilodatainputdir, "month_ty_trans_gm_1990s.nc"))
ψᵢsubmeso_ds = open_dataset(joinpath(Tilodatainputdir, "month_tx_trans_submeso_1990s.nc"))
ψⱼsubmeso_ds = open_dataset(joinpath(Tilodatainputdir, "month_ty_trans_submeso_1990s.nc"))

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



yearly_T = Dict{Int64, SparseMatrixCSC{Float64, Int64}}()
yearly_Tadv = Dict{Int64, SparseMatrixCSC{Float64, Int64}}()
yearly_TκH = Dict{Int64, SparseMatrixCSC{Float64, Int64}}()
yearly_TκVML = Dict{Int64, SparseMatrixCSC{Float64, Int64}}()
yearly_TκVdeep = Dict{Int64, SparseMatrixCSC{Float64, Int64}}()
yearly_weights = Dict{Int64, Float64}()

years = 1990:1999
months = 1:12
# years = 1990:1991
# months = 1:2

for year in years

    println("year $year")

    monthly_T = Dict{Int64, SparseMatrixCSC{Float64, Int64}}()
    monthly_Tadv = Dict{Int64, SparseMatrixCSC{Float64, Int64}}()
    monthly_TκH = Dict{Int64, SparseMatrixCSC{Float64, Int64}}()
    monthly_TκVML = Dict{Int64, SparseMatrixCSC{Float64, Int64}}()
    monthly_TκVdeep = Dict{Int64, SparseMatrixCSC{Float64, Int64}}()
    monthly_weights = Dict{Int64, Float64}()

    for month in months

        println("  month $month")

        # Load variables in memory
        mlotst = readcubedata(mlotst_ds.mlotst[Ti=Date(year, month, 10)..Date(year, month, 20)])
        umo = readcubedata(umo_ds.umo[Ti=Date(year, month, 10)..Date(year, month, 20)])
        vmo = readcubedata(vmo_ds.vmo[Ti=Date(year, month, 10)..Date(year, month, 20)])

        ψᵢGM = readcubedata(ψᵢGM_ds.tx_trans_gm[Ti=Date(year, month, 10)..Date(year, month, 20)])
        ψⱼGM = readcubedata(ψⱼGM_ds.ty_trans_gm[Ti=Date(year, month, 10)..Date(year, month, 20)])
        ψᵢsubmeso = readcubedata(ψᵢsubmeso_ds.tx_trans_submeso[Ti=Date(year, month, 10)..Date(year, month, 20)])
        ψⱼsubmeso = readcubedata(ψⱼsubmeso_ds.ty_trans_submeso[Ti=Date(year, month, 10)..Date(year, month, 20)])

        # Replace missing values and convert to arrays
        # I think latest YAXArrays converts _FillValues to missing
        ψᵢGM = replace(ψᵢGM |> Array, missing => 0) .|> Float64
        ψⱼGM = replace(ψⱼGM |> Array, missing => 0) .|> Float64
        ψᵢsubmeso = replace(ψᵢsubmeso |> Array, missing => 0) .|> Float64
        ψⱼsubmeso = replace(ψⱼsubmeso |> Array, missing => 0) .|> Float64

        # Also remove missings in umo and vmo
        umo = replace(umo, missing => 0)
        vmo = replace(vmo, missing => 0)

        # Take the vertical diff of zonal/meridional transport diagnostics to get their mass transport
        (nx, ny, _) = size(ψᵢGM)
        ϕᵢGM = diff([fill(0.0, nx, ny, 1);;; ψᵢGM |> Array], dims=3)
        ϕⱼGM = diff([fill(0.0, nx, ny, 1);;; ψⱼGM |> Array], dims=3)
        ϕᵢsubmeso = diff([fill(0.0, nx, ny, 1);;; ψᵢsubmeso |> Array], dims=3)
        ϕⱼsubmeso = diff([fill(0.0, nx, ny, 1);;; ψⱼsubmeso |> Array], dims=3)

        # TODO fix incompatible dimensions betwewen umo and ϕᵢGM/ϕᵢsubmeso Dim{:i} and Dim{:xu_ocean}
        ϕ = let umo = umo + ϕᵢGM + ϕᵢsubmeso, vmo = vmo + ϕⱼGM + ϕⱼsubmeso
            facefluxesfrommasstransport(; umo, vmo, gridmetrics, indices)
        end

        (; T, Tadv, TκH, TκVML, TκVdeep) = transportmatrix(; ϕ, mlotst, gridmetrics, indices, ρ, κH, κVML, κVdeep)

        monthly_T[month] = T
        monthly_Tadv[month] = Tadv
        monthly_TκH[month] = TκH
        monthly_TκVML[month] = TκVML
        monthly_TκVdeep[month] = TκVdeep
        monthly_weights[month] = daysinmonth(Date(year, month, 15))

    end

    @show monthly_weights
    yearly_weights[year] = sum(w for w in values(monthly_weights))

    @info "summing monthly matrices for year $year"
    @time yearly_T[year], yearly_Tadv[year], yearly_TκH[year], yearly_TκVML[year], yearly_TκVdeep[year] = let
        T = monthly_weights[first(months)] * monthly_T[first(months)]
        Tadv = monthly_weights[first(months)] * monthly_Tadv[first(months)]
        TκH = monthly_weights[first(months)] * monthly_TκH[first(months)]
        TκVML = monthly_weights[first(months)] * monthly_TκVML[first(months)]
        TκVdeep = monthly_weights[first(months)] * monthly_TκVdeep[first(months)]
        for month in months[2:end]
            T += monthly_weights[month] * monthly_T[month]
            Tadv += monthly_weights[month] * monthly_Tadv[month]
            TκH += monthly_weights[month] * monthly_TκH[month]
            TκVML += monthly_weights[month] * monthly_TκVML[month]
            TκVdeep += monthly_weights[month] * monthly_TκVdeep[month]
        end
        T /= yearly_weights[year]
        Tadv /= yearly_weights[year]
        TκH /= yearly_weights[year]
        TκVML /= yearly_weights[year]
        TκVdeep /= yearly_weights[year]
        T, Tadv, TκH, TκVML, TκVdeep
    end

end

@show yearly_weights

@info "summing yearly matrices for years $years"
total_weights = sum(w for w in values(yearly_weights))
@time T, Tadv, TκH, TκVML, TκVdeep = let
    T = yearly_weights[first(years)] * yearly_T[first(years)]
    Tadv = yearly_weights[first(years)] * yearly_Tadv[first(years)]
    TκH = yearly_weights[first(years)] * yearly_TκH[first(years)]
    TκVML = yearly_weights[first(years)] * yearly_TκVML[first(years)]
    TκVdeep = yearly_weights[first(years)] * yearly_TκVdeep[first(years)]
    for year in years[2:end]
        T += yearly_weights[year] * yearly_T[year]
        Tadv += yearly_weights[year] * yearly_Tadv[year]
        TκH += yearly_weights[year] * yearly_TκH[year]
        TκVML += yearly_weights[year] * yearly_TκVML[year]
        TκVdeep += yearly_weights[year] * yearly_TκVdeep[year]
    end
    T /= total_weights
    Tadv /= total_weights
    TκH /= total_weights
    TκVML /= total_weights
    TκVdeep /= total_weights
    T, Tadv, TκH, TκVML, TκVdeep
end

WIP

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
outputfile = joinpath(CMIP6outputdir, "ideal_mean_age_$str.nc")
@info "Saving ideal mean age as netCDF file:\n  $(outputfile)"
savedataset(Γin_ds, path = outputfile, driver = :netcdf, overwrite = true)

# Save Γoutyr3D to netCDF file
outputfile = joinpath(CMIP6outputdir, "mean_reemergence_time_$str.nc")
@info "Saving mean reemergence time as netCDF file:\n  $(outputfile)"
savedataset(Γout_ds, path = outputfile, driver = :netcdf, overwrite = true)



