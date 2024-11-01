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

# Making monthly matrices for every MONTH from 1990s and averaging them into a single matrix, and computing the age

model = "ACCESS-ESM1-5"
member = "r1i1p1f1"
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

YEARS = 1990:1999
MONTHS = 1:12

for YEAR in YEARS

    println("YEAR $YEAR")

    monthly_T = Dict{Int64, SparseMatrixCSC{Float64, Int64}}()
    monthly_Tadv = Dict{Int64, SparseMatrixCSC{Float64, Int64}}()
    monthly_TκH = Dict{Int64, SparseMatrixCSC{Float64, Int64}}()
    monthly_TκVML = Dict{Int64, SparseMatrixCSC{Float64, Int64}}()
    monthly_TκVdeep = Dict{Int64, SparseMatrixCSC{Float64, Int64}}()
    monthly_weights = Dict{Int64, Float64}()

    for MONTH in MONTHS

        println("  MONTH $MONTH")

        # Load variables in memory
        mlotst = readcubedata(mlotst_ds.mlotst[Ti=Where(x -> month(x) == MONTH && year(x) == YEAR)])
        umo = readcubedata(umo_ds.umo[Ti=Where(x -> month(x) == MONTH && year(x) == YEAR)])
        vmo = readcubedata(vmo_ds.vmo[Ti=Where(x -> month(x) == MONTH && year(x) == YEAR)])

        ψᵢGM = readcubedata(ψᵢGM_ds.tx_trans_gm[Ti=Where(x -> month(x) == MONTH && year(x) == YEAR)])
        ψⱼGM = readcubedata(ψⱼGM_ds.ty_trans_gm[Ti=Where(x -> month(x) == MONTH && year(x) == YEAR)])
        ψᵢsubmeso = readcubedata(ψᵢsubmeso_ds.tx_trans_submeso[Ti=Where(x -> month(x) == MONTH && year(x) == YEAR)])
        ψⱼsubmeso = readcubedata(ψⱼsubmeso_ds.ty_trans_submeso[Ti=Where(x -> month(x) == MONTH && year(x) == YEAR)])

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

        monthly_T[MONTH] = T
        monthly_Tadv[MONTH] = Tadv
        monthly_TκH[MONTH] = TκH
        monthly_TκVML[MONTH] = TκVML
        monthly_TκVdeep[MONTH] = TκVdeep
        monthly_weights[MONTH] = daysinmonth(Date(YEAR, MONTH, 15))

    end

    @show monthly_weights
    yearly_weights[YEAR] = sum(w for w in values(monthly_weights))

    @info "summing monthly matrices for YEAR $YEAR"
    @time yearly_T[YEAR], yearly_Tadv[YEAR], yearly_TκH[YEAR], yearly_TκVML[YEAR], yearly_TκVdeep[YEAR] = let
        T = monthly_weights[first(MONTHS)] * monthly_T[first(MONTHS)]
        Tadv = monthly_weights[first(MONTHS)] * monthly_Tadv[first(MONTHS)]
        TκH = monthly_weights[first(MONTHS)] * monthly_TκH[first(MONTHS)]
        TκVML = monthly_weights[first(MONTHS)] * monthly_TκVML[first(MONTHS)]
        TκVdeep = monthly_weights[first(MONTHS)] * monthly_TκVdeep[first(MONTHS)]
        for MONTH in MONTHS[2:end]
            T += monthly_weights[MONTH] * monthly_T[MONTH]
            Tadv += monthly_weights[MONTH] * monthly_Tadv[MONTH]
            TκH += monthly_weights[MONTH] * monthly_TκH[MONTH]
            TκVML += monthly_weights[MONTH] * monthly_TκVML[MONTH]
            TκVdeep += monthly_weights[MONTH] * monthly_TκVdeep[MONTH]
        end
        T /= yearly_weights[YEAR]
        Tadv /= yearly_weights[YEAR]
        TκH /= yearly_weights[YEAR]
        TκVML /= yearly_weights[YEAR]
        TκVdeep /= yearly_weights[YEAR]
        T, Tadv, TκH, TκVML, TκVdeep
    end

end

@show yearly_weights

@info "summing yearly matrices for YEARS $YEARS"
total_weights = sum(w for w in values(yearly_weights))
@time T, Tadv, TκH, TκVML, TκVdeep = let
    T = yearly_weights[first(YEARS)] * yearly_T[first(YEARS)]
    Tadv = yearly_weights[first(YEARS)] * yearly_Tadv[first(YEARS)]
    TκH = yearly_weights[first(YEARS)] * yearly_TκH[first(YEARS)]
    TκVML = yearly_weights[first(YEARS)] * yearly_TκVML[first(YEARS)]
    TκVdeep = yearly_weights[first(YEARS)] * yearly_TκVdeep[first(YEARS)]
    for YEAR in YEARS[2:end]
        T += yearly_weights[YEAR] * yearly_T[YEAR]
        Tadv += yearly_weights[YEAR] * yearly_Tadv[YEAR]
        TκH += yearly_weights[YEAR] * yearly_TκH[YEAR]
        TκVML += yearly_weights[YEAR] * yearly_TκVML[YEAR]
        TκVdeep += yearly_weights[YEAR] * yearly_TκVdeep[YEAR]
    end
    T /= total_weights
    Tadv /= total_weights
    TκH /= total_weights
    TκVML /= total_weights
    TκVdeep /= total_weights
    T, Tadv, TκH, TκVML, TκVdeep
end

# Save matrices
outputfile = joinpath(CMIP6outputdir, "transportmatrix_frommonthlymatrices.jld2")
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
        "note" => "built from averaging monthly matrices",
        "OceanTransportMatrixBuilder" => "v$(pkgversion(OceanTransportMatrixBuilder))",
    )
)
