# using Pkg
# Pkg.activate(".")
# Pkg.instantiate()

# using OceanTransportMatrixBuilder
# using NetCDF
# using YAXArrays
# using DataFrames
# using DimensionalData
# using SparseArrays
# using LinearAlgebra
# using Unitful
# using Unitful: s, yr
# using NaNStatistics
# using Format
# using CairoMakie
# using GeoMakie
# using Interpolations
# using OceanBasins
# using NaNStatistics
# using FileIO
# using Makie.StructArrays

# # Load functions for GM terms
# include("GentMcWilliams.jl")
# include("plotting_functions.jl")

# model = "ACCESS-ESM1-5"
# member = "r1i1p1f1"
# CMIP_version = "CMIP5"
# experiment = "historical"
# time_window = "Jan1990-Dec1999"

# # Gadi directory for input files
# inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)"

# # Load umo, vmo, mlotst, volcello, and areacello
# volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
# areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))

# # Load variables in memory
# areacello = readcubedata(areacello_ds.areacello)
# volcello = readcubedata(volcello_ds.volcello)
# lon = readcubedata(volcello_ds.lon)
# lat = readcubedata(volcello_ds.lat)
# lev = volcello_ds.lev

# # Identify the vertices keys (vary across CMIPs / models)
# volcello_keys = propertynames(volcello_ds)
# lon_vertices_key = volcello_keys[findfirst(x -> occursin("lon", x) & occursin("vert", x), string.(volcello_keys))]
# lat_vertices_key = volcello_keys[findfirst(x -> occursin("lat", x) & occursin("vert", x), string.(volcello_keys))]
# lon_vertices = readcubedata(getproperty(volcello_ds, lon_vertices_key))
# lat_vertices = readcubedata(getproperty(volcello_ds, lat_vertices_key))

# strs = ("resolved_GM_submeso", "frommonthlymatrices")

# # Load ideal mean age and reemergence time
# Ts = [load(joinpath(inputdir, "transportmatrix_$k.jld2")) for k in strs]

fig = Figure(size=(1000, 1000))

ga = fig[1, 1] = GridLayout()
gb = fig[2, 1] = GridLayout()

function makesparseentriescomparable(Tx, Ty, nskip)
    ix, jx, vx = findnz(Tx)
    Tx0 = sparse(ix, jx, 1e-30, size(Tx)...)
    iy, jy, vy = findnz(Ty)
    Ty0 = sparse(iy, jy, 1e-30, size(Ty)...)
    Tx = Tx + Ty0
    Ty = Ty + Tx0
    ix, jx, vx = findnz(Tx)
    iy, jy, vy = findnz(Ty)
    @assert isequal(ix, iy)
    @assert isequal(jx, jy)
    return vx[1:nskip:end], vy[1:nskip:end]
end
nskip = 19
ax = Axis(ga[1,1], title = "T")
vx, vy = makesparseentriescomparable(Ts[1]["T"], Ts[2]["T"], nskip)
ablines!(ax, 0, 1, color = :red)
scatter!(ax, vx, vy, color = :black, markersize = 3)
# points = StructArray{Point2f}((vx, vy))
# datashader!(ax, points, colormap=[:white, :black])

axes = [
    "Tadv"    "TκH"
    "TκVdeep" "TκVML"
]

for I in eachindex(IndexCartesian(), axes)
    irow, icol  = Tuple(I)
    Tstr = axes[irow, icol]
    local ax = Axis(gb[irow, icol], title = Tstr)
    local vx, vy = makesparseentriescomparable(Ts[1][Tstr], Ts[2][Tstr], nskip)
    ablines!(ax, 0, 1, color = :red)
    scatter!(ax, vx, vy, color = :black, markersize = 3)
end

# save plot
outputfile = joinpath(inputdir, "matrices_monthlycomparison.png")
@info "Saving matrices comparsion as image file:\n  $(outputfile)"
save(outputfile, fig)

