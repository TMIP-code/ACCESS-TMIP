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
using CairoMakie
using GeoMakie
using Interpolations
using OceanBasins
using NaNStatistics
using FileIO
using Makie.StructArrays

# Load functions for GM terms
include("GentMcWilliams.jl")
include("plotting_functions.jl")

model = "ACCESS-ESM1-5"
member = "r1i1p1f1"
CMIP_version = "CMIP5"
experiment = "historical"
time_window = "Jan1990-Dec1999"

# Gadi directory for input files
inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)"

# Load umo, vmo, mlotst, volcello, and areacello
volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))

# Load variables in memory
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

strs = ("resolved_GM_submeso", "frommonthlymatrices")
# strs = ("resolved_GM_submeso_maxMLD", "frommonthlymatrices")

# Load ideal mean age and reemergence time
Ts = [load(joinpath(inputdir, "transportmatrix_$k.jld2")) for k in strs]

fig = Figure(size=(1000, 2000))

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
nskip = 1

myscale(x) = Makie.pseudolog10(1e7x)
powersoften = -6:-2
ticks = [-reverse(exp10.(powersoften)); 0; exp10.(powersoften)]
signedstr(x) = x > 0 ? "+$x" : "−$(-x)"
ticklabels = [rich("10", superscript(signedstr(i))) for i in powersoften]
ticklabels = [[rich("−", x) for x in reverse(ticklabels)]; rich("0"); [rich("+", x) for x in ticklabels]]
xticks = yticks = (myscale.(ticks), ticklabels)
ax = Axis(ga[1,1]; title = "T", xlabel = strs[1], ylabel = strs[2], xticks, yticks, aspect = 1)
vx, vy = makesparseentriescomparable(Ts[1]["T"], Ts[2]["T"], nskip)
# scatter!(ax, myscale.(vx), myscale.(vy), color = :black, markersize = 3)
points = StructArray{Point2f}((myscale.(vx), myscale.(vy)))
colormap = cgrad([:white; collect(cgrad(:managua))])
ds = datashader!(ax, points; colormap, async = false)
ablines!(ax, 0, 1, color = (:black, 0.1), linewidth = 10)
vlines!(ax, 0, color = (:black, 0.1), linewidth = 10)
hlines!(ax, 0, color = (:black, 0.1), linewidth = 10)

subplots = [
    "Tadv"    "TκH"
    "TκVdeep" "TκVML"
]

for I in eachindex(IndexCartesian(), subplots)
    irow, icol  = Tuple(I)
    Tstr = subplots[irow, icol]
    subax = Axis(gb[irow, icol]; title = Tstr, xlabel = strs[1], ylabel = strs[2], xticks, yticks, aspect = 1)
    local vx, vy = makesparseentriescomparable(Ts[1][Tstr], Ts[2][Tstr], nskip)
    # scatter!(subax, vx, vy, color = :black, markersize = 3)
    # scatter!(subax, myscale2.(vx), myscale2.(vy), color = :black, markersize = 3)
    local points = StructArray{Point2f}((myscale.(vx), myscale.(vy)))
    datashader!(subax, points; colormap, async = false)
    ablines!(subax, 0, 1, color = (:black, 0.1), linewidth = 10)
    vlines!(subax, 0, color = (:black, 0.1), linewidth = 10)
    hlines!(subax, 0, color = (:black, 0.1), linewidth = 10)
    linkaxes!(subax, ax)
    hidexdecorations!(subax,
        label = irow < 2, ticklabels = irow < 2,
        ticks = irow < 2, grid = false)
    hideydecorations!(subax,
        label = icol > 1, ticklabels = icol > 1,
        ticks = icol > 1, grid = false)
end

cb = Colorbar(fig[end+1, :], ds, label="Density of transport matrix coefficients", vertical = false, flipaxis = false)
cb.width = Relative(2/3)

# save plot
outputfile = joinpath(inputdir, "matrices_$(strs[1])_vs_$(strs[2]).png")
@info "Saving matrices comparsion as image file:\n  $(outputfile)"
save(outputfile, fig)

