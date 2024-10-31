using Pkg
Pkg.activate(".")
Pkg.instantiate()

using OceanTransportMatrixBuilder
using NetCDF
using YAXArrays
using DataFrames
using DimensionalData
# using SparseArrays
# using LinearAlgebra
using Unitful
using Unitful: s, yr, m, kg, K, °C
using CairoMakie
using GeoMakie
using PairPlots
using Interpolations
using OceanBasins
using NaNStatistics
using Format
using GibbsSeaWater
using Statistics

include("plotting_functions.jl")

# 1. Load data from Matt's ACCESS1-3 = the data used for the GM terms in previous papers.
GMmasstransport_fromMattChamberlain = let
    # Extra stuff for GM terms
    MattChamberlaindirectory = "/g/data/v19/mtc599/matrix/access1.3"
    GMtermfilepath = joinpath(MattChamberlaindirectory, "historical/grid1990/transport.nc")
    # Below is required because MC's directory names
    Sv = 1e6m^3/s * 1000kg/m^3
    @info """
    Loading GM velocities from Matt Chamberlain's ACCESS1-3 data
    from file: $GMtermfilepath
    """
    uGM = ncread(GMtermfilepath, "tx_trans_gm") |> Array .|> Float64
    vGM = ncread(GMtermfilepath, "ty_trans_gm") |> Array .|> Float64
    # convert from Sverdrups to kg/s
    uGM = ustrip.(kg/s, uGM * Sv)
    vGM = ustrip.(kg/s, vGM * Sv)
    # the GM fields are summed vertically, so we need their diff.
    (nx, ny, _) = size(uGM)
    uGM = diff([fill(0.0, nx, ny, 1);;; uGM], dims=3)
    vGM = diff([fill(0.0, nx, ny, 1);;; vGM], dims=3)
    (; uGM, vGM)
end


# 2. Load the GM velocities that I build from CMIP output
# begin
    # model = "ACCESS-ESM1-5"
    # model = "ACCESS-CM2"
    model = "ACCESS1-3"
    experiment = "historical"
    member = "r1i1p1"
    time_window = "Jan1990-Dec1999"


    # Directory for data
    if gethostname() == "benoits-MacBook-Pro.local"
        DATADIR = "/Users/benoitpasquier/Data/TMIP/data"
    else # on Gadi. TODO: change this to the correct path
        DATADIR = "/scratch/xv83/TMIP/data"
    end

    inputdir = joinpath(DATADIR, "$model/$experiment/$member/$(time_window)")
    outputdir = joinpath("output", "plots")
    mkpath(outputdir)

    # Load datasets lazily
    volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
    areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))
    umo_ds = open_dataset(joinpath(inputdir, "umo.nc"))
    vmo_ds = open_dataset(joinpath(inputdir, "vmo.nc"))

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
    umo = readcubedata(umo_ds.umo)
    vmo = readcubedata(vmo_ds.vmo)
    umo_lon = readcubedata(umo_ds.lon)
    umo_lat = readcubedata(umo_ds.lat)
    vmo_lon = readcubedata(vmo_ds.lon)
    vmo_lat = readcubedata(vmo_ds.lat)

    # Make makegridmetrics
    gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
    (; lon_vertices, lat_vertices, v3D, zt, Z3D, thkcello) = gridmetrics


    # Load thetato and so to compute density
    @info """
    Loading thetao and so to compute density
    from files: $(joinpath(inputdir, "thetao.nc")), $(joinpath(inputdir, "so.nc"))
    """
    thetao_ds = open_dataset(joinpath(inputdir, "thetao.nc"))
    so_ds = open_dataset(joinpath(inputdir, "so.nc"))
    # Load variables in memory
    thetao = readcubedata(thetao_ds.thetao)
    @show thetao_mean = nanmean(thetao)
    if thetao_mean > 200 # If Temperature is in K (i.e., typically >200)
        thetao = ustrip.(°C, thetao * K) # convert to °C
    end
    so = readcubedata(so_ds.so)
    @show nanmean(so)
    # Convert thetao and so to density
    ct = gsw_ct_from_pt.(so, thetao)
    @show nanmean(ct)
    ρ = gsw_rho.(so, ct, Z3D)

    indices = makeindices(v3D)

    κGM = 600 # m^2/s
    maxslope = 0.01
    uGM, vGM = OceanTransportMatrixBuilder.bolus_GM_velocity(ρ, gridmetrics, indices; κGM, maxslope)
    # Turn uGM into a YAXArray by rebuilding from uo
    uGM_YAXArray = rebuild(umo;
        data = uGM,
        dims = dims(umo),
        metadata = Dict(
            "origin" => "u GM bolus velocity computed from thetao and so",
            "kGM" => κGM,
            "units" => "m/s",
        )
    )
    arrays = Dict(:uo_GM => uGM_YAXArray, :lat => umo_lat, :lon => umo_lon)
    uGM_ds = Dataset(; umo_ds.properties, arrays...)
    # Save to netCDF file
    # outputfile = joinpath(outputdir, "uo_GM.nc")
    # @info "uo_GM as netCDF file:\n  $(outputfile)"
    # savedataset(uGM_ds, path = outputfile, driver = :netcdf, overwrite = true)
    # Turn vGM into a YAXArray by rebuilding from uo
    vGM_YAXArray = rebuild(vmo;
        data = vGM,
        dims = dims(vmo),
        metadata = Dict(
            "origin" => "v GM bolus velocity computed from thetao and so",
            "kGM" => κGM,
            "units" => "m/s",
        )
    )
    arrays = Dict(:vo_GM => vGM_YAXArray, :lat => vmo_lat, :lon => vmo_lon)
    vGM_ds = Dataset(; vmo_ds.properties, arrays...)
    # Save to netCDF file
    # outputfile = joinpath(outputdir, "vo_GM.nc")
    # @info "Saving vo_GM as netCDF file:\n  $(outputfile)"
    # savedataset(vGM_ds, path = outputfile, driver = :netcdf, overwrite = true)
    umo_GM, vmo_GM = velocity2fluxes(uGM, umo_lon, umo_lat, vGM, vmo_lon, vmo_lat, gridmetrics, ρ)
# end

# 3. Load the GM terms for ACCESS-ESM1-5
@info "Loading GM velocities from ACCESS-ESM1-5 data"
ACCESSESM15directory = "/g/data/p73/archive/CMIP6/ACCESS-ESM1-5/"
GMtermfilepaths = [joinpath(ACCESSESM15directory, "HI-05/history/ocn/ocean_month.nc-199$(i)1231") for i in 0:9]
FILLVALUE = -1f20
@info "  u GM"
ΔϕGMᵢmonthly = replace(ncread(first(GMtermfilepaths), "tx_trans_gm") |> Array .|> Float64, FILLVALUE => 0.0)
ΔϕGMᵢmean = dropdims(nanmean(ΔϕGMᵢmonthly, dims=4), dims=4)
for filepath in GMtermfilepaths[2:end]
    ΔϕGMᵢmonthly .= replace(ncread(filepath, "tx_trans_gm") |> Array .|> Float64, FILLVALUE => 0.0)
    ΔϕGMᵢmean .+= dropdims(nanmean(ΔϕGMᵢmonthly, dims=4), dims=4)
end
ΔϕGMᵢmean ./= length(GMtermfilepaths)
@info "  v GM"
ΔϕGMⱼmonthly = replace(ncread(first(GMtermfilepaths), "ty_trans_gm") |> Array .|> Float64, FILLVALUE => 0.0)
ΔϕGMⱼmean = dropdims(nanmean(ΔϕGMⱼmonthly, dims=4), dims=4)
for filepath in GMtermfilepaths[2:end]
    ΔϕGMⱼmonthly .= replace(ncread(filepath, "ty_trans_gm") |> Array .|> Float64, FILLVALUE => 0.0)
    ΔϕGMⱼmean .+= dropdims(nanmean(ΔϕGMⱼmonthly, dims=4), dims=4)
end
ΔϕGMⱼmean ./= length(GMtermfilepaths)
@info "  Taking vertical diff"
(nx, ny, _) = size(ΔϕGMᵢmean)
ϕGMᵢmean = diff([fill(0.0, nx, ny, 1);;; ΔϕGMᵢmean], dims=3)
ϕGMⱼmean = diff([fill(0.0, nx, ny, 1);;; ΔϕGMⱼmean], dims=3)
GMmasstransport_fromACCESSESM15 = (; ϕGMᵢmean, ϕGMⱼmean)


# make indices
indices = makeindices(gridmetrics.v3D)
(; Lwet, wet3D) = indices

# Plot pair plots of u and v
u1 = GMmasstransport_fromMattChamberlain.uGM[wet3D] * 1e-9
u2 = umo_GM[wet3D] * 1e-9
u3 = GMmasstransport_fromACCESSESM15.ϕGMᵢmean[wet3D] * 1e-9
v1 = GMmasstransport_fromMattChamberlain.vGM[wet3D] * 1e-9
v2 = vmo_GM[wet3D] * 1e-9
v3 = GMmasstransport_fromACCESSESM15.ϕGMⱼmean[wet3D] * 1e-9

idx = findall(@. (u1 ≠ 0) & (u2 ≠ 0) & (u3 ≠ 0) & (v1 ≠ 0) & (v2 ≠ 0) & (v3 ≠ 0) & !isnan(u1) & !isnan(u2) & !isnan(u3) & !isnan(v1) & !isnan(v2) & !isnan(v3))
df = DataFrame(uMC = u1[idx], uBP = u2[idx], uESM15 = u3[idx], vMC = v1[idx], vBP = v2[idx], vESM15 = v3[idx])
# filter zeros out
lims = (low = -0.1, high = +0.1)
axis = Dict(:uMC => (; lims), :uBP => (; lims), :uESM15 => (; lims), :vMC => (; lims), :vBP => (; lims), :vESM15 => (; lims))
bins = -0.1:0.01:0.1
bins = Dict(:uMC => bins, :uBP => bins, :uESM15 => bins, :vMC => bins, :vBP => bins, :vESM15 => bins)
fig = pairplot(df; bins, axis)
# save plot
outputfile = joinpath(outputdir, "GM_masstransport_MC_vs_BP.png")
@info "Saving GM ϕ comparison as image file:\n  $(outputfile)"
save(outputfile, fig)

@info "Relative sizes of GM terms"
@show [mean(abs.(u1[idx])), mean(abs.(u2[idx])), mean(abs.(u3[idx]))]
@show [mean(abs.(v1[idx])), mean(abs.(v2[idx])), mean(abs.(v3[idx]))]














foo

inputdirfun(member) = joinpath(DATADIR, "$model/$experiment/$member/$(time_window)")

# find all members for which the inputdir contains umo.nc, vmo.nc, mlotst.nc, volcello.nc, and areacello.nc
requiredvariables = ["umo", "vmo", "mlotst", "volcello", "areacello", "agessc"]
hasrequireddata(member, variable_name) = isfile(joinpath(inputdirfun(member), "$variable_name.nc"))
hasrequireddata(member) = all(variable_name -> hasrequireddata(member, variable_name), requiredvariables)
members = readdir(joinpath(DATADIR, "$model/$experiment"))
members = [m for m in members if m ≠ ".DS_Store"]

# sort members by r, i, p[, f]
memmber_regex = CMIP_version == "CMIP6" ? r"r(\d+)i(\d+)p(\d+)f(\d+)" : r"r(\d+)i(\d+)p(\d+)"
parse_member(member) = parse.(Int, match(memmber_regex, member).captures)
members = sort(members, by = x -> parse_member(x))
dataavailability = DataFrame(
    :member => members,
    :has_it_all => hasrequireddata.(members),
    [Symbol(var) => [hasrequireddata(member, var) for member in members] for var in requiredvariables]...,
)
show(dataavailability; allrows = true)
println()


# for member in members[dataavailability.has_it_all]
# for member in [last(members)]
member = first(members)

    inputdir = inputdirfun(member)

    # Load umo, vmo, mlotst, volcello, and areacello
    umo_ds = open_dataset(joinpath(inputdir, "umo.nc"))
    vmo_ds = open_dataset(joinpath(inputdir, "vmo.nc"))
    uo_ds = open_dataset(joinpath(inputdir, "uo.nc"))
    vo_ds = open_dataset(joinpath(inputdir, "vo.nc"))
    mlotst_ds = open_dataset(joinpath(inputdir, "mlotst.nc"))
    volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
    areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))

    # Load the required variables
    # Note: mass transport is expected to live on the (u,v) coordinates of an Arakawa C-grid
    umo = umo_ds.umo # mass transport across "east" cell face
    vmo = vmo_ds.vmo # mass transport across "north" cell face
    uo = uo_ds.uo # velocities
    vo = vo_ds.vo # velocities
    mlotst = mlotst_ds.mlotst # MLD
    areacello = areacello_ds.areacello # top area of cells
    volcello = volcello_ds.volcello # volume of cells
    lon = volcello_ds.lon # longitude of cell centers
    lat = volcello_ds.lat # latitude of cell centers
    lev = volcello_ds.lev # depth of cell centers
    if :lon_verticies in propertynames(volcello_ds)
        lon_vertices = volcello_ds.lon_verticies # cell vertices
        lat_vertices = volcello_ds.lat_verticies # cell vertices
    else
        lon_vertices = volcello_ds.lon_vertices # cell vertices
        lat_vertices = volcello_ds.lat_vertices # cell vertices
    end
    # "verticies" name is an xmip bug: https://github.com/jbusecke/xMIP/issues/369
    # (my local data was preprocessed with xmip)

    # Make the required data from grid geometry
    gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)

    # Make the indices for going back and forth between 3D and 1D
    indices = makeindices(gridmetrics.v3D)

    # Some parameter values
    ρ = 1035.0    # density (kg/m^3)

    # Make arrays of the flux on each face for each grid cell
    ϕ = facefluxesfrommasstransport(; umo, vmo, gridmetrics, indices)

    # Make arrays of the flux from velocities as well
    uo_lon = modify(Array, uo_ds.lon)
    uo_lat = modify(Array, uo_ds.lat)
    vo_lon = modify(Array, vo_ds.lon)
    vo_lat = modify(Array, vo_ds.lat)
    ϕbis = facefluxesfromvelocities(; uo, uo_lon, uo_lat, vo, vo_lon, vo_lat, gridmetrics, ρ)

    # unpack model grid
    (; lon, lat, zt, v3D,) = gridmetrics
    lev = zt
    # unpack indices
    (; wet3D, N) = indices



    # basins
    basin_keys = (:ATL, :PAC, :IND)
    basin_strs = ("Atlantic", "Pacific", "Indian")
    basin_functions = (isatlantic, ispacific, isindian)
    basin_values = (reshape(f(lat[:], lon[:], OCEANS), size(lat)) for f in basin_functions)
    basins = (; (basin_keys .=> basin_values)...)
    basin_latlims_values = [clamp.((-5, +5) .+ extrema(lat[.!isnan.(v3D[:,:,1]) .& basin[:,:,1]]), -80, 80) for basin in basins]
    basin_latlims = (; (basin_keys .=> basin_latlims_values)...)




    # Plot meriodional overturning circulation for each basin

    levels = -20:20
    colormap = cgrad(:curl, length(levels) + 1; categorical=true)
    extendlow = colormap[1]
    extendhigh = colormap[end]
    colormap = cgrad(colormap[2:end-1]; categorical=true)

    fig = Figure(size = (1200, 600), fontsize = 18)
    axs = Array{Any, 2}(undef, (2, 3))
    contours = Array{Any, 2}(undef, (2, 3))
    for (icol, (basin_key, basin)) in enumerate(pairs(basins))

        for (irow, x3D) in enumerate((ϕ.north, ϕbis.north))

            x2D = dropdims(reverse(nancumsum(reverse(nansum(basin .* x3D, dims = 1), dims=3), dims = 3), dims=3), dims = 1) # kg/s
            x2Dmask = zonalaverage(1, gridmetrics; mask = basin) .> 0
            x2D[.!x2Dmask] .= NaN

            # convert to Sv
            x2D = x2D / 1e6 / ρ # Sv

            local ax = Axis(fig[irow, icol],
                backgroundcolor=:lightgray,
                xgridvisible=false, ygridvisible=false,
                ylabel = "depth (m)")

            X = dropdims(maximum(lat, dims=1), dims=1)
            Y = zt
            Z = x2D
            co = contourf!(ax, X, Y, Z;
                levels,
                colormap,
                nan_color = :lightgray,
                extendlow,
                extendhigh,
            )
            translate!(co, 0, 0, -100)
            contours[irow, icol] = co

            xlim = basin_latlims[basin_key]
            # basin2 = LONGTEXT[basin]

            ax.yticks = (ztick, zticklabel)
            xticks = -90:30:90
            ax.xticks = (xticks, latticklabel.(xticks))
            ylims!(ax, zlim)
            # xlims!(ax, (-90, 90))
            xlims!(ax, xlim)


            hidexdecorations!(ax,
                label = irow < size(axs, 1), ticklabels = irow < size(axs, 1),
                ticks = irow < size(axs, 1), grid = false)
            hideydecorations!(ax,
                label = icol > 1, ticklabels = icol > 1,
                ticks = icol > 1, grid = false)


            axs[irow, icol] = ax
        end

        Label(fig[0, icol], basin_strs[icol], tellwidth = false)

    end

    cb = Colorbar(fig[1:size(axs, 1), 4], contours[1, 1];
        vertical = true, flipaxis = true,
        # ticks = (, cbarticklabelformat.(levels)),
        tickformat = x -> map(t -> replace(format("{:+d}", t), "-" => "−"), x),
        label = "MOC (Sv)",
        )
    cb.height = Relative(0.666)

    rowgap!(fig.layout, 15)
    # rowgap!(fig.layout, 4, 15)
    # # # rowgap!(fig.layout, 5, 10)
    colgap!(fig.layout, 15)
    title = "$model $experiment $member $(time_window) MOC"
    Label(fig[-1, 1:3]; text = title, fontsize = 20, tellwidth = false)
    Label(fig[1, 0]; text = "from mass transports", rotation = π/2, tellheight = false)
    Label(fig[2, 0]; text = "from velocities", rotation = π/2, tellheight = false)
    fig
    # save plot
    outputfile = joinpath(inputdir, "MOC.png")
    @info "Saving MOC as image file:\n  $(outputfile)"
    save(outputfile, fig)






# end


