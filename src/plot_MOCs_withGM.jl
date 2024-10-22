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
using Unitful: s, yr, m, kg
using CairoMakie
using GeoMakie
using Interpolations
using OceanBasins
using NaNStatistics
using Format

include("plotting_functions.jl")

model = "ACCESS-ESM1-5"
# model = "ACCESS-CM2"
# model = "ACCESS1-3"

CMIP_version = model == "ACCESS1-3" ? "CMIP5" : "CMIP6"

experiment = "historical"
# experiment = "piControl"

time_window = "Jan1990-Dec1999"
# time_window = "Jan1071-Dec1100" # <- last 30 years of ACCESS-ESM1-5 piControl
# time_window = "Jan1420-Dec1449" # <- last 30 years of ACCESS-CM2 piControl


# Directory for data
if gethostname() == "benoits-MacBook-Pro.local"
    DATADIR = "/Users/benoitpasquier/Data/TMIP/data"
else # on Gadi. TODO: change this to the correct path
    DATADIR = "/scratch/xv83/TMIP/data"
end

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

# tx_trans_gm is located at the bottom of the east face of the cell. It's an "overturning diagnostic"


# for member in members[dataavailability.has_it_all]
# for member in [last(members)]
# member = first(members)
member = "r5i1p1f1" # <- presumably the ACCESS-ESM1-5 member which has GM transport available ("HI-05")
# member = "r1i1p1" # <- presumably the ACCESS1-3 member that Matt Chamberlain used for GM transport

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
    umo = readcubedata(umo_ds.umo) # mass transport across "east" cell face
    vmo = readcubedata(vmo_ds.vmo) # mass transport across "north" cell face
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




    # if model == "ACCESS-ESM1-5"

    #     # Load the GM terms for ACCESS-ESM1-5
    #     @info "Loading GM velocities from ACCESS-ESM1-5 data"
    #     ACCESSESM15directory = "/g/data/p73/archive/CMIP6/ACCESS-ESM1-5/"
    #     GMtermfilepaths = [joinpath(ACCESSESM15directory, "HI-05/history/ocn/ocean_month.nc-199$(i)1231") for i in 0:9]
    #     FILLVALUE = -1f20
    #     @info "  u GM"
    #     ψGMᵢmonthly = replace(ncread(first(GMtermfilepaths), "tx_trans_gm") |> Array .|> Float64, FILLVALUE => 0.0)
    #     ψGMᵢmean = dropdims(nanmean(iψGMᵢmonthly, dims=4), dims=4)
    #     for filepath in GMtermfilepaths[2:end]
    #         ψGMᵢmonthly .= replace(ncread(filepath, "tx_trans_gm") |> Array .|> Float64, FILLVALUE => 0.0)
    #         ψGMᵢmean .+= dropdims(nanmean(ψGMᵢmonthly, dims=4), dims=4)
    #     end
    #     ψGMᵢmean ./= length(GMtermfilepaths)
    #     @info "  v GM"
    #     ψGMⱼmonthly = replace(ncread(first(GMtermfilepaths), "ty_trans_gm") |> Array .|> Float64, FILLVALUE => 0.0)
    #     ψGMⱼmean = dropdims(nanmean(iψGMⱼmonthly, dims=4), dims=4)
    #     for filepath in GMtermfilepaths[2:end]
    #         ψGMⱼmonthly .= replace(ncread(filepath, "ty_trans_gm") |> Array .|> Float64, FILLVALUE => 0.0)
    #         ψGMⱼmean .+= dropdims(nanmean(ψGMⱼmonthly, dims=4), dims=4)
    #     end
    #     ψGMⱼmean ./= length(GMtermfilepaths)
    #     @info "  Taking vertical diff"
    #     (nx, ny, _) = size(iψGMᵢmean)
    #     ϕGMᵢmean = diff([fill(0.0, nx, ny, 1);;; ψGMᵢmean], dims=3)
    #     ϕGMⱼmean = diff([fill(0.0, nx, ny, 1);;; ψGMⱼmean], dims=3)
    #     GMmasstransport_fromACCESSESM15 = (; ϕGMᵢmean, ϕGMⱼmean)
    #     @info "  u submeso"
    #     ψsubmesoᵢmonthly = replace(ncread(first(GMtermfilepaths), "tx_trans_submeso") |> Array .|> Float64, FILLVALUE => 0.0)
    #     ψsubmesoᵢmean = dropdims(nanmean(ψsubmesoᵢmonthly, dims=4), dims=4)
    #     for filepath in GMtermfilepaths[2:end]
    #         ψsubmesoᵢmonthly .= replace(ncread(filepath, "tx_trans_submeso") |> Array .|> Float64, FILLVALUE => 0.0)
    #         ψsubmesoᵢmean .+= dropdims(nanmean(ψsubmesoᵢmonthly, dims=4), dims=4)
    #     end
    #     ψsubmesoᵢmean ./= length(GMtermfilepaths)
    #     @info "  v submeso"
    #     ψsubmesoⱼmonthly = replace(ncread(first(GMtermfilepaths), "ty_trans_submeso") |> Array .|> Float64, FILLVALUE => 0.0)
    #     ψsubmesoⱼmean = dropdims(nanmean(ψsubmesoⱼmonthly, dims=4), dims=4)
    #     for filepath in GMtermfilepaths[2:end]
    #         ψsubmesoⱼmonthly .= replace(ncread(filepath, "ty_trans_submeso") |> Array .|> Float64, FILLVALUE => 0.0)
    #         ψsubmesoⱼmean .+= dropdims(nanmean(ψsubmesoⱼmonthly, dims=4), dims=4)
    #     end
    #     ψsubmesoⱼmean ./= length(GMtermfilepaths)
    #     @info "  Taking vertical diff"
    #     (nx, ny, _) = size(iψsubmesoᵢmean)
    #     ϕsubmesoᵢmean = diff([fill(0.0, nx, ny, 1);;; ψsubmesoᵢmean], dims=3)
    #     ϕsubmesoⱼmean = diff([fill(0.0, nx, ny, 1);;; ψsubmesoⱼmean], dims=3)
    #     submesomasstransport_fromACCESSESM15 = (; ϕsubmesoᵢmean, ϕsubmesoⱼmean)

    # elseif model == "ACCESS1-3"

    #     GMmasstransport_fromMattChamberlain = begin
    #         # Extra stuff for GM terms
    #         MattChamberlaindirectory = "/g/data/v19/mtc599/matrix/access1.3"
    #         GMtermfilepath = joinpath(MattChamberlaindirectory, "historical/grid1990/transport.nc")
    #         # Below is required because MC's directory names
    #         @info """
    #         Loading GM velocities from Matt Chamberlain's ACCESS1-3 data
    #         from file: $GMtermfilepath
    #         """
    #         ϕGMᵢmean = ncread(GMtermfilepath, "tx_trans_gm") |> Array .|> Float64
    #         ϕGMⱼmean = ncread(GMtermfilepath, "ty_trans_gm") |> Array .|> Float64
    #         # Replace fill values with 0
    #         Sv = 1e6m^3/s * 1000kg/m^3
    #         ϕGMᵢmean = replace(ϕGMᵢmean, -1.0e34 => 0)
    #         ϕGMⱼmean = replace(ϕGMⱼmean, -1.0e34 => 0)
    #         # convert from Sverdrups to kg/s
    #         ϕGMᵢmean = ustrip.(kg/s, ϕGMᵢmean * Sv)
    #         ϕGMⱼmean = ustrip.(kg/s, ϕGMⱼmean * Sv)
    #         # the GM fields are summed vertically, so we need their diff.
    #         (nx, ny, _) = size(ϕGMᵢmean)
    #         ϕGMᵢmean = diff([fill(0.0, nx, ny, 1);;; ϕGMᵢmean], dims=3)
    #         ϕGMⱼmean = diff([fill(0.0, nx, ny, 1);;; ϕGMⱼmean], dims=3)
    #         (; ϕGMᵢmean, ϕGMⱼmean)
    #     end

    # end


    # Compute the fluxes without the GM terms
    ϕ = facefluxesfrommasstransport(; umo, vmo)
    # Compute the fluxes with the GM terms
    umo = umo .+ ϕGMᵢmean
    vmo = vmo .+ ϕGMⱼmean
    ϕbis = facefluxesfrommasstransport(; umo, vmo)
    # Compute the fluxes with the submseo terms
    umo = umo .+ ϕsubmesoᵢmean
    vmo = vmo .+ ϕsubmesoⱼmean
    ϕter = facefluxesfrommasstransport(; umo, vmo)


    # Plot meriodional overturning circulation for each basin

    levels = -50:5:50
    colormap = cgrad(:RdBu, length(levels) + 1; categorical=true, rev=true)
    extendlow = colormap[1]
    extendhigh = colormap[end]
    colormap = cgrad(colormap[2:end-1]; categorical=true)

    fig = Figure(size = (1200, 900), fontsize = 18)
    axs = Array{Any, 2}(undef, (3, 3))
    contours = Array{Any, 2}(undef, (3, 3))
    for (icol, (basin_key, basin)) in enumerate(pairs(basins))

        for (irow, x3D) in enumerate((ϕ.north, ϕbis.north, ϕter.north))

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
    ∬ = "∬"
    dx = rich("d", rich("x", font=:italic))
    dz = rich("d", rich("z", font=:italic))
    Label(fig[1, 0]; text = rich(∬, " vmo ", dx, " ", dz), rotation = π/2, tellheight = false)
    Label(fig[2, 0]; text = rich(∬, " vmo + GM ", dx, " ", dz), rotation = π/2, tellheight = false)
    Label(fig[3, 0]; text = rich(∬, " vmo + GM + submeso ", dx, " ", dz), rotation = π/2, tellheight = false)
    fig
    # save plot
    outputfile = joinpath(inputdir, "MOC_withGM_and_submeso.png")
    @info "Saving MOC as image file:\n  $(outputfile)"
    save(outputfile, fig)



# end


