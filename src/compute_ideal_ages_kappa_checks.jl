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
using CairoMakie
using GeoMakie
using Interpolations
using OceanBasins
using NaNStatistics
using Format

include("plotting_functions.jl")

# model = "ACCESS-ESM1-5"
# model = "ACCESS-CM2"
model = "ACCESS1-3"

# CMIP_version = "CMIP6"
CMIP_version = "CMIP5"

experiment = "historical"
# experiment = "piControl"

time_window = "Jan1990-Dec1999"
# time_window = "Jan1071-Dec1100" # <- last 30 years of ACCESS-ESM1-5 piControl
# time_window = "Jan1420-Dec1449" # <- last 30 years of ACCESS-CM2 piControl


# Gadi directory for input files
inputdirfun(member) = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)"

# find all members for which the inputdir contains umo.nc, vmo.nc, mlotst.nc, volcello.nc, and areacello.nc
requiredvariables = ["umo", "vmo", "mlotst", "volcello", "areacello", "agessc"]
hasrequireddata(member, variable_name) = isfile(joinpath(inputdirfun(member), "$variable_name.nc"))
hasrequireddata(member) = all(variable_name -> hasrequireddata(member, variable_name), requiredvariables)
members = readdir("/scratch/xv83/TMIP/data/$model/$experiment")
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

κHs = 500.0 .* [1/2, 1, 2] # m^2/s
κVdeeps = 3e-5 .* [1/2, 1, 2] # m^2/s
fe = FormatExpr("{1:0.0e}")

# for member in members[dataavailability.has_it_all]
# for member in members[dataavailability.has_it_all], κH in κHs, κVdeep in κVdeeps
for member in members[dataavailability.has_it_all][1:3], κH in κHs, κVdeep in κVdeeps
    # member = last(members)

    inputdir = inputdirfun(member)
    κVdeep_str = "kVdeep" * format(fe, κVdeep)
    κH_str = "kH" * format(fe, κH)
    outputdir = joinpath(inputdir, κH_str, κVdeep_str)
    mkpath(outputdir)

    # Load umo, vmo, mlotst, volcello, and areacello
    umo_ds = open_dataset(joinpath(inputdir, "umo.nc"))
    vmo_ds = open_dataset(joinpath(inputdir, "vmo.nc"))
    mlotst_ds = open_dataset(joinpath(inputdir, "mlotst.nc"))
    volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
    areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))

    mlotst = mlotst_ds["mlotst"] |> Array{Float64}

    # Make ualldirs
    u = makeualldirections(; umo_ds, vmo_ds)

    # Make makemodelgrid
    modelgrid = makemodelgrid(; areacello_ds, volcello_ds, mlotst_ds)

    # Make indices
    indices = makeindices(modelgrid.v3D)

    # Make transport matrix
    @warn "using κVdeep = $κVdeep, κH = $κH"
    (; T, Tadv, TκH, TκVML, TκVdeep) = transportmatrix(; u, mlotst, modelgrid, indices,
        ρ = 1025.0,
        κH, # m^2/s
        κVML = 0.1, # m^2/s
        κVdeep, # m^2/s
    )

    # unpack model grid
    (; lon, lat, zt, v3D,) = modelgrid
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
            "kVdeep" => κVdeep,
            "kH" => κH,
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
            "kVdeep" => κVdeep,
            "kH" => κH,
            "units" => "yr",
        )
    )
    arrays = Dict(:age => Γoutyr_YAXArray, :lat => volcello_ds.lat, :lon => volcello_ds.lon)
    Γout_ds = Dataset(; volcello_ds.properties, arrays...)

    # Save Γinyr3D to netCDF file
    outputfile = joinpath(outputdir, "ideal_mean_age.nc")
    @info "Saving ideal mean age as netCDF file:\n  $(outputfile)"
    savedataset(Γin_ds, path = outputfile, driver = :netcdf, overwrite = true)

    # Save Γoutyr3D to netCDF file
    outputfile = joinpath(outputdir, "mean_reemergence_time.nc")
    @info "Saving mean reemergence time as netCDF file:\n  $(outputfile)"
    savedataset(Γout_ds, path = outputfile, driver = :netcdf, overwrite = true)

    # Plot the ideal mean

    depth = 1000
    # Interpolate `Γinyr3D` to the given `depth`
    itp = interpolate((lev, ), [Γinyr3D[:,:,i] for i in axes(Γinyr3D, 3)], Gridded(Linear()))
    Γinyr2D = itp(depth)
    title = "$model $experiment $member $(time_window) (κH=$κH, κVdeep=$κVdeep) ideal mean age (yr) at $depth m"
    # plot options
    colorrange = (0, 1500)
    colormap = :viridis
    # plot
    fig = Figure(size = (1200, 600), fontsize = 18)
    ax = Axis(fig[1,1]; title, xtickformat, ytickformat)
    plt = plotmap!(ax, Γinyr2D, modelgrid; colorrange, colormap)
    Colorbar(fig[1,2], plt, label="Ideal mean age (yr)")
    # save plot
    outputfile = joinpath(outputdir, "ideal_mean_age_v2.png")
    @info "Saving ideal mean age as image file:\n  $(outputfile)"
    save(outputfile, fig)

    # Plot comparison with agessc
    agessc_ds = open_dataset(joinpath(inputdir, "agessc.nc"))
    agessc3D = agessc_ds["agessc"] |> Array{Float64}
    itp = interpolate((lev, ), [agessc3D[:,:,i] for i in axes(agessc3D, 3)], Gridded(Linear()))
    agessc2D = itp(depth)
    fig = Figure(size = (1200, 1800), fontsize = 18)
    Γdown = rich("Γ", superscript("↓"))
    Γup = rich("Γ", superscript("↑"))
    title = rich("$model $experiment $member $(time_window) (κH=$κH, κVdeep=$κVdeep) ", Γdown, " (yr) at $depth m")
    colorrange = (0, 1500)
    colormap = :viridis
    ax = Axis(fig[1,1]; title, xtickformat, ytickformat)
    plt1 = plotmap!(ax, Γinyr2D, modelgrid; colorrange, colormap)
    title = "$model $experiment $member $(time_window) (κH=$κH, κVdeep=$κVdeep) agessc (yr) at $depth m"
    ax = Axis(fig[2,1]; title, xtickformat, ytickformat)
    plt2 = plotmap!(ax, agessc2D, modelgrid; colorrange, colormap)
    Colorbar(fig[1:2,2], plt1, label="Ideal mean age (yr)")
    ax = Axis(fig[3,1]; title, xtickformat, ytickformat)
    colorrange = (-500, 500)
    colormap = :RdBu
    plt3 = plotmap!(ax, Γinyr2D - agessc2D, modelgrid; colorrange, colormap)
    Colorbar(fig[3,2], plt3, label=rich("Δ", Γdown, " (yr)"))
    # save plot
    outputfile = joinpath(outputdir, "ideal_mean_age_maps_vs_agessc_$(depth)m.png")
    @info "Saving ideal mean age as image file:\n  $(outputfile)"
    save(outputfile, fig)




    modulo_cmap = [
        0.   30.    0.   30.
        15.  100.    0.  100.
        25.   30.    0.  100.
        26.    0.    0.  100.
        34.    0.    0.   60.
        50.    0.  100.  100.
        60.    0.   50.    0.
        75.  100.  100.    0.
        90.  100.    0.    0.
       100.   40.    0.    0.
    ]


    # Plot zonal averages

    basin_keys = (:ATL, :PAC, :IND)
    basin_strs = ("Atlantic", "Pacific", "Indian")
    basin_functions = (isatlantic, ispacific, isindian)
    basin_values = (reshape(f(lat[:], lon[:], OCEANS), size(lat)) for f in basin_functions)
    basins = (; (basin_keys .=> basin_values)...)
    basin_latlims_values = [clamp.((-5, +5) .+ extrema(lat[.!isnan.(v3D[:,:,1]) .& basin[:,:,1]]), -80, 80) for basin in basins]
    basin_latlims = (; (basin_keys .=> basin_latlims_values)...)

    levels = 0:100:1500
    colormap = cgrad(:viridis, length(levels); categorical=true)
    extendlow = nothing
    extendhigh = colormap[end]
    colormap = cgrad(colormap[1:end-1]; categorical=true)

    Δlevels = -500:100:500
    Δcolormap = cgrad(:RdBu, length(Δlevels) + 1; categorical=true)
    Δextendlow = Δcolormap[1]
    Δextendhigh = Δcolormap[end]
    Δcolormap = cgrad(Δcolormap[2:end-1]; categorical=true)

    # Plot Γ↓ zonal averages

    fig = Figure(size = (1200, 800), fontsize = 18)
    axs = Array{Any,2}(undef, (3, 3))
    contours = Array{Any,2}(undef, (3, 3))
    for (icol, (basin_key, basin)) in enumerate(pairs(basins))

        for (irow, x3D) in enumerate((Γinyr3D, agessc3D))

            x2D = zonalaverage(x3D, modelgrid; mask = basin)

            ax = Axis(fig[irow, icol],
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
                label = irow < 3, ticklabels = irow < 3,
                ticks = irow < 3, grid = false)
            hideydecorations!(ax,
                label = icol > 1, ticklabels = icol > 1,
                ticks = icol > 1, grid = false)


            axs[irow, icol] = ax
        end
    end

    cb = Colorbar(fig[1:2, 4], contours[1, 1];
        vertical = true, flipaxis = true,
        # ticks = (, cbarticklabelformat.(levels)),
        label = rich(Γdown, " (yr)"),
        )
    cb.height = Relative(0.666)

    for (icol, (basin_key, basin)) in enumerate(pairs(basins))

        irow = 3
        x2D = zonalaverage(Γinyr3D - agessc3D, modelgrid; mask = basin)

        ax = Axis(fig[irow, icol],
            backgroundcolor=:lightgray,
            xgridvisible=false, ygridvisible=false,
            ylabel = "depth (m)")

        X = dropdims(maximum(lat, dims=1), dims=1)
        Y = zt
        Z = x2D
        co = contourf!(ax, X, Y, Z;
            levels = Δlevels,
            colormap = Δcolormap,
            nan_color = :lightgray,
            extendlow = Δextendlow,
            extendhigh = Δextendhigh,
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
            label = irow < 3, ticklabels = irow < 3,
            ticks = irow < 3, grid = false)
        hideydecorations!(ax,
            label = icol > 1, ticklabels = icol > 1,
            ticks = icol > 1, grid = false)

        axs[irow, icol] = ax
    end

    cb = Colorbar(fig[3, 4], contours[3, 1];
        vertical = true, flipaxis = true,
        # ticks = (, cbarticklabelformat.(levels)),
        label = rich("Δ", Γdown, " (yr)"),
        )
    cb.height = Relative(1)

    for (icol, (basin_str, xlims)) in enumerate(zip(basin_strs, basin_latlims))
        Label(fig[0, icol], basin_str, fontsize=20, tellwidth=false)
        colsize!(fig.layout, icol, Auto(xlims[2] - xlims[1]))
    end
    Label(fig[1, 0], text = "Transport matrix", fontsize=20, tellheight=false, rotation=π/2)
    Label(fig[2, 0], text = "agessc", fontsize=20, tellheight=false, rotation=π/2)
    Label(fig[3, 0], text = "Difference", fontsize=20, tellheight=false, rotation=π/2)

    title = "$model $experiment $member $(time_window) (κH=$κH, κVdeep=$κVdeep) ideal age"
    Label(fig[-1, 1:3], text = title, fontsize=20, tellwidth=false)

    # text = rich("Upstream sweeping time, ", ΓupΩ, ", for Ω = $(LONGTEXT[Ωz]) $(LONGTEXT[Ωbasin])")
    # Label(f[-1, :]; text, fontsize=20)

    # for (irun, run) in enumerate(runs)
    #     prefix = irun > 1 ? "future " : ""
    #     # Label(f[irun + 1, 0], text="future $(LONGTEXT[run])", fontsize=20, tellheight=false, rotation=π/2)
    #     # Label(f[irun + 2, 0], text=LONGTEXT[run], fontsize=20, tellheight=false, rotation=π/2)
    #     (irun == 1) && continue
    #     Label(f[irun + 3, 0], text=LONGTEXT[run], fontsize=20, tellheight=false, rotation=π/2)
    # end
    # Label(f[5:6, -1], text="preindustrial-to-future change", fontsize=20, tellheight=false, rotation=π/2)


    # # rowgap!(fig.layout, 1, 10)
    # # rowgap!(fig.layout, 2, 5)
    # # rowgap!(fig.layout, 5)
    rowgap!(fig.layout, 10)
    rowgap!(fig.layout, 4, 15)
    # # rowgap!(fig.layout, 5, 10)

    colgap!(fig.layout, 10)
    # save plot
    outputfile = joinpath(outputdir, "ideal_age_ZAVGs.png")
    @info "Saving ideal age ZAVGs as image file:\n  $(outputfile)"
    save(outputfile, fig)







    # Plot Γ↑ zonal averages

    fig = Figure(size = (1200, 300), fontsize = 18)
    axs = Array{Any,2}(undef, (1, 3))
    contours = Array{Any,2}(undef, (1, 3))
    for (icol, (basin_key, basin)) in enumerate(pairs(basins))

        x2D = zonalaverage(Γoutyr3D, modelgrid; mask = basin)
        irow = 1

        ax = Axis(fig[irow, icol],
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
            label = irow < 1, ticklabels = irow < 1,
            ticks = irow < 1, grid = false)
        hideydecorations!(ax,
            label = icol > 1, ticklabels = icol > 1,
            ticks = icol > 1, grid = false)

        axs[irow, icol] = ax
    end

    cb = Colorbar(fig[1, 4], contours[1, 1];
        vertical = true, flipaxis = true,
        # ticks = (, cbarticklabelformat.(levels)),
        label = rich(Γup, " (yr)"),
        )
    cb.height = Relative(1)


    for (icol, (basin_str, xlims)) in enumerate(zip(basin_strs, basin_latlims))
        Label(fig[0, icol], basin_str, fontsize=20, tellwidth=false)
        colsize!(fig.layout, icol, Auto(xlims[2] - xlims[1]))
    end
    Label(fig[1, 0], text = "Transport matrix", fontsize=20, tellheight=false, rotation=π/2)

    title = "$model $experiment $member $(time_window) (κH=$κH, κVdeep=$κVdeep) reemergence time"
    Label(fig[-1, 1:3], text = title, fontsize=20, tellwidth=false)

    rowgap!(fig.layout, 10)
    colgap!(fig.layout, 10)

    # save plot
    outputfile = joinpath(outputdir, "reemergence_time_ZAVGs.png")
    @info "Saving reemergence time ZAVGs as image file:\n  $(outputfile)"
    save(outputfile, fig)


    # Plot mean age at the seafloor level
    Γinyrseafloor = seafloorvalue(Γinyr3D, wet3D)
    title = "$model $experiment $member $(time_window) (κH=$κH, κVdeep=$κVdeep) mean age at seafloor"
    # plot options
    colorrange = (0, 1500)
    colormap = :viridis
    # plot
    fig = Figure(size = (1200, 600), fontsize = 18)
    ax = Axis(fig[1,1]; title, xtickformat, ytickformat)
    plt = plotmap!(ax, Γinyrseafloor, modelgrid; colorrange, colormap)
    Colorbar(fig[1,2], plt, label=rich(Γup, " at seafloor (yr)"))
    # save plot
    outputfile = joinpath(outputdir, "mean_age_at_seafloor.png")
    @info "Saving ideal mean age at sea floor as image file:\n  $(outputfile)"
    save(outputfile, fig)



    # Plot reemergence time at the seafloor level
    Γoutyrseafloor = seafloorvalue(Γoutyr3D, wet3D)
    title = "$model $experiment $member $(time_window) (κH=$κH, κVdeep=$κVdeep) reemergence time at seafloor"
    # plot options
    colorrange = (0, 1500)
    colormap = :viridis
    # plot
    fig = Figure(size = (1200, 600), fontsize = 18)
    ax = Axis(fig[1,1]; title, xtickformat, ytickformat)
    plt = plotmap!(ax, Γoutyrseafloor, modelgrid; colorrange, colormap)
    Colorbar(fig[1,2], plt, label=rich(Γup, " at seafloor (yr)"))
    # save plot
    outputfile = joinpath(outputdir, "reemergence_time_at_seafloor.png")
    @info "Saving mean reemergence time at seafloor as image file:\n  $(outputfile)"
    save(outputfile, fig)


end


