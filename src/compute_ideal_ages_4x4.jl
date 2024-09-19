# using Pkg
# Pkg.activate(".")
# Pkg.instantiate()

# using OceanTransportMatrixBuilder
# using NetCDF
# using YAXArrays
# using DimensionalData
# using SparseArrays
# using LinearAlgebra
# using Unitful
# using Unitful: s, yr
# using CairoMakie
# using GeoMakie
# using Interpolations

# member = "r1i1p1f1"

# # Gadi directory for input files
# inputdir = "/scratch/xv83/TMIP/data/ACCESS-ESM1-5/historical/$member/Jan1990-Dec1999"

# # Load umo, vmo, mlotst, volcello, and areacello
# umo_ds = open_dataset(joinpath(inputdir, "umo.nc"))
# vmo_ds = open_dataset(joinpath(inputdir, "vmo.nc"))
# mlotst_ds = open_dataset(joinpath(inputdir, "mlotst.nc"))
# volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
# areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))

# mlotst = mlotst_ds["mlotst"] |> Array{Float64}

# # Make ualldirs
# u = makeualldirections(; umo_ds, vmo_ds)

# # Make makemodelgrid
# modelgrid = makemodelgrid(; areacello_ds, volcello_ds, mlotst_ds)

# # Make indices
# indices = makeindices(modelgrid.v3D)

# # Make transport matrix
# (; T, Tadv, TκH, TκVML, TκVdeep) = transportmatrix(; u, mlotst, modelgrid, indices,
#     ρ = 1025.0,
#     κH = 500.0, # m^2/s
#     κVML = 0.1, # m^2/s
#     κVdeep = 1e-5, # m^2/s
# )

# # unpack model grid
# (; lon, lat, zt, v3D,) = modelgrid
# lev = zt
# # unpack indices
# (; wet3D, N) = indices

# v = v3D[wet3D]

# @info "coarsening grid"
# di = dj = 4
# LUMP, SPRAY, wet3D_c, v_c = OceanTransportMatrixBuilder.lump_and_spray(wet3D, v; di, dj, dk=1)

# # surface mask
# issrf3D = copy(wet3D)
# issrf3D[:,:,2:end] .= false
# issrf = issrf3D[wet3D]
# # Ideal mean age Γ is governed by
# # 	∂Γ/∂t + T Γ = 1 - M Γ
# # where M is matrix mask of surface with short timescale (1s)
# sΓ = ones(size(v))
# T_c = LUMP * T * SPRAY
# issrf_c = LUMP * issrf .> 0
# M_c = sparse(Diagonal(issrf_c))
# sΓ_c = LUMP * sΓ
# @info "Solving ideal mean age"
# Γ_c = (T_c + M_c) \ sΓ_c
# Γ = SPRAY * Γ_c
# Γyr = ustrip.(yr, Γ .* s)
# Γyr3D = OceanTransportMatrixBuilder.as3D(Γyr, wet3D)

# (v' * Γyr) / sum(v)

# # Turn Γ into a YAXArray by rebuilding from volcello
# Γyr_YAXArray = rebuild(volcello_ds["volcello"];
#     data = Γyr3D,
#     dims = dims(volcello_ds["volcello"]),
#     metadata = Dict(
#         "origin" => "ideal_mean_age_$(di)x$(dj) computed from ACCESS-ESM1-5 historical $member Jan1990-Dec1999",
#         "units" => "yr",
#     )
# )
# arrays = Dict(:age => Γyr_YAXArray, :lat => volcello_ds.lat, :lon => volcello_ds.lon)
# Γ_ds = Dataset(; volcello_ds.properties, arrays...)

# # Save Γyr3D to netCDF file
# outputfile = joinpath(inputdir, "ideal_mean_age_$(di)x$(dj).nc")
# @info "Saving ideal mean age as netCDF file:\n  $(outputfile)"
# savedataset(Γ_ds, path = outputfile, driver = :netcdf, overwrite = true)

# # Plot the ideal mean

# depth = 1000
# # Interpolate `Γyr3D` to the given `depth`
itp = interpolate((lev, ), [Γyr3D[:,:,i] for i in axes(Γyr3D, 3)], Gridded(Linear()))
Γyr2D = itp(depth)

title = "ACCESS-ESM1-5 $(di)x$(dj) historical $member Jan1990-Dec1999 ideal mean age (yr) at $depth m"
colorrange = (0, 1500)

fig = Figure(size = (1200,600), fontsize=24)
ax = GeoAxis(fig[1,1]; title, xlabel="Longitude", ylabel="Latitude")
# plt = surface!(ax, Γyr3D[:,:,30]; colormap = :viridis)
plt = surface!(ax, lon, lat, Γyr2D; colormap = :viridis, colorrange)
# save figure
outputfile = joinpath(inputdir, "ideal_mean_age_$(di)x$(dj).png")
@info "Saving ideal mean age as image file:\n  $(outputfile)"
save(outputfile, fig)



# ilev = findfirst(lev .≥ depth)

fig = let
    lonv = mlotst_ds.lon_verticies |> Array
    # make sure quads are not too distorted
    loninsamewindow(l1, l2) = mod(l1 - l2 + 180, 360) + l2 - 180
    lon = mod.(mlotst_ds.lon[:,:] .+ 180, 360) .- 180
    lonv = loninsamewindow.(lonv, reshape(lon, (1, size(lon)...)))
    latv = mlotst_ds.lat_verticies |> Array
    quad_points = vcat([Point2{Float64}.(lonv[:, i, j], latv[:, i, j]) for i in axes(lonv, 2), j in axes(lonv, 3)]...)
    quad_faces = vcat([begin; j = (i-1) * 4 + 1; [j j+1 j+2; j+2 j+3 j]; end for i in 1:length(quad_points)÷4]...)
    colors_per_point = vcat(fill.(vec(Γyr2D), 4)...)

    fig = Figure(size = (1200, 600), fontsize = 18)
    ax = Axis(fig[1,1]; title)
    plt = mesh!(ax, quad_points, quad_faces; color = colors_per_point, shading = NoShading, colorrange)
    xlims!(ax, (-180, 180))
    ylims!(ax, (-90, 90))
    cl=lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=0.85)
    translate!(cl, 0, 0, 1000)

    Colorbar(fig[1,2], plt, label="Ideal mean age (yr)")

    # fig = Figure(resolution = (1200,600))
    # ax = GeoAxis(fig[1,1]; dest = "+proj=moll")
    # mesh!(ax, quad_points, quad_faces; color = colors_per_point, shading = NoShading)
    # cl=lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=0.85)
    # translate!(cl, 0, 0, 1000)
    fig
end
outputfile2 = joinpath(inputdir, "ideal_mean_age_$(di)x$(dj)_v2.png")
@info "Saving ideal mean age as image file:\n  $(outputfile2)"
save(outputfile2, fig)



