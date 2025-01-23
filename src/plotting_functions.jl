
function lonticklabel(lon)
    lon = mod(lon + 180, 360) - 180
    lon = isinteger(lon) ? Int(lon) : lon
    if lon == 0
        "0°"
    elseif (lon ≈ 180) || (lon ≈ -180)
        "180°"
    elseif lon > 0
        "$(string(lon))°E"
    else
        "$(string(-lon))°W"
    end
end
xtickformat(x) = lonticklabel.(x)

function latticklabel(lat)
    lat = isinteger(lat) ? Int(lat) : lat
    if lat == 0
        "0°"
    elseif lat > 0
        "$(string(lat))°N"
    else
        "$(string(-lat))°S"
    end
end
ytickformat(y) = latticklabel.(y)

loninsamewindow(l1, l2) = mod(l1 - l2 + 180, 360) + l2 - 180

function plotmap!(ax, x2D, gridmetrics; colorrange, colormap, levels=nothing)

    # unpack gridmetrics
    lonv = gridmetrics.lon_vertices
    latv = gridmetrics.lat_vertices
    lon = gridmetrics.lon

    # make sure quads are not too distorted
    lon = mod.(lon .+ -20, 360) .- -20
    ilon = sortperm(lon[:, 1])
    x2D = x2D[ilon, :]
    lon = lon[ilon, :]
    lonv = lonv[:, ilon, :]
    latv = latv[:, ilon, :]


    lonv = loninsamewindow.(lonv, reshape(lon, (1, size(lon)...)))

    # create quads
    quad_points = reduce(vcat, Point2{Float64}.(lonv[:, i, j], latv[:, i, j]) for i in axes(lonv, 2), j in axes(lonv, 3))
    quad_faces = reduce(vcat, begin; j = (i-1) * 4 + 1; [j j+1 j+2; j+2 j+3 j]; end for i in 1:length(quad_points)÷4)
    # quad_faces = [
    #     1 2 3;
    #     3 4 1;
    # ]
    colors_per_point = vcat(fill.(vec(x2D), 4)...)
    colors_per_face = vec(x2D)
    quad_faces = vcat([begin; j = (i-1) * 4 + 1; QuadFace(j, j+1, j+2, j+4); end for i in 1:length(quad_points)÷4]...)
    FT = eltype(quad_faces)
    N = length(quad_faces)
    cs = FaceView(colors_per_face, [FT(i) for i in 1:N])
    # generate normals per face (this creates a FaceView as well)
    ns = face_normals(quad_points, quad_faces)
    mymesh = Makie.GeometryBasics.mesh(quad_points, quad_faces, normal = ns, color = cs)


    # create plot
    # plt = mesh!(ax, quad_points, quad_faces; color = colors_per_point, shading = NoShading, colormap, colorrange, rasterize = 2)
    plt = mesh!(ax, mymesh; shading = NoShading, colormap, colorrange, rasterize = 2)
    xlims!(ax, (-180, 180))
    ylims!(ax, (-90, 90))

    # Add contourlines if levels is present
    if !isnothing(levels)
        lon2 = mod.(gridmetrics.lon .- 80, 360) .+ 80
        contourlevels = Contour.contours(lon2, gridmetrics.lat, x2D, levels)
        for cl in Contour.levels(contourlevels)
            lvl = level(cl) # the z-value of this contour level
            for line in Contour.lines(cl)
                xs, ys = coordinates(line) # coordinates of this line segment
                ls = lines!(ax, xs, ys; color = (:black, 0.5), linewidth = 1)
                translate!(ls, 0, 0, 110) # draw the contours above all
            end
        end
    end


    # add coastlines
    cl = if ax isa GeoAxis
        lines!(ax, GeoMakie.coastlines(); color = :black)
    else
        poly!(ax, GeoMakie.land(); color = :lightgray, strokecolor = :black, strokewidth = 1)
    end
    translate!(cl, 0, 0, 100)

    # move the plot behind the grid so we can see them
    translate!(plt, 0, 0, -100)

    return plt
end

function plotcontourfmap!(ax, x2D, gridmetrics; kwargs...)

    # unpack gridmetrics
    lon = gridmetrics.lon
    lat = gridmetrics.lat

    # extend west, east, and north
    extend(x2D) = [[transpose(x2D[end,:]); x2D; transpose(x2D[1, :])] [x2D[1,end]; x2D[end:-1:1, end]; x2D[end, end]]]
    extendlon(x2D) = [[transpose(x2D[end,:]) .- 360; x2D; transpose(x2D[1, :]) .+ 360] [x2D[1,end] - 360; x2D[end:-1:3end÷4+1, end] .- 360; x2D[3end÷4:-1:end÷4+1, end]; x2D[end÷4:-1:1, end] .+ 360 ; x2D[end, end] + 360]]
    # extend west and east only
    # extend(x2D) = [transpose(x2D[end,:]); x2D; transpose(x2D[1, :])]
    # extendlon(x2D) = [transpose(x2D[end,:]) .- 360; x2D; transpose(x2D[1, :]) .+ 360]
    # TODO delete: helper functions
    # extendnorth(x2D) = [x2D x2D[end:-1:1, end]]
    # extendeastwest(x2D) = [transpose(x2D[end,:]); x2D; transpose(x2D[1, :])]

    x2D = extend(x2D)
    lon = mod.(lon .- 80, 360) .+ 80
    lon = extendlon(lon)
    lat = extend(lat)


    plt = contourf!(ax, lon, lat, x2D; kwargs...)
    # xlims!(ax, (-180, 180))
    # ylims!(ax, (-90, 90))

    # add coastlines
    cl=lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=0.85)
    translate!(cl, 0, 0, 100)

    # move the plot behind the grid so we can see them
    translate!(plt, 0, 0, -100)

    return plt
end

function kinkline(A, B, c=2/3)
    # A -- C
    #       \
    #        B
    return [(0.9A[1] + 0.1B[1], A[2]), (c * B[1] + (1-c) * A[1], A[2]), B]
end



function volumeintegral(x, gridmetrics; mask = 1, dim = Tuple(1:ndims(x)))
    # unpack gridmetrics
    (; v3D) = gridmetrics

    mask3D = mask .* v3D

    ∫xdV = nansum(x .* mask3D; dim) |> Array

    # re-add NaNs where there is no water?
    ∫dV = nansum(mask3D; dim) |> Array
    ∫xdV[∫dV .== 0] .= NaN

    return ∫xdV
end

function average(x, gridmetrics; mask = 1, dims = Tuple(1:ndims(x)))
    # unpack gridmetrics
    (; v3D) = gridmetrics

    # TODO explain what purpose the dim stuff below serves
    # average x over dims
    # (Note special care for summing over the mask * v3D)
    mask3D = mask .* v3D
    maskdims = Tuple(dims ∩ Tuple(1:ndims(mask3D)))
    repeateddims = setdiff(dims, maskdims)
    multiplier = prod(size(x)[repeateddims])

    xmean = nansum(x .* mask3D; dims) ./ (multiplier * nansum(mask3D; dims = maskdims))

    return dropdims(xmean; dims) |> Array
end

zonalaverage(x, gridmetrics; mask = 1) = average(x, gridmetrics; mask, dims = 1)
horizontalaverage(x, gridmetrics; mask = 1) = average(x, gridmetrics; mask, dims = (1, 2))

zlim = (6000, 0)
ztick = 0:1000:6000
zticklabel = map(x -> string.(x), ztick)

function myhidexdecorations!(ax, condition)
    hidexdecorations!(ax,
        label = condition, ticklabels = condition,
        ticks = condition, grid = false
    )
end
function myhideydecorations!(ax, condition)
    hideydecorations!(ax,
        label = condition, ticklabels = condition,
        ticks = condition, grid = false
    )
end


# find the seafloor level
function seafloorvalue(x3D, wet3D, i, j)
    k = findlast(wet3D[i,j,:])
    isnothing(k) ? NaN : x3D[i,j,k]
end
function seafloorvalue(x3D, wet3D)
    [seafloorvalue(x3D, wet3D, i, j) for i in axes(x3D, 1), j in axes(x3D, 2)]
end
function seafloorvalue1D(x3D, wet3D)
    sfv = seafloorvalue(x3D, wet3D)
    sfv[.!isnan.(sfv)]
end
function seafloormask(wet3D)
    mask = falses(size(wet3D))
    for i in axes(wet3D, 1), j in axes(wet3D, 2)
        k = findlast(wet3D[i,j,:])
        isnothing(k) && continue
        mask[i,j,k] = true
    end
    mask
end

function seafloorvalue(x3D, wet3D, gridmetrics)
    (; Z3D, thkcello) = gridmetrics
    zseafloor2D = nansum(thkcello, dim = 3)
    map(eachindex(IndexCartesian(), zseafloor2D)) do I
        i, j = Tuple(I)
        zseafloor = zseafloor2D[i, j]
        (zseafloor == 0) && return NaN
        ks = wet3D[i,j,:]
        z = Z3D[i,j,ks]
        vals = x3D[i,j,ks]
        etp = linear_interpolation(z, vals; extrapolation_bc=Line())
        return etp(zseafloor)
    end
end



labeloptions = (
    font = :bold,
    align = (:left, :top),
    offset = (5, -2),
    space = :relative,
    fontsize = 24
)

