
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

function plotmap!(ax, x2D, modelgrid; kwargs...)

    # unpack modelgrid
    lonv = modelgrid.lon_vertices
    latv = modelgrid.lat_vertices
    lon = modelgrid.lon

    # make sure quads are not too distorted
    lon = mod.(modelgrid.lon .+ 180, 360) .- 180
    lonv = loninsamewindow.(lonv, reshape(lon, (1, size(lon)...)))

    # create quads
    quad_points = vcat([Point2{Float64}.(lonv[:, i, j], latv[:, i, j]) for i in axes(lonv, 2), j in axes(lonv, 3)]...)
    quad_faces = vcat([begin; j = (i-1) * 4 + 1; [j j+1 j+2; j+2 j+3 j]; end for i in 1:length(quad_points)÷4]...)
    colors_per_point = vcat(fill.(vec(x2D), 4)...)

    # create plot
    plt = mesh!(ax, quad_points, quad_faces; color = colors_per_point, shading = NoShading, kwargs...)
    xlims!(ax, (-180, 180))
    ylims!(ax, (-90, 90))

    # add coastlines
    cl=lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=0.85)
    translate!(cl, 0, 0, 100)

    # move the plot behind the grid so we can see them
    translate!(plt, 0, 0, -100)

    return plt
end

OCEANS = oceanpolygons()

function zonalaverage(x3D, modelgrid; mask = 1)
    # unpack modelgrid
    (; v3D) = modelgrid

    # create zonal average
    x2D = nansum(x3D .* v3D .* mask, dims = 1) ./ nansum(mask .* v3D, dims = 1)

    return dropdims(x2D, dims = 1)
end

zlim = (6000, 0)
ztick = 0:1000:6000
zticklabel = map(x -> string.(x), ztick)