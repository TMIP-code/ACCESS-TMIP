"""
    (; lon, lat, areacello, lon_vertices, lat_vertices) = supergrid(model; dims)

returns the longitude, latitude, area, and vertex coordinates from the model's supergrid.
"""
function supergrid(model::String; dims)

    # ACCESS-OM2-1 => mom1deg.nc
    # ACCESS-OM2-025 => mom025deg.nc
    # ACCESS-OM2-01 => mom01deg.nc
    modelsupergridfile = "mom$(split(model, "-")[end])deg.nc"
    gridfile = joinpath("/g/data/xp65/public/apps/access_moppy_data/grids", modelsupergridfile)

    # Load data
    supergrid_ds = open_dataset(gridfile)
    superarea = readcubedata(supergrid_ds.area)
    lon = readcubedata(supergrid_ds.x)[2:2:end, 2:2:end]
    lat = readcubedata(supergrid_ds.y)[2:2:end, 2:2:end]
    areacello = YAXArray(
        dims,
        [sum(superarea[i:(i + 1), j:(j + 1)]) for i in 1:2:size(superarea, 1) , j in 1:2:size(superarea, 2)],
        Dict("name" => "areacello", "units" => "m^2"),
    )

    # Build vertices from supergrid
    # Dimensions of vertices ar (vertex, x, y)
    # Note to self: NCO shows it as (y, x, vertex)
    SW(x) = x[1:2:(end - 2), 1:2:(end - 2)]
    SE(x) = x[3:2:end, 1:2:(end - 2)]
    NE(x) = x[3:2:end, 3:2:end]
    NW(x) = x[1:2:(end - 2), 3:2:end]
    (nx, ny) = size(lon)
    vertices(x) = [
        reshape(SW(x), (1, nx, ny))
        reshape(SE(x), (1, nx, ny))
        reshape(NE(x), (1, nx, ny))
        reshape(NW(x), (1, nx, ny))
    ]
    lon_vertices = vertices(supergrid_ds.x)
    lat_vertices = vertices(supergrid_ds.y)

    return (; lon, lat, areacello, lon_vertices, lat_vertices)
end