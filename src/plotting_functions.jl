
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