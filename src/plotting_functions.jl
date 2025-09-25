
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

function divergingcbarticklabel(x)
    isinteger(x) && (x = Int(x))
    if x == 0
        "0"
    elseif x > 0
        "+" * string(x)
    else
        "−" * string(-x)
    end
end
divergingcbarticklabelformat(x) = divergingcbarticklabel.(x)

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

"""
To be used as `colorscale` of Makie's `contourf` and `heatmap`.
Picewise linear mapping such that
func.([v1,v2, . . . , vn]) == [0, 1, . . . , n-1].
Outside the range [v1, vn], it's the simple lienar extrapolation.
vs = [v1,v2, . . . , vn] must be strictly increasing.
from https://discourse.julialang.org/t/makie-nonlinear-color-levels-in-colorbar/118056/5
"""
function mk_piecewise_linear(vs)
    @assert length(vs) > 1
    function is_increasing(ss)
        prev = ss[1]
        for s in ss[2:end]
            (s ≤ prev) && return false
            prev = s
        end
        return true
    end
    @assert is_increasing(vs)
    d1 = vs[2] - vs[1]
    d2 = vs[end] - vs[end-1]
    un = size(vs, 1) - 1
    function piecewise_linear(v)
        if v <= vs[1]
            (v - vs[1]) / d1
        elseif v ≥ vs[end]
            (v - vs[end]) / d2 + un
        else
            i = findfirst(q -> v < q, vs) - 1
            d = vs[i + 1] - vs[i]
            (v - vs[i]) / d + i - 1
        end
    end
    function its_inverse(u)
        if u ≤ 0
            u * d1 + vs[1]
        elseif u ≥ un
            (u - un) * d2 + vs[end]
        else
            iu = floor(Int, u)
            i = iu + 1
            d = vs[i + 1] - vs[i]
            (u - iu) * d + vs[i]
        end
    end
    return ReversibleScale(piecewise_linear, its_inverse)
end


land1 = GeoMakie.land()
mapwindow = GeometryBasics.Polygon([
    Point2f(20, -90),
    Point2f(20 + 360, -90),
    Point2f(20 + 360, +90),
    Point2f(20, +90),
    Point2f(20, -90),
])
land1cut = [LibGEOS.intersection(p, mapwindow) for p in land1]
land1cut = [p for p in land1cut if !LibGEOS.isEmpty(p)]
land2 = GeometryOps.transform(P -> P + Point2f(360, 0), land1)
land2cut = [LibGEOS.intersection(p, mapwindow) for p in land2]
land2cut = [p for p in land2cut if !LibGEOS.isEmpty(p)]

function plotmap!(ax, x2D, gridmetrics; colorrange, colormap, levels=nothing, highclip = automatic, lowclip = automatic, colorscale = identity)

    # unpack gridmetrics
    lonv = gridmetrics.lon_vertices
    latv = gridmetrics.lat_vertices
    lon = gridmetrics.lon

    # make sure quads are not too distorted
    lon = mod.(gridmetrics.lon .+ -20, 360) .- -20
    lonv = loninsamewindow.(lonv, reshape(lon, (1, size(lon)...)))

    # create quads
    quad_points = vcat([Point2{Float64}.(lonv[:, i, j], latv[:, i, j]) for i in axes(lonv, 2), j in axes(lonv, 3)]...)
    quad_faces = vcat([begin; j = (i-1) * 4 + 1; [j j+1 j+2; j+2 j+3 j]; end for i in 1:length(quad_points)÷4]...)
    colors_per_point = vcat(fill.(vec(x2D), 4)...)

    # create plot
    plt = mesh!(ax, quad_points, quad_faces; color = colors_per_point, shading = NoShading, colormap, colorrange, rasterize = 2, highclip, lowclip, colorscale)
    xlims!(ax, (20, 20 + 360))
    ylims!(ax, (-90, 90))

    # Add contourlines if levels is present
    if !isnothing(levels)
        ilon = sortperm(lon[:,1])
        contour!(ax, lon[ilon, :], lat[ilon, :], x2D[ilon, :]; levels, color=:black, labels = true) # <- looks terrible
        # lon2 = mod.(gridmetrics.lon .- 80, 360) .+ 80
        # contourlevels = Contour.contours(lon2, gridmetrics.lat, x2D, levels)
        # for cl in Contour.levels(contourlevels)
        #     lvl = level(cl) # the z-value of this contour level
        #     for line in Contour.lines(cl)
        #         xs, ys = coordinates(line) # coordinates of this line segment
        #         ls = lines!(ax, xs, ys; color = (:black, 0.5), linewidth = 1)
        #         translate!(ls, 0, 0, 110) # draw the contours above all
        #     end
        # end
    end


    # add coastlines

    cl1 = poly!(ax, land1cut; color = :lightgray, strokecolor = :black, strokewidth = 1)
    translate!(cl1, 0, 0, 100)
    cl2 = poly!(ax, land2cut; color = :lightgray, strokecolor = :black, strokewidth = 1)
    translate!(cl2, 0, 0, 100)

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
        etp = linear_interpolation(z, vals; extrapolation_bc=Interpolations.Line())
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

# taken from MakieExtra (not using Pkg because it has outdated Makie dep)
using Makie.IntervalSets
Makie.project(s, r::HyperRectangle) = HyperRectangle(Makie.project(s, r.origin), Makie.project(s, r.origin + r.widths) - Makie.project(s, r.origin))
corner(r::HyperRectangle{2}, which::NTuple{2,Integer}) = Makie.Point(extrema(r)[_which_to_ix(which[1])][1], extrema(r)[_which_to_ix(which[2])][2])
_which_to_ix(which::Integer) = which == -1 ? 1 : which == 1 ? 2 : error("which must be -1 or 1, got $which")
fullproject(ax, p) = Makie.project(Makie.get_scene(ax), Makie.apply_transform(Makie.transform_func(ax), p)) + viewport(ax)[].origin
Base.:(⊆)(a::HyperRectangle, b::HyperRectangle) = all(map(⊆, intervals(a), intervals(b)))
intervals(r::HyperRectangle) = Interval.(r.origin, r.origin + r.widths)
function zoom_lines!(ax1, ax2; strokewidth=1.5, strokecolor=:black, color=(:black, 0), rectattrs=(;), lineattrs=(;))
    pscene = parent(parent(Makie.parent_scene(ax1)))
    @assert parent(parent(Makie.parent_scene(ax2))) === pscene
    obs = lift(ax1.finallimits, ax2.finallimits, ax1.scene.viewport, ax2.scene.viewport, ax1.scene.camera.projectionview, ax2.scene.camera.projectionview, Makie.transform_func(ax1), Makie.transform_func(ax2)) do _...
        lims = [ax1.finallimits[], ax2.finallimits[]]
        axs = lims[1] ⊆ lims[2] ? (ax1, ax2) :
              lims[2] ⊆ lims[1] ? (ax2, ax1) :
              nothing
        slines = if isnothing(axs)
            nothing
        else
            r1 = fullproject(axs[1], axs[1].finallimits[])
            r2 = fullproject(axs[2], axs[1].finallimits[])
            # cornsets = [
            #     ((corner(r1, (1,1)), corner(r2, (-1,1))), (corner(r1, (1,-1)), corner(r2, (-1,-1)))),
            #     ((corner(r1, (1,-1)), corner(r2, (1,1))), (corner(r1, (-1,-1)), corner(r2, (-1,1)))),
            #     ((corner(r1, (-1,-1)), corner(r2, (1,-1))), (corner(r1, (-1,1)), corner(r2, (1,1)))),
            #     ((corner(r1, (-1,1)), corner(r2, (-1,-1))), (corner(r1, (1,1)), corner(r2, (1,-1)))),
            # ]
            # argmin(cornsets) do ((a1, a2), (b1, b2))
            #     min(norm(a1-a2), norm(b1-b2))
            # end
            # BP: below is my zoom lines that does not work in general
            # cornsets2 = [
            #     ((corner(r1, (1,1)), corner(r2, (1,1))), (corner(r1, (-1,-1)), corner(r2, (-1,-1)))),
            #     ((corner(r1, (1,-1)), corner(r2, (1,-1))), (corner(r1, (-1,1)), corner(r2, (-1,1)))),
            # ]
            # argmin(cornsets2) do ((a1, a2), (b1, b2))
            #     max(norm(a1-a2), norm(b1-b2))
            # end
            [
                corner(r1, (1,1)), corner(r2, (1,1)),
                corner(r1, (-1,-1)), corner(r2, (-1,-1)),
                corner(r1, (1,-1)), corner(r2, (1,-1)),
                corner(r1, (-1,1)), corner(r2, (-1,1)),
            ]
            # BP: below is my zoom lines that does not work in general
            # (corner(r1, (1,-1)), corner(r2, (1,-1))), (corner(r1, (-1,1)), corner(r2, (-1,1)))
            # (corner(r1, (1,1)), corner(r2, (1,1))), (corner(r1, (-1,-1)), corner(r2, (-1,-1)))
        end
        (
            rect1=ax2.finallimits[],
            rect2=ax1.finallimits[],
            # slines=isnothing(slines) ? Point2{Float32}[] : Point2{Float32}[slines[1]..., slines[2]...],
            slines=isnothing(slines) ? Point2{Float32}[] : Point2{Float32}[slines...],
        )
    end

    rectattrs = (; strokewidth, strokecolor, color, xautolimits=false, yautolimits=false, rectattrs...)
    p1 = poly!(ax1, (@lift $obs.rect1); rectattrs...)
    p2 = poly!(ax2, (@lift $obs.rect2); rectattrs...)
    translate!(p1, 0, 0, -200)
    translate!(p2, 0, 0, -200)
    plt = linesegments!(pscene, (@lift $obs.slines); color=strokecolor, linewidth=strokewidth, linestyle=:dot, lineattrs...)
    translate!(plt, 0, 0, -200)
    return nothing
end


plotgridcell!(ax, lonv, latv; kwargs...) = poly!(ax, GeometryBasics.Polygon([Point2f(lon, lat) for (lon, lat) in zip(lonv, latv)]); kwargs...)
