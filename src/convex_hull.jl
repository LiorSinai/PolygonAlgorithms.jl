
abstract type ConvexHullAlgorithm end

struct GiftWrappingAlg <: ConvexHullAlgorithm end
struct GrahamScanAlg   <: ConvexHullAlgorithm end

"""
    convex_hull(points, alg=GiftWrappingAlg(); atol=default_atol)

Determine the indices of the convex hull for a set of points.

`alg` can either be `GiftWrappingAlg()` or `GrahamScanAlg()`.

For  `n` input vertices and `h` resultant vertices on the convex hull:
- `GiftWrappingAlg` runs in `O(nh)` time.
- `GrahamScanAlg` runs in `O(n*log(n))` time.
"""
convex_hull(points::Polygon2D; options...) = convex_hull(points, GiftWrappingAlg(); options...)

function convex_hull(points::Polygon2D, ::GiftWrappingAlg; atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_rtol)
    # https://www.geeksforgeeks.org/convex-hull-using-jarvis-algorithm-or-wrapping/    
    topleft = left_topmost(points)
    hull_idxs = Int[]
    n = length(points)

    p = topleft
    q = 1
    while true
        push!(hull_idxs, p)
        q = (p % n + 1)

        point_p = points[p]
        point_q = points[q]
        # find more counter-clockwise point than q
        for i in eachindex(points)
            point_i = points[i]
            turn = get_orientation(point_p, point_i, point_q; atol=atol)
            if turn == COUNTER_CLOCKWISE
                q = i
                point_q = point_i
            elseif turn == COLINEAR
                dq = norm2(point_q, point_p)
                di = norm2(point_i, point_p)
                if di > dq
                    q = i
                    point_q = point_i
                end
            end
        end
        p = q

        if (p == topleft)
            break
        end
    end
    hull_idxs
end

function left_topmost(points::Polygon2D)
    idx = 1
    for i in eachindex(points)
        if points[i][1] < points[idx][1]
            idx = i
        elseif (points[i][1] == points[idx][1]) && (points[i][2] > points[idx][2])
            idx = i
        end
    end
    idx
end

function convex_hull(points::Polygon2D, ::GrahamScanAlg; atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_rtol)
    # https://www.geeksforgeeks.org/convex-hull-using-graham-scan/
    idx = bottom_leftmost(points)
    p0 = points[idx]
    idxs = sortperm(points, lt=(p, q) -> isless_orientation(p, q, p0; atol=atol))

    hull = idxs[[1, 2, 3]]
    if get_orientation(points[hull[1]], points[hull[2]], points[hull[3]]; atol=atol) == COLINEAR
        deleteat!(hull, 2)
    end
    for idx in idxs[4:end]
        while length(hull) > 1 &&
            get_orientation(points[hull[end-1]], points[hull[end]], points[idx]; atol=atol) != COUNTER_CLOCKWISE
            pop!(hull)
        end
        push!(hull, idx)
    end
    hull
end

function bottom_leftmost(points::Polygon2D)
    idx = 1
    for i in eachindex(points)
        if points[i][2] < points[idx][2]
            idx = i
        elseif (points[i][2] == points[idx][2]) && (points[i][1] < points[idx][1])
            idx = i
        end
    end
    idx
end
