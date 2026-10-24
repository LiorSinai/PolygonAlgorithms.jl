"""
    do_intersect(segment1, segment2; rtol=default_rtol(T))

Determine if the segments intersect. For the intersection point, use `intersect_geometry`. 
"""
function do_intersect(segment1::Segment2D{T}, segment2::Segment2D{T}; rtol::AbstractFloat=default_rtol(T)) where T
    o1 = get_orientation(segment1[1], segment1[2], segment2[1]; rtol=rtol)
    o2 = get_orientation(segment1[1], segment1[2], segment2[2]; rtol=rtol)
    o3 = get_orientation(segment2[1], segment2[2], segment1[1]; rtol=rtol)
    o4 = get_orientation(segment2[1], segment2[2], segment1[2]; rtol=rtol)

    # general case
    # Colinear points within tolerance are excluded
    # They are handled by the special cases instead.
    if (o1 != o2) && (o3 != o4) &&
        (o1 != COLINEAR) && (o2 != COLINEAR) && (o3 != COLINEAR) && (o4 != COLINEAR)
        return true
    end

    # special cases -> co-linear points on the segment
    return ((o1 == COLINEAR && on_segment(segment2[1], segment1; rtol=rtol)) ||
        (o2 == COLINEAR && on_segment(segment2[2], segment1; rtol=rtol)) ||
        (o3 == COLINEAR && on_segment(segment1[1], segment2; rtol=rtol)) || 
        (o4 == COLINEAR && on_segment(segment1[2], segment2; rtol=rtol))
        )
end
    
"""
    intersect_geometry(segment1, segment2; rtol=default_rtol(T))

Returns the intersection point if it exists or nothing otherwise.
Use `do_intersect` for a quicker boolean test.

Equation of each segment:
```
f(t) = P1 + t(P2 - P1)
g(u) = Q1 + u(Q2 - Q1)
```

Intersection is when `f(t) = g(u)`:
```
(P2-P1)t + (Q2-Q1)u = Q1 - P1
```
where both `0≤ t ≤1` and `0≤ u ≤1`. This can be solved as a linear equation:
```
[t;u] = ΔP \\ (Q1 - P1)
```

Whether the segments intersect is decided by the orientation predicates, as in `do_intersect`,
so that the two functions always agree:
1. If both end points of one segment are strictly on the same side of the other, there is no intersection.
2. If both end points of one segment are colinear with the other, the segments are parallel
and `nothing` is returned. Use the orientations to check if they overlap.
3. If an end point is colinear with the other segment, it is the intersection and is returned exactly.
4. Otherwise the segments cross and `t` is used to calculate the point. It is clamped to `[0, 1]`,
so the returned point lies on `segment1`.
"""
function intersect_geometry(
    segment1::Segment2D{T}, segment2::Segment2D{T};
    rtol::AbstractFloat=default_rtol(T)
    ) where T
    o1 = get_orientation(segment1[1], segment1[2], segment2[1]; rtol=rtol)
    o2 = get_orientation(segment1[1], segment1[2], segment2[2]; rtol=rtol)
    o3 = get_orientation(segment2[1], segment2[2], segment1[1]; rtol=rtol)
    o4 = get_orientation(segment2[1], segment2[2], segment1[2]; rtol=rtol)
    if (o1 == o2 && o1 != COLINEAR) || (o3 == o4 && o3 != COLINEAR)
        return nothing # both end points are on the same side
    end
    if (o1 == COLINEAR && o2 == COLINEAR) || (o3 == COLINEAR && o4 == COLINEAR)
        return nothing # parallel lines or coincident
    end
    # an end point on the other segment is the intersection
    o1 == COLINEAR && return segment2[1]
    o2 == COLINEAR && return segment2[2]
    o3 == COLINEAR && return segment1[1]
    o4 == COLINEAR && return segment1[2]
    # A = [P2-P1  Q1-Q2]
    A11 = segment1[2][1] - segment1[1][1] # ΔPx
    A12 = segment2[1][1] - segment2[2][1] # -ΔQx
    A21 = segment1[2][2] - segment1[1][2] # ΔPy
    A22 = segment2[1][2] - segment2[2][2] # -ΔQy
    determinant = A11 * A22 - A12 * A21 # = cross_product((ΔPx, ΔPy), (0.0, 0.0), (ΔQx, ΔQy))
    b1 = segment2[1][1] - segment1[1][1]
    b2 = segment2[1][2] - segment1[1][2]
    # tu = A \ b
    t = (A22 * b1 - A12 * b2) / determinant
    # clamp to on segment
    t = clamp(t, zero(t), one(t))
    x = segment1[1][1] + t * A11
    y = segment1[1][2] + t * A21
    (x, y)
end

"""
    classify_intersection(segment, point; rtol=default_rtol(T))

Classify where `point` lies on `segment`. Returns a tuple of booleans `(at_start, at_end, along)`,
where `along` is true if the point is strictly between the end points.

Assumes `point` is on `segment` and the segment is ordered left to right (see `_compare_points`).
An intersection point calculated on another segment can be a few ulps past the end of this one.
So a point at or before the start is classified as the start and a point at or after the end as the end.
"""
function classify_intersection(segment::Segment2D{T}, point::Point2D{T}; rtol::AbstractFloat=default_rtol(T)) where T
    at_start = is_same_point(segment[1], point; rtol=rtol) || _compare_points(point, segment[1]; rtol=rtol) <= 0
    at_end = !at_start &&
        (is_same_point(segment[2], point; rtol=rtol) || _compare_points(point, segment[2]; rtol=rtol) >= 0)
    along = !(at_start || at_end)
    (at_start, at_end, along)
end

"""
    intersect_geometry(line1::Line2D, line2::Line2D)

Returns the intersection point if it exists or nothing if they are parallel or coincident.

Lines are parallel if the determinant is within floating point error: `|det| ≤ default_rtol(T) * (|a1 * b2| + |a2 * b1|)`.

Equation of line is `ax + by = c`.

Intersection is then the solution of:
```
|a1 b1 | | x | = | c1 |
|a2 b2 | | y |   | c2 |
```
"""
function intersect_geometry(line1::Line2D{T}, line2::Line2D{T}) where T
    a1, b1, c1 = line1.a, line1.b, line1.c
    a2, b2, c2 = line2.a, line2.b, line2.c
    determinant = a1 * b2 - a2 * b1
    if abs(determinant) <= default_rtol(T) * (abs(a1 * b2) + abs(a2 * b1))
        return nothing # parallel lines or coincident
    end
    x = (b2 * c1 - b1 * c2) / determinant
    y = (a1 * c2 - c1 * a2) / determinant
    (x, y)
end

"""
    intersect_edges(polygon1, polygon2; rtol=default_rtol(T))

Find all points which lie on the intersection of the edges of the vertices given by `polygon1` and `polygon2`.

Time complexity is `O(nm)` where `n` and `m` are the number of vertices of polygon 1 and 2 respectively.
"""
function intersect_edges(
    polygon1::Path2D{T}, polygon2::Path2D{T}
    ; rtol::AbstractFloat=default_rtol(T)
    ) where T
    points = Point2D{T}[]
    n = length(polygon1)
    m = length(polygon2)
    for i in 1:n
        edge1 = (polygon1[i], polygon1[i % n + 1])
        for j in 1:m
            edge2 = (polygon2[j], polygon2[j % m + 1])
            p = intersect_geometry(edge1, edge2; rtol=rtol)
            if !isnothing(p)
                push!(points, (p[1], p[2]))
            end
        end
    end
    points
end
