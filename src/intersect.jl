"""
    do_intersect(segment1, segment2; atol=default_atol)

Determine if the segments intersect. For the intersection point, use `intersect_geometry`. 
"""
function do_intersect(segment1::Segment2D, segment2::Segment2D; rtol::AbstractFloat=default_rtol, atol::AbstractFloat=default_atol)
    o1 = get_orientation(segment1[1], segment1[2], segment2[1]; atol=atol)
    o2 = get_orientation(segment1[1], segment1[2], segment2[2]; atol=atol)
    o3 = get_orientation(segment2[1], segment2[2], segment1[1]; atol=atol)
    o4 = get_orientation(segment2[1], segment2[2], segment1[2]; atol=atol)

    # general case
    if (o1 != o2) && (o3 != o4)
        return true
    end

    # special cases -> co-linear points on the segment
    return ((o1 == COLINEAR && on_segment(segment2[1], segment1; atol=atol)) ||
        (o2 == COLINEAR && on_segment(segment2[2], segment1; atol=atol)) ||
        (o3 == COLINEAR && on_segment(segment1[1], segment2; atol=atol)) || 
        (o4 == COLINEAR && on_segment(segment1[2], segment2; atol=atol))
        )
end
    
"""
    intersect_geometry(segment1, segment2; atol=default_atol; rtol=default_rtol)

Returns the intersection point if it exists or nothing otherwise.
Use `do_intersect` for a quicker boolean test.

Equation of each segment:
```
f(t) = P1 + t(P2 - P1)
g(u) = P3 + u(P4 - P3)
```

Intersection is when `f(t) = g(u)`:
```
(P2-P1)t + (P3-P4)u = P3 - P1
```
where both `0≤ t ≤1` and `0≤ u ≤1`. This can be solved as a linear equation:
```
[t;u] = ΔP \\ (P3 - P1)
```
"""
function intersect_geometry(
    segment1::Segment2D, segment2::Segment2D; 
    atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_rtol
    )
    A11 = segment1[2][1] - segment1[1][1]
    A12 = segment2[1][1] - segment2[2][1]
    A21 = segment1[2][2] - segment1[1][2]
    A22 = segment2[1][2] - segment2[2][2]
    determinant = A11 * A22 - A12 * A21
    if abs(determinant) < atol
        return nothing # parallel lines or coincident
    end
    b1 = segment2[1][1] - segment1[1][1]
    b2 = segment2[1][2] - segment1[1][2]
    # tu = A \ b
    t = (A22 * b1 - A12 * b2) / determinant
    u = (A11 * b2 - A21 * b1) / determinant
    if -rtol <= t <= 1 && -rtol <= u <= 1
        x = segment1[1][1] + t * A11
        y = segment1[1][2] + t * A21
        return (x, y)
    else
        return nothing
    end
end

function classify_intersection(segment::Segment2D, point::Point2D; atol::AbstractFloat=default_atol)
    # assumes point is on segment
    at_start = is_same_point(segment[1], point; atol=atol)
    at_end = is_same_point(segment[2], point; atol=atol)
    along = !(at_start || at_end)
    (at_start, at_end, along)
end

"""
    intersect_geometry(line1::Line2D, line2::Line2D; atol=default_atol)

Returns the intersection point if it exists or nothing if they are parallel or coincident.

Equation of line is `ax + by = c`.

Intersection is then the solution of:
```
|a1 b1 | | x | = | c1 |
|a2 b2 | | y |   | c2 |
```
"""
function intersect_geometry(line1::Line2D, line2::Line2D; atol::AbstractFloat=default_atol)
    a1, b1, c1 = line1.a, line1.b, line1.c
    a2, b2, c2 = line2.a, line2.b, line2.c
    determinant = a1 * b2 - a2 * b1
    if abs(determinant) < atol
        return nothing # parallel lines or coincident
    end
    x = (b2 * c1 - b1 * c2) / determinant
    y = (a1 * c2 - c1 * a2) / determinant
    (x, y)
end

"""
    intersect_edges(polygon1, polygon2; atol=default_atol)

Find all points which lie on the intersection of the edges of the vertices given by `polygon1` and `polygon2`.

Time complexity is `O(nm)` where `n` and `m` are the number of vertices of polygon 1 and 2 respectively.
"""
function intersect_edges(
    polygon1::Path2D{T}, polygon2::Path2D{T}
    ; atol::AbstractFloat=default_atol
    ) where T
    points = Point2D{T}[]
    n = length(polygon1)
    m = length(polygon2)
    for i in 1:n
        edge1 = (polygon1[i], polygon1[i % n + 1])
        for j in 1:m
            edge2 = (polygon2[j], polygon2[j % m + 1])
            p = intersect_geometry(edge1, edge2; atol=atol)
            if !isnothing(p)
                push!(points, (p[1], p[2]))
            end
        end
    end
    points
end
