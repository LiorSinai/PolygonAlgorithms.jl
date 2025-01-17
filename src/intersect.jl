"""
    do_intersect(segment1, segment2; rtol=1e-4)

Determine if the segments intersect. For the intersection point, use `intersect_geometry`. 
"""
function do_intersect(segment1::Segment2D, segment2::Segment2D; rtol::AbstractFloat=1e-4)
    o1 = get_orientation(segment1[1], segment1[2], segment2[1]; rtol=rtol)
    o2 = get_orientation(segment1[1], segment1[2], segment2[2]; rtol=rtol)
    o3 = get_orientation(segment2[1], segment2[2], segment1[1]; rtol=rtol)
    o4 = get_orientation(segment2[1], segment2[2], segment1[2]; rtol=rtol)

    # general case
    if (o1 != o2) && (o3 != o4)
        return true
    end

    # special cases -> co-linear points on the segment
    return ((o1 == COLINEAR && on_segment(segment2[1], segment1)) ||
        (o2 == COLINEAR && on_segment(segment2[2], segment1)) ||
        (o3 == COLINEAR && on_segment(segment1[1], segment2)) || 
        (o4 == COLINEAR && on_segment(segment1[2], segment2))
        )
end
    
"""
    intersect_geometry(segment1, segment2; atol=1e-6)

Returns the intersection point if it exists or nothing otherwise.
Use `do_intersect` for a quicker boolean test.
"""
function intersect_geometry(
        segment1::Segment2D, segment2::Segment2D; atol::AbstractFloat=1e-6
    )
    line1 = line_from_segment(segment1)
    line2 = line_from_segment(segment2)

    point = intersect_geometry(line1, line2; atol=atol)
    if isnothing(point)
        return nothing
    end
    
    on_segment1 = on_segment(point, segment1, true; atol=atol)
    on_segment2 = on_segment(point, segment2, true; atol=atol)
    if on_segment1 && on_segment2
        return point 
    end
    nothing        
end

function classify_intersection(segment::Segment2D, point::Point2D; atol::AbstractFloat=1e-6)
    # assumes point is on segment
    at_start = is_same_point(segment[1], point; atol=atol)
    at_end = is_same_point(segment[2], point; atol=atol)
    along = !(at_start || at_end)
    (at_start, at_end, along)
end

"""
    intersect_geometry(line1::Line2D, line2::Line2D; atol=1e-6)

Returns the intersection point if it exists or nothing if they are parallel or coincident.

Equation of line is `ax + by = c`.

Intersection is then the solution of:
```
|a1 b1 | | x | = | c1 |
|a2 b2 | | y |   | c2 |
```
"""
function intersect_geometry(line1::Line2D, line2::Line2D; atol::AbstractFloat=1e-6)
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
    intersect_edges(polygon1::Vector{<:Point2D}, polygon2::Vector{<:Point2D})

Find all points which lie on the intersection of the edges of the vertices given by `polygon1` and `polygon2`.

Time complexity is `O(nm)` where `n` and `m` are the number of vertices of polygon 1 and 2 respectively.
"""
function intersect_edges(
    polygon1::Polygon2D{T}, polygon2::Polygon2D{T}
    ; atol::AbstractFloat=1e-6
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
