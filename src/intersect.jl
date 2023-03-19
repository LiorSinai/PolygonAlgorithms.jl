"""
    do_intersect(segment1, segment2)

Determine if the segments intersect. For the intersection point, use `intersect_segments`. 
"""
function do_intersect(segment1::Segment2D, segment2::Segment2D)
    o1 = get_orientation(segment1[1], segment1[2], segment2[1])
    o2 = get_orientation(segment1[1], segment1[2], segment2[2])
    o3 = get_orientation(segment2[1], segment2[2], segment1[1])
    o4 = get_orientation(segment2[1], segment2[2], segment1[2])

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
    intersect_geometry(segment1, segment2)

Returns the intersection point if it exists or nothing otherwise.
Use `do_intersect` for a quicker boolean test.
"""
function intersect_geometry(
        segment1::Segment2D, segment2::Segment2D; atol=1e-6
    )
    line1 = line_from_segment(segment1)
    line2 = line_from_segment(segment2)

    point = intersect_geometry(line1, line2; atol=atol)
    if isnothing(point)
        return nothing
    end
    
    on_segment1 = on_segment(point, segment1; atol=atol, on_line=true)
    on_segment2 = on_segment(point, segment2; atol=atol, on_line=true)
    if on_segment1 && on_segment2
        return point 
    end
    nothing        
end

"""
    line_from_segment(segment)

Equation of line is `ax + by = c`.
"""
function line_from_segment(segment::Segment2D)
    a = segment[1][2] - segment[2][2]
    b = segment[2][1] - segment[1][1] 
    c = a * segment[1][1] + b * segment[1][2]
    Line2D(a, b, c)
end

"""
    intersect_geometry(line1, line2)

Returns the intersection point if it exists or nothing if they are parallel.

Equation of line is `ax + by = c`.

Intersection is then the solution of:
```
|a1 b1 | | x | = | c1 |
|a2 b2 | | y |   | c2 |
```
"""
function intersect_geometry(line1::Line2D, line2::Line2D; atol=1e-6)
    a1, b1, c1 = line1.a, line1.b, line1.c
    a2, b2, c2 = line2.a, line2.b, line2.c
    determinant = a1 * b2 - a2 * b1
    if abs(determinant) < atol
        return nothing # parallel lines
    end
    x = (b2 * c1 - b1 * c2) / determinant
    y = (a1 * c2 - c1 * a2) / determinant
    (negative_zero_to_zero(x), negative_zero_to_zero(y))
end

function negative_zero_to_zero(x::T) where T <: AbstractFloat
    # in the edge case that x is zero this will set comparisons to fail
    # Set([0.0]) == Set([-0.0]) is false
    # so just set it zero
    x == -0.0 ? zero(x) : x
end

"""
    intersect_edges(polygon1, polygon2)

Find all points which lie on the intersection of the edges of the vertices given by `polygon1` and `polygon2`.

Time complexity is `O(nm)` where `n` and `m` are the number of vertices of polygon 1 and 2 respectively.
"""
function intersect_edges(polygon1::Polygon{T}, polygon2::Polygon{T}) where T
    points = Point2D{T}[]
    n = length(polygon1)
    m = length(polygon2)
    for i in 1:n
        edge1 = (polygon1[i], polygon1[i % n + 1])
        for j in 1:m
            edge2 = (polygon2[j], polygon2[j % m + 1])
            p = intersect_geometry(edge1, edge2)
            if !isnothing(p)
                push!(points, (p[1], p[2]))
            end
        end
    end
    points
end
