"""
    area_polygon(vertices::Vector{<:Point2D})

The (absolute) area of the polygon.

Uses the shoelace formula: `½|Σ(xᵢyᵢ₊₁ - yᵢxᵢ₊₁)|`.
This is only valid for polygons that do not self-intersect.

For the signed area, use `first_moment`.
"""
function area_polygon(vertices::Path2D)
    abs(first_moment(vertices))
end

"""
    first_moment(vertices::Vector{<:Point2D})
    first_moment(polygon::Polygon)

The signed area of the polygon using the shoelace formula: `½Σ(xᵢyᵢ₊₁ - yᵢxᵢ₊₁)`.

The result is positive if the vertices are ordered counter-clockwise and negative if they are clockwise.

For a `Polygon`, the first moments of the holes are added to that of the exterior.
If the holes have the opposite orientation to the exterior, they will reduce the area.
"""
function first_moment(vertices::Path2D)
    # first moment of a simple polygon
    n = length(vertices)
    moment = 0.0
    for i in 1:n
        pt1 = vertices[i]
        pt2 = vertices[(i % n) + 1]
        moment += pt1[1] * pt2[2] - pt1[2] * pt2[1]
    end
    moment *= 0.5
    moment
end

"""
    is_clockwise(vertices::Vector{<:Point2D})
    is_clockwise(polygon::Polygon)

Return `true` if the vertices are ordered clockwise, determined by the sign of `first_moment`.
For a `Polygon`, only the exterior is checked.

A polygon with zero area is considered both clockwise and counter-clockwise.
"""
is_clockwise(vertices::Path2D) = first_moment(vertices) <= 0.0

"""
    is_counter_clockwise(vertices::Vector{<:Point2D})
    is_counter_clockwise(polygon::Polygon)

Return `true` if the vertices are ordered counter-clockwise, determined by the sign of `first_moment`.
For a `Polygon`, only the exterior is checked.

A polygon with zero area is considered both clockwise and counter-clockwise.
"""
is_counter_clockwise(vertices::Path2D) = first_moment(vertices) >= 0.0

"""
    centroid_polygon(vertices::Vector{<:Point2D})

Uses the following formulas:
- Cx = (1/6A)Σ(xᵢ + xᵢ₊₁)(xᵢyᵢ₊₁ - yᵢxᵢ₊₁)
- Cy = (1/6A)Σ(yᵢ + yᵢ₊₁)(xᵢyᵢ₊₁ - yᵢxᵢ₊₁)
"""
function centroid_polygon(vertices::Path2D)
    area = first_moment(vertices)
    Cx, Cy = _centroid_loop(vertices)
    centroid = (Cx/(6 * area), Cy/(6 * area))
    centroid
end

function _centroid_loop(vertices::Path2D)
    Cx = 0.0
    Cy = 0.0
    n = length(vertices)
    for i in 1:n
        pt1 = vertices[i]
        pt2 = vertices[(i % n) + 1]
        Cx += (pt1[1] + pt2[1]) * (pt1[1] * pt2[2] - pt1[2] * pt2[1])
        Cy += (pt1[2] + pt2[2]) * (pt1[1] * pt2[2] - pt1[2] * pt2[1])
    end
    (Cx, Cy)
end
