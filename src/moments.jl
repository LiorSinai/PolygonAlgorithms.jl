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
function first_moment(vertices::Path2D{T}) where T
    # first moment of a simple polygon
    n = length(vertices)
    moment = zero(T)
    for i in 1:n
        pt1 = vertices[i]
        pt2 = vertices[(i % n) + 1]
        moment += pt1[1] * pt2[2] - pt1[2] * pt2[1]
    end
    moment *= 0.5
    moment
end

"""
    _first_moment_with_error(vertices; atol=default_atol)

The `first_moment` and the size of its error. The polygon has zero area if `|moment| ≤ error`.

The area of a thin line is `area = L * w ≈ (perimeter/2) * w`, 
where the perimeter is divided by 2 because the line loops back on itself.
Using `atol` as the width (thickness), the error term is `area = perimeter/ 2 * atol`

The perimeter is calculated with the Manhattan distance `|Δx| + |Δy|` for each edge.
This avoids square roots and is between 1 and √2 times the Euclidean perimeter.
"""
function _first_moment_with_error(vertices::Path2D{T}; atol::AbstractFloat=default_atol) where T
    n = length(vertices)    
    moment = zero(T)
    perimeter = zero(T)
    for i in 1:n
        pt1 = vertices[i]
        pt2 = vertices[(i % n) + 1]
        moment += pt1[1] * pt2[2] - pt1[2] * pt2[1]
        perimeter += abs(pt2[1] - pt1[1]) + abs(pt2[2] - pt1[2])
    end
    moment / 2, atol * perimeter / 2
end

"""
    is_clockwise(vertices::Vector{<:Point2D}; atol=default_atol)
    is_clockwise(polygon::Polygon; atol=default_atol)

Return `true` if the vertices are ordered clockwise, determined by the sign of `first_moment`.
For a `Polygon`, only the exterior is checked.

A polygon with zero area is not considered clockwise.
The area is zero if the polygon is thinner than `atol`. See `_first_moment_with_error`.
To only check the direction, use a very small `atol` such as `eps(T)`.
"""
function is_clockwise(vertices::Path2D; atol::AbstractFloat=default_atol)
    moment, error = _first_moment_with_error(vertices; atol=atol)
    moment < -error
end

"""
    is_counter_clockwise(vertices::Vector{<:Point2D}; atol=default_atol)
    is_counter_clockwise(polygon::Polygon; atol=default_atol)

Return `true` if the vertices are ordered counter-clockwise, determined by the sign of `first_moment`.
For a `Polygon`, only the exterior is checked.

A polygon with zero area is not considered counter-clockwise.
The area is zero if the polygon is thinner than `atol`. See `_first_moment_with_error`.
To only check the direction, use a very small `atol` such as `eps(T)`.
"""
function is_counter_clockwise(vertices::Path2D; atol::AbstractFloat=default_atol)
    moment, error = _first_moment_with_error(vertices; atol=atol)
    moment > error
end

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

function _centroid_loop(vertices::Path2D{T}) where T
    Cx = zero(T)
    Cy = zero(T)
    n = length(vertices)
    for i in 1:n
        pt1 = vertices[i]
        pt2 = vertices[(i % n) + 1]
        Cx += (pt1[1] + pt2[1]) * (pt1[1] * pt2[2] - pt1[2] * pt2[1])
        Cy += (pt1[2] + pt2[2]) * (pt1[1] * pt2[2] - pt1[2] * pt2[1])
    end
    (Cx, Cy)
end
