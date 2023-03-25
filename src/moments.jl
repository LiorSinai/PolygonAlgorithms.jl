"""
    area_polygon(vertices::Vector{<:Point2D})

Uses the shoelace formula: Σ(xᵢyᵢ₊₁ - yᵢxᵢ₊₁)
"""
function area_polygon(vertices::Polygon2D)
    abs(first_moment(vertices))
end

function first_moment(vertices::Polygon2D)
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

is_clockwise(vertices::Polygon2D) = first_moment(vertices) <= 0.0
is_counter_clockwise(vertices::Polygon2D) = first_moment(vertices) >= 0.0

"""
  centroid_polygon(vertices::Vector{<:Point2D})

Uses the following formulas:
- Cx = (1/6A)Σ(xᵢ + xᵢ₊₁)(xᵢyᵢ₊₁ - yᵢxᵢ₊₁)
- Cy = (1/6A)Σ(yᵢ + yᵢ₊₁)(xᵢyᵢ₊₁ - yᵢxᵢ₊₁)   
"""
function centroid_polygon(vertices::Polygon2D)
    area = first_moment(vertices)
    n = length(vertices)
    Cx = 0.0
    Cy = 0.0
    for i in 1:n
        pt1 = vertices[i]
        pt2 = vertices[(i % n) + 1]
        Cx += (pt1[1] + pt2[1]) * (pt1[1] * pt2[2] - pt1[2] * pt2[1])
        Cy += (pt1[2] + pt2[2]) * (pt1[1] * pt2[2] - pt1[2] * pt2[1])
    end
    centroid = (Cx/(6 * area), Cy/(6 * area))
    centroid 
end
