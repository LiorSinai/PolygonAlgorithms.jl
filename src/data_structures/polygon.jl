import Base: ==

"""
    Polygon(exterior, holes)
    Polygon(exterior; holes=holes)

A representation of a Polygon with holes.
"""
struct Polygon{T<:AbstractFloat}
    exterior::Path2D{T}
    holes::Vector{<:Path2D{T}}
end

==(polygon1::Polygon, polygon2::Polygon) = 
    (polygon1.exterior == polygon2.exterior) && (polygon1.holes == polygon2.holes)

"""
    validate_polygon(polygon::Polygon)
    validate_polygon(exterior::Path2D; holes=Path2D{T}[]; atol=default_atol)
    
The following checks are done:
1. The exterior and each hole has 3 or more points.
2. The exterior and none of its holes self-intersect.
3. All the holes are contained within the exterior. They can touch at the vertex.

The following additional checks are not done:
1. Holes should not intersect each other.
2. No nested holes: holes should not be inside other holes.
"""
function validate_polygon(exterior::Path2D{T}; holes::Vector{<:Path2D}=Path2D{T}[], atol::AbstractFloat=default_atol) where T
    @assert length(exterior) > 2 "exterior requires at least 3 points."
    queue = convert_to_event_queue(exterior)
    @assert !any_intersect(queue; atol=atol, include_vertices=false) "Invalid exterior: edges self-intersect."
    for (idx, hole) in enumerate(holes)
        @assert length(hole) > 2 "Invalid hole $(idx): requires at least 3 points."
        queue_h = convert_to_event_queue(hole)
        @assert !any_intersect(queue_h; atol=atol, include_vertices=false) "Invalid hole $(idx): edges self-intersect."
        for event in queue
            insert_in_order!(queue_h, event; lt=compare_events)
        end
        @assert !any_intersect(queue_h; atol=atol, include_vertices=false) "Hole $(idx) intersects with the exterior."
        j = 1
        while (j < length(hole)) && on_border(exterior, hole[j])
            j += 1
        end
        @assert contains(exterior, hole[j]; atol=atol, on_border_is_inside=false) "Hole $(idx) is outside the polygon."
    end
    true
end

validate_polygon(polygon::Polygon) = validate_polygon(polygon.exterior, polygon.holes)

function Polygon(
    exterior::Path2D{T}
    ; holes::Vector{<:Path2D}=Path2D{T}[],
    validate::Bool=false,
    atol::AbstractFloat=default_atol) where T
    if validate
        validate_polygon(exterior, holes; atol=atol)
    end
    Polygon{T}(exterior, holes)
end


"""
    fully_contains(polygon1::Path2D, polygon2::Path2D)

A `polygon1` fully contains another `polygon2` if:
    1. None of their segments intersect. However, they can touch.
    2. At least one point of `polygon2` is inside `polygon1`.

The polygons are assumed to not self-intersect.

The algorithm runs in `O((n+m+k)log(n+m))` time where `n` and `m` are the number of vertices of `polygon1` 
and `polygon2` respectively and `k` is the total number of intersections 
"""
function fully_contains(polygon1::Path2D, polygon2::Path2D; atol::AbstractFloat=default_atol)
    queue = convert_to_event_queue(polygon1; atol=atol)
    events2 = convert_to_event_queue(polygon2; primary=false, atol=atol)
    for event in events2
        insert_in_order!(queue, event; lt=compare_events)
    end
    if any_intersect(queue; include_vertices=false)
        return false
    end
    j = 1
    while (j < length(polygon2)) && on_border(polygon1, polygon2[j])
        j += 1
    end
    contains(polygon1, polygon2[j]; atol=atol, on_border_is_inside=false)
end

"""
    bounds(polygon::Polygon)

A rectangle which bounds the polygon, given as `(xmin, ymin, xmax, ymax)`.
"""
bounds(polygon::Polygon) = bounds(polygon.exterior)

"""
    bounds(polygons::Polygon)

A rectangle which bounds all the polygons, given as `(xmin, ymin, xmax, ymax)`.
"""
bounds(polygons::Vector{<:Polygon}) = bounds(map(p -> p.exterior), polygons)

first_moment(polygon) = first_moment(polygon.exterior) + sum(first_moment, polygon.holes, init=0.0)

"""
    area_polygon(polygon::Polygon)

Uses the shoelace formula: Σ(xᵢyᵢ₊₁ - yᵢxᵢ₊₁). 

Holes are counted as negative area.
"""
function area_polygon(polygon::Polygon)
    area_polygon(polygon.exterior) - sum(map(area_polygon, polygon.holes), init=0.0)
end

is_clockwise(polygon::Polygon) = is_clockwise(polygon.exterior)
is_counter_clockwise(polygon::Polygon) = is_counter_clockwise(polygon.exterior)

"""
    centroid_polygon(polygon::Polygon)

Uses the following formulas:
- Cxⱼ = (1/6Aⱼ)Σ(xᵢ + xᵢ₊₁)(xᵢyᵢ₊₁ - yᵢxᵢ₊₁)
- Cyⱼ = (1/6Aⱼ)Σ(yᵢ + yᵢ₊₁)(xᵢyᵢ₊₁ - yᵢxᵢ₊₁)
- Cx = ΣCxⱼAⱼ/ΣAⱼ
- Cy = ΣCyⱼAⱼ/ΣAⱼ

Holes are counted as negative area.
"""
function centroid_polygon(polygon::Polygon)
    e_area = first_moment(polygon.exterior)
    Cex, Cey = _centroid_loop(polygon.exterior)
    Cx = Cex/(6 * e_area) * abs(e_area)
    Cy = Cey/(6 * e_area) * abs(e_area)
    area = abs(e_area)
    for hole in polygon.holes
        Chx, Chy = _centroid_loop(hole)
        h_area = first_moment(hole)
        Cx -= Chx/(6 * h_area) * abs(h_area)
        Cy -= Chy/(6 * h_area) * abs(h_area)
        area -= abs(h_area)
    end
    Cx /= area
    Cy /= area
    (Cx, Cy)
end

"""
    contains(polygon::Polygon, point; options...)

Inside the polygon's exterior but not in any of its holes. 
"""
function contains(polygon::Polygon, point::Point2D; on_border_is_inside::Bool=true, options...)
    hole_border = !on_border_is_inside
    contains(polygon.exterior, point; on_border_is_inside=on_border_is_inside, options...) && 
        !any(
            hole -> contains(hole, point; on_border_is_inside=hole_border, options...),
            polygon.holes
        )
end

"""
    translate(polygon::Polygon, t::Point2D)

Translate a set of polygon by delta `t`.
"""
function translate(polygon::Polygon, t::Point2D)
    # no type inference if holes is empty
    holes = isempty(polygon.holes) ? polygon.holes : map(p->translate(p, t), polygon.holes)
    Polygon(
        translate(polygon.exterior, t),
        holes,
    )
end

"""
    rotate(polygon::Polygon, θ, p0=(0, 0))

Rotate a polygon by `θ` radians about `p0`.
"""
function rotate(polygon::Polygon, θ::AbstractFloat)
    # no type inference if holes is empty
    holes = isempty(polygon.holes) ? polygon.holes : map(p->rotate(p, θ), polygon.holes)
    Polygon(
        rotate(polygon.exterior, θ),
        holes
    )
end

function rotate(polygon::Polygon, θ::AbstractFloat, p0::Point2D)
    holes = isempty(polygon.holes) ? polygon.holes : map(p->rotate(p, θ, p0), polygon.holes)
    Polygon(
        rotate(polygon.exterior, θ, p0),
        holes
    )
end

x_coords(polygon::Polygon) = 
    vcat([x_coords(polygon.exterior)], map(x_coords, polygon.holes))
y_coords(polygon::Polygon) = 
    vcat([y_coords(polygon.exterior)], map(y_coords, polygon.holes))
