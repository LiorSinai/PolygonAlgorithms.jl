@enum Orientation COLINEAR=0 CLOCKWISE=1 COUNTER_CLOCKWISE=2

"""
    get_orientation(p, q, r; rtol=1e-4, atol=1e-6)

Determine orientation of three points. 

Colinear is returned if `cross(pq, qr) <= tol`, where `tol` is the 
- `atol` if either magnitude is less than `atol`.
- `rtol * |pq||qr|` otherwise, where `cross(pq, qr) ≈ |pq||qr|θ` for small `θ`.

Clockwise is returned if `cross(pq, qr)` is positive, else counter-clockwise.
"""
function get_orientation(p::Point2D, q::Point2D, r::Point2D; rtol::AbstractFloat=1e-4, atol::AbstractFloat=1e-6)
    pq = (q[1] - p[1], q[2] - p[2])
    qr = (r[1] - q[1], r[2] - q[2])
    cross_product = pq[2] * qr[1] - qr[2] * pq[1]
    mag_pq = sqrt(pq[1] * pq[1] + pq[2] * pq[2])
    mag_qr = sqrt(qr[1] * qr[1] + qr[2] * qr[2])
    tol = (mag_qr >= atol && mag_pq >= atol) ? (mag_pq * mag_qr) * rtol : atol
    orientation = abs(cross_product) <= tol ? COLINEAR : 
        cross_product >= 0 ?  CLOCKWISE : 
        COUNTER_CLOCKWISE
    orientation
end

"""
    on_segment(point, segment, [on_line]; atol=1e-6)

Determine if a point lies on the segment. 
"""
function on_segment(q::Point2D, segment::NTuple{2, Point2D}; atol::AbstractFloat=1e-6)
    on_line = get_orientation(q, segment[1], segment[2]) == COLINEAR
    on_segment(q, segment, on_line; atol=atol)
end

function on_segment(
    q::Point2D, segment::NTuple{2, Point2D}, on_line::Bool
    ; atol::AbstractFloat=1e-6
    )
    p, r = segment
    return on_line && (
        (q[1] <= max(p[1] + atol, r[1] + atol)) &&
        (q[1] >= min(p[1] - atol, r[1] - atol)) &&
        (q[2] <= max(p[2] + atol, r[2] + atol)) &&
        (q[2] >= min(p[2] - atol, r[2] - atol))
        )
end

function isless_orientation(p::Point2D, q::Point2D, p0::Point2D; rtol::AbstractFloat=1e-4)
    # a point p is "less than" another if it has a smaller angle from p0 in a counter-clockwise direction
    # or if the angle is the same, if is closer
    # instead of calculating the angle atan(p[2]-p0[2], p[1]-p0[1]), determine if (p0, p, q) is counter-clockwise
    # doesn't seem to work properly for on a circle
    ori = get_orientation(p0, p, q; rtol=rtol)
    if ori == COLINEAR
        dp = norm2(p, p0)
        dq = norm2(q, p0)
        val = dp <= dq
    else
        val = ori == COUNTER_CLOCKWISE
    end
    val
end

function isless_polar_angle(p::Point2D, q::Point2D, p0::Point2D)
    angle_p = atan(p[2] - p0[2], p[1] - p0[1])
    angle_q = atan(q[2] - p0[2], q[1] - p0[1])
    if angle_p == angle_q
        dp = norm2(p, p0)
        dq = norm2(q, p0)
        val = dp <= dq
    else
        val = angle_p < angle_q
    end
    val
end

function sort_counter_clockwise(points::Vector{<:Point2D})
    middle = reduce(.+, points, init=(0.0, 0.0)) ./ length(points)
    sort(points, lt=(p, q) -> isless_polar_angle(p, q, middle))
end

function sort_clockwise(points::Vector{<:Point2D})
    middle = reduce(.+, points, init=(0.0, 0.0)) ./ length(points)
    sort(points, lt=(p, q) -> !isless_polar_angle(p, q, middle))
end

# orientation between two lines
function cross_product(edge1::Segment2D, edge2::Segment2D)
    a = (edge1[2][1] - edge1[1][1], edge1[2][2] - edge1[1][2])
    b = (edge2[2][1] - edge2[1][1], edge2[2][2] - edge2[1][2])
    a[1] * b[2] - a[2] * b[1]
end

# angle between two lines
function dot_product(p::Point2D, q::Point2D, r::Point2D)
    (p[1] - q[1]) * (r[1] - q[1]) + (p[2] - q[2]) * (r[2] - q[2])
end

function inner_angle(p::Point2D, q::Point2D, r::Point2D)
    dot = dot_product(p, q, r)
    mag1 = dot_product(p, q, p)
    mag2 = dot_product(r, q, r)
    acos(dot / (sqrt(mag1 * mag2)))
end

# the half-plane is drawn by extending the edge to ±infinity and back into the polygon
function in_half_plane(edge::Segment2D, x::Point2D, is_counter_clockwise::Bool=true; on_border_is_inside=true)
    c = cross_product(edge, (edge[1], x))
    (c == 0.0 && on_border_is_inside) || (is_counter_clockwise ? c > 0 : c < 0)
end

"""
    is_above_or_on(point, segment)

Is `point` above or on `segment`?

In the general case, checks if the point `(xp, yp)` is above or on the line through the segment
`((x₁, y₁), (x₂, y₂))`:
```
yp ≥ (y₂-y₁)/(x₂-x₁) * (xp-x₁) + y₁
```

This is equivalent to checking if the three points `((x₁, y₁), (x₂, y₂), (xp, yp))` are orientated counter-clockwise or are co-linear:
```
0  ≥ (y₂-y₁) * (xp-x₁) - (yp-y₁) * (x₂-x₁)
```

In the special case of a vertical segment (`x₂=x₁`), this compares `y` values:
```
yp ≥ max(y₂, y₁)
```
"""
function is_above_or_on(point::Point2D, segment::Segment2D; atol::AbstractFloat=1e-6, rtol::AbstractFloat=1e-4)
    if abs(segment[2][1] - segment[1][1]) <= atol # vertical segment
        return point[2] >= max(segment[1][2], segment[2][2])
    end
    cmp = get_orientation(segment[1], segment[2], point; atol=atol, rtol=rtol)
    cmp != CLOCKWISE # is counter-clockwise or co-linear
end
