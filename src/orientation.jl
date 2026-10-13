"""
    Orientation

The orientation of three ordered points `(p, q, r)`, as returned by `get_orientation`:
- `COLINEAR`: the points lie on a straight line (within tolerance).
- `CLOCKWISE`: the path `p → q → r` turns right.
- `COUNTER_CLOCKWISE`: the path `p → q → r` turns left.
"""
@enum Orientation COLINEAR=0 CLOCKWISE=1 COUNTER_CLOCKWISE=2

"""
    get_orientation(p, q, r; atol=default_atol, rtol=default_rtol)

Determine orientation of three points. 

Colinear if `|cross(pq, qr)| < atol + rtol * scale` where `scale = |pq_y * qr_x| + |qr_y * pq_x|`.

Clockwise if `cross(pq, qr)` is positive, else counter-clockwise.
"""
function get_orientation(p::Point2D, q::Point2D, r::Point2D; rtol::AbstractFloat=default_rtol, atol::AbstractFloat=default_atol)
    pq = (q[1] - p[1], q[2] - p[2])
    qr = (r[1] - q[1], r[2] - q[2])
    cross_product = pq[2] * qr[1] - qr[2] * pq[1]
    scale = abs(pq[2] * qr[1]) + abs(qr[2] * pq[1])
    tol = atol + rtol * scale
    orientation = abs(cross_product) < tol ? COLINEAR : 
        cross_product >= 0 ?  CLOCKWISE : 
        COUNTER_CLOCKWISE
    orientation
end

"""
    on_segment(point, segment, [on_line]; atol=default_atol, rtol=default_rtol)

Determine if a point lies on the segment.

The point must be colinear with the segment (see `get_orientation`) and lie within the
bounding box of the segment expanded by `atol` on all sides.

If `on_line::Bool` is given, it is used instead of calculating the colinearity.
This is useful if the orientation is already known.
"""
function on_segment(
    q::Point2D, segment::NTuple{2, Point2D}
    ; atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_rtol
    )
    on_line = get_orientation(q, segment[1], segment[2]; atol=atol, rtol=rtol) == COLINEAR
    on_segment(q, segment, on_line; atol=atol, rtol=rtol)
end

function on_segment(
    q::Point2D, segment::NTuple{2, Point2D}, on_line::Bool
    ; atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_rtol
    )
    p, r = segment
    return on_line && (
        (q[1] <= max(p[1] + atol, r[1] + atol)) &&
        (q[1] >= min(p[1] - atol, r[1] - atol)) &&
        (q[2] <= max(p[2] + atol, r[2] + atol)) &&
        (q[2] >= min(p[2] - atol, r[2] - atol))
        )
end

"""
    isless_orientation(p, q, p0; atol=default_atol, rtol=default_rtol)

Comparison function for sorting points by angle around `p0` without calculating angles.
"""
function isless_orientation(
    p::Point2D, q::Point2D, p0::Point2D
    ; atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_rtol
    )
    # a point p is "less than" another if it has a smaller angle from p0 in a counter-clockwise direction
    # or if the angle is the same, if is closer
    # instead of calculating the angle atan(p[2]-p0[2], p[1]-p0[1]), determine if (p0, p, q) is counter-clockwise
    # doesn't seem to work properly for on a circle
    ori = get_orientation(p0, p, q; atol=atol, rtol=rtol)
    if ori == COLINEAR
        dp = norm2(p, p0)
        dq = norm2(q, p0)
        val = dp < dq
    else
        val = ori == COUNTER_CLOCKWISE
    end
    val
end

"""
    isless_polar_angle(p, q, p0)

Comparison function for sorting points by their polar angle around `p0`,
from `-π` to `π` (counter-clockwise starting from the negative x-axis).
If the angles are exactly equal, the point closer to `p0` is smaller.
"""
function isless_polar_angle(p::Point2D, q::Point2D, p0::Point2D)
    angle_p = atan(p[2] - p0[2], p[1] - p0[1])
    angle_q = atan(q[2] - p0[2], q[1] - p0[1])
    if angle_p == angle_q
        dp = norm2(p, p0)
        dq = norm2(q, p0)
        val = dp < dq
    else
        val = angle_p < angle_q
    end
    val
end

"""
    sort_counter_clockwise!(points)

Sort points in place counter-clockwise by polar angle around their mean.
"""
function sort_counter_clockwise!(points::Vector{<:Point2D})
    middle = reduce(.+, points, init=(0.0, 0.0)) ./ length(points)
    sort!(points, lt=(p, q) -> isless_polar_angle(p, q, middle))
end

"""
    sort_counter_clockwise(points)

Return a copy of `points` sorted counter-clockwise. See `sort_counter_clockwise!`.
"""
sort_counter_clockwise(points::Vector{<:Point2D}) = sort_counter_clockwise!(copy(points))

"""
    sort_clockwise!(points)

Sort points in place clockwise by polar angle around their mean.
"""
function sort_clockwise!(points::Vector{<:Point2D})
    middle = reduce(.+, points, init=(0.0, 0.0)) ./ length(points)
    sort!(points, lt=(p, q) -> isless_polar_angle(q, p, middle))
end

"""
    sort_clockwise(points)

Return a copy of `points` sorted clockwise. See `sort_clockwise!`.
"""
sort_clockwise(points::Vector{<:Point2D}) = sort_clockwise!(copy(points))

"""
    cross_product(edge1, edge2)

The z-component of the cross product of the direction vectors of `edge1` and `edge2`.

Positive if `edge2` points counter-clockwise relative to `edge1`, negative if clockwise
and zero if they are parallel.
"""
function cross_product(edge1::Segment2D, edge2::Segment2D)
    a = (edge1[2][1] - edge1[1][1], edge1[2][2] - edge1[1][2])
    b = (edge2[2][1] - edge2[1][1], edge2[2][2] - edge2[1][2])
    a[1] * b[2] - a[2] * b[1]
end

"""
    dot_product(p, q, r)

The dot product of the vectors `qp` and `qr`, both starting at the vertex `q`.
"""
function dot_product(p::Point2D, q::Point2D, r::Point2D)
    (p[1] - q[1]) * (r[1] - q[1]) + (p[2] - q[2]) * (r[2] - q[2])
end

"""
    inner_angle(p, q, r)

The angle in radians at vertex `q` between the vectors `qp` and `qr`, in the range `[0, π]`.
"""
function inner_angle(p::Point2D, q::Point2D, r::Point2D)
    dot = dot_product(p, q, r)
    mag1 = dot_product(p, q, p)
    mag2 = dot_product(r, q, r)
    acos(dot / (sqrt(mag1 * mag2)))
end

"""
    in_half_plane(edge, x, [is_counter_clockwise=true]; on_border_is_inside=true)

Return `true` if the point `x` lies in the inner half-plane of `edge`.

The half-plane is drawn by extending the edge to ±infinity and back into the polygon.
For a counter-clockwise polygon this is the half-plane to the left of the edge;
for a clockwise polygon it is the half-plane to the right.

Uses exact comparisons with zero, so points very close to the line may be misclassified.
"""
function in_half_plane(edge::Segment2D, x::Point2D, is_counter_clockwise::Bool=true; on_border_is_inside=true)
    c = cross_product(edge, (edge[1], x))
    (c == 0.0 && on_border_is_inside) || (is_counter_clockwise ? c > 0 : c < 0)
end

"""
    is_vertical(segment)
    is_vertical(pt1, pt2)

Return `true` if the x co-ordinates of the two points differ by at most `eps(T)`.

This deliberately uses machine epsilon rather than `atol`.
"""
is_vertical(seg::Segment2D{T}) where T = abs(seg[2][1] - seg[1][1]) <= eps(T)
is_vertical(pt1::Point2D{T}, pt2::Point2D{T}) where T = abs(pt2[1] - pt1[1]) <= eps(T)

"""
    is_above_or_on(point, segment; atol=default_atol, rtol=default_rtol)

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
function is_above_or_on(point::Point2D{T}, segment::Segment2D{T}
    ; atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_rtol
    ) where T
    if is_vertical(segment)
        return point[2] >= max(segment[1][2], segment[2][2])
    end
    cmp = get_orientation(segment[1], segment[2], point; atol=atol, rtol=rtol)
    cmp != CLOCKWISE # is counter-clockwise or co-linear
end
