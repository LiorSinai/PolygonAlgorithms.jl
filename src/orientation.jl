"""
    Orientation

The orientation of three ordered points `(p, q, r)`, as returned by `get_orientation`:
- `COLINEAR`: the points lie on a straight line (within tolerance).
- `CLOCKWISE`: the path `p → q → r` turns right.
- `COUNTER_CLOCKWISE`: the path `p → q → r` turns left.
"""
@enum Orientation COLINEAR=0 CLOCKWISE=1 COUNTER_CLOCKWISE=2

"""
    get_orientation(p, q, r; rtol=default_rtol(T))

Determine orientation of three points.

The cross product `cross(pq, qr) = L - R` is the difference of two products.
The points are colinear if this is within floating point error:
```
|L - R| ≤ rtol * (|L| + |R| + scale * side)
```
where `scale` is the largest absolute co-ordinate and `side` is the longest side of the triangle `pqr`,
measured with the Manhattan distance `|Δx| + |Δy|` to avoid square roots.
Otherwise counter-clockwise if `cross(pq, qr)` is positive (see `cross_product`), else clockwise.

Tolerance derivation:
1. Arithmetic error. Each floating point operation has a relative error of at most `u = eps(T)/2`.
Shewchuk showed that the computed cross product differs from the exact value by at most
`(3u + 16u²) * (|L| + |R|) ≈ 1.5eps(T) * (|L| + |R|)`.
The bound is relative to `|L| + |R|`, not `|L - R|`, so it is meaningful for nearly colinear points,
where `L` and `R` are large and almost cancel.
2. Representation error. Shewchuk's bound assumes the points are exact.
But points calculated from intersections are rounded, so each co-ordinate is only known to within
about `u * scale`. Moving one point a distance `δ` changes the cross product by at most `δ * side`.
This is important when the co-ordinates are large compared to the triangle,
for example a short segment far from the origin.

The test is symmetric: all permutations of `(p, q, r)` agree on colinearity, and
the orientation flips if two points are swapped.

Reference: "Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates",
Jonathan Richard Shewchuk (1997).
"""
function get_orientation(p::Point2D{T}, q::Point2D{T}, r::Point2D{T}; rtol::AbstractFloat=default_rtol(T)) where T
    pq = (q[1] - p[1], q[2] - p[2])
    qr = (r[1] - q[1], r[2] - q[2])
    L = pq[2] * qr[1]
    R = qr[2] * pq[1]
    cross_product = L - R
    scale = max(abs(p[1]), abs(p[2]), abs(q[1]), abs(q[2]), abs(r[1]), abs(r[2]))
    side = max(
        abs(pq[1]) + abs(pq[2]),
        abs(qr[1]) + abs(qr[2]),
        abs(r[1] - p[1]) + abs(r[2] - p[2]),
    )
    is_colinear = abs(cross_product) <= rtol * (abs(L) + abs(R) + scale * side)
    orientation = is_colinear ? COLINEAR :
        cross_product >= 0 ?  CLOCKWISE :
        COUNTER_CLOCKWISE
    orientation
end

"""
    on_segment(point, segment, [on_line]; rtol=default_rtol(T))

Determine if a point lies on the segment.

The point must be colinear with the segment (see `get_orientation`) and lie within the
bounding box of the segment, expanded by the floating point error of the co-ordinates:
`rtol * max(|x|)` horizontally and `rtol * max(|y|)` vertically.

If `on_line::Bool` is given, `false` returns early without calculating the distance.
This is useful if the orientation is already known.
"""
function on_segment(
    q::Point2D{T}, segment::NTuple{2, Point2D{T}}
    ; rtol::AbstractFloat=default_rtol(T)
    ) where T
    # the bounding box check is cheaper than the orientation, so do it first
    _in_expanded_bounds(q, segment, rtol) &&
        get_orientation(segment[1], segment[2], q; rtol=rtol) == COLINEAR
end

function on_segment(
    q::Point2D{T}, segment::NTuple{2, Point2D{T}}, on_line::Bool
    ; rtol::AbstractFloat=default_rtol(T)
    ) where T
    on_line && _in_expanded_bounds(q, segment, rtol)
end

function _in_distance(p::Point2D, segment::NTuple{2, Point2D}, atol::AbstractFloat)
    # This places a rectangular bounds of half-width atol around the segment,
    # with two half-circles at either end
    # s = s1 + t * (s2 - s1)  # equation for the segment
    # dot(p - s, s2 - s1) = 0 # shortest at perpendicular
    # ⟹ t = dot(p - s1, s2 - s1) / |s2 - s1|², then clamp t to [0, 1]
    dx = segment[2][1] - segment[1][1]
    dy = segment[2][2] - segment[1][2]
    length2 = dx * dx + dy * dy
    if iszero(length2) # degenerate segment, avoid 0/0
        return norm2(p, segment[1]) <= atol * atol
    end
    t = ((p[1] - segment[1][1]) * dx + (p[2] - segment[1][2]) * dy) / length2
    t = clamp(t, zero(t), one(t))
    # on segment if |p - s| ≤ atol
    px = segment[1][1] + t * dx - p[1]
    py = segment[1][2] + t * dy - p[2]
    px * px + py * py <= atol * atol 
end

function _in_expanded_bounds(q::Point2D, segment::NTuple{2, Point2D}, rtol::AbstractFloat)
    # This places an axis aligned rectangular bounds around the segment,
    # expanded by the floating point error of the co-ordinates.
    p, r = segment
    δx = rtol * max(abs(p[1]), abs(r[1]), abs(q[1]))
    δy = rtol * max(abs(p[2]), abs(r[2]), abs(q[2]))
    (q[1] <= max(p[1], r[1]) + δx) &&
    (q[1] >= min(p[1], r[1]) - δx) &&
    (q[2] <= max(p[2], r[2]) + δy) &&
    (q[2] >= min(p[2], r[2]) - δy)
end

"""
    isless_orientation(p, q, p0; rtol=default_rtol(T))

Comparison function for sorting points by angle around `p0` without calculating angles.
"""
function isless_orientation(
    p::Point2D{T}, q::Point2D{T}, p0::Point2D{T}
    ; rtol::AbstractFloat=default_rtol(T)
    ) where T
    # a point p is "less than" another if it has a smaller angle from p0 in a counter-clockwise direction
    # or if the angle is the same, if is closer
    # instead of calculating the angle atan(p[2]-p0[2], p[1]-p0[1]), determine if (p0, p, q) is counter-clockwise
    # doesn't seem to work properly for on a circle
    ori = get_orientation(p0, p, q; rtol=rtol)
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
    cross_product(p, q, r)

The z-component of the cross product of the direction vectors of `edge1` and `edge2`.
Or the cross product of `pq` and `qr`, which is the same as `cross_product((p, q), (q, r))`.

Positive if `edge2` points counter-clockwise relative to `edge1`, negative if clockwise
and zero if they are parallel.
For three points, positive if `p → q → r` turns left (counter-clockwise).
"""
function cross_product(edge1::Segment2D, edge2::Segment2D)
    a = (edge1[2][1] - edge1[1][1], edge1[2][2] - edge1[1][2])
    b = (edge2[2][1] - edge2[1][1], edge2[2][2] - edge2[1][2])
    a[1] * b[2] - a[2] * b[1]
end

function cross_product(p::Point2D, q::Point2D, r::Point2D)
    pq = (q[1] - p[1], q[2] - p[2])
    qr = (r[1] - q[1], r[2] - q[2])
    pq[1] * qr[2] - pq[2] * qr[1]
end

"""
    dot_product(p, q, r)

The dot product of the vectors `qp` and `qr`, both starting at the vertex `q`.
"""
function dot_product(p::Point2D, q::Point2D, r::Point2D)
    (p[1] - q[1]) * (r[1] - q[1]) + (p[2] - q[2]) * (r[2] - q[2])
end

function dot_product(edge1::Segment2D, edge2::Segment2D)
    a = (edge1[2][1] - edge1[1][1], edge1[2][2] - edge1[1][2])
    b = (edge2[2][1] - edge2[1][1], edge2[2][2] - edge2[1][2])
    a[1] * b[1] + a[2] * b[2]
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
    is_vertical(segment; rtol=default_rtol(T))
    is_vertical(pt1, pt2; rtol=default_rtol(T))

Return `true` if the x co-ordinates of the two points are equal within floating point error:
`|x₂ - x₁| ≤ rtol * max(|x₁|, |x₂|)`.
"""
is_vertical(seg::Segment2D{T}; rtol::AbstractFloat=default_rtol(T)) where T = is_vertical(seg[1], seg[2]; rtol=rtol)
function is_vertical(pt1::Point2D{T}, pt2::Point2D{T}; rtol::AbstractFloat=default_rtol(T)) where T
    abs(pt2[1] - pt1[1]) <= rtol * max(abs(pt1[1]), abs(pt2[1]))
end

"""
    is_above_or_on(point, segment; rtol=default_rtol(T))

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
    ; rtol::AbstractFloat=default_rtol(T)
    ) where T
    if is_vertical(segment; rtol=rtol)
        return point[2] >= max(segment[1][2], segment[2][2])
    end
    cmp = get_orientation(segment[1], segment[2], point; rtol=rtol)
    cmp != CLOCKWISE # is counter-clockwise or co-linear
end
