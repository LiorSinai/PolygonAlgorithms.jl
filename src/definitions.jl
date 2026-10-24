"""
    Point2D{T} = NTuple{2,T}

A point `(x, y)` in the plane.
"""
const Point2D{T<:AbstractFloat} = NTuple{2,T}

"""
    Segment2D{T} = NTuple{2,Point2D{T}}

A line segment between two points `(p1, p2)`.
"""
const Segment2D{T<:AbstractFloat} = NTuple{2,Point2D{T}}

"""
    Path2D{T} = Vector{<:Point2D{T}}

A polygon represented as a list of vertices.
The last point implicitly connects back to the first point, so it should not be repeated.
"""
const Path2D{T<:AbstractFloat} = Vector{<:Point2D{T}}
"""
    default_rtol(T)

The default tolerance for floats of type `T`: `64eps(T)`.
This is `1.4e-14` for `Float64` and `7.6e-6` for `Float32`.

`rtol` is a relative tolerance for floating point error only.
A computed value `x` is treated as zero if `|x| ≤ rtol * Σ|terms|`, where the terms are the products
that were added or subtracted to calculate `x`.
For example, Shewchuk's error bound for the orientation determinant is `≈ 1.5eps * Σ|terms|`.

The default is larger than this bound to allow for inputs which were themselves calculated,
for example by rotation. Their rounding error is relative to the size of the numbers in that calculation,
which can be larger than the co-ordinates of the point.
"""
default_rtol(::Type{T}) where T<:AbstractFloat = 64eps(T)

"""
    decimal_tolerance(rtol)

The number of decimal places to round based on the tolerance.
"""
decimal_tolerance(rtol::AbstractFloat) = rtol > 1 ? 0 : ceil(Int, min(18, -log10(rtol)))

"""
    Line2D(a, b, c)

Equation of line is `ax + by = c`.
"""
struct Line2D{T<:AbstractFloat}
    a::T
    b::T
    c::T
end

function Base.getindex(line::Line2D, ind::Int)
    if ind > 3 || ind < 1
        throw(DomainError(Line2D, "Attempted to access 3-element Line2D at ind=$ind"))
    end
    ind == 1 ? line.a : ind == 2 ? line.b : line.c
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

### basic operations on points
norm2(p::Point2D, q::Point2D) = (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2])
norm(p::Point2D, q::Point2D) = sqrt(norm2(p, q))

"""
    is_same_point(p, q; rtol=default_rtol(T))

Return `true` if the distance between `p` and `q` is within floating point error of their magnitudes:
`|p - q|² ≤ rtol² * max(|p|², |q|²)`.
"""
function is_same_point(p::Point2D{T}, q::Point2D{T}; rtol::AbstractFloat=default_rtol(T)) where T
    scale2 = max(p[1] * p[1] + p[2] * p[2], q[1] * q[1] + q[2] * q[2])
    norm2(p, q) <= rtol * rtol * scale2
end

"""
    translate(points::Vector{<:Point2D}, t::Point2D)

Translate a set of points by delta `t`.
"""
translate(points::Vector{<:Point2D}, t::Point2D) = [(p[1] + t[1], p[2] + t[2]) for p in points]

"""
    rotate(points, θ, p0=(0, 0))

Rotate a set of points `θ` radians about `p0`.
"""
function rotate(points::Vector{<:Point2D}, θ::AbstractFloat)
    s = sin(θ)
    c = cos(θ)
    [(p[1] * c - p[2] * s, p[1] * s + p[2] * c) for p in points]
end

rotate(points::Vector{<:Point2D}, θ::AbstractFloat, p0::Point2D) = 
    translate(rotate(translate(points, (-p0[1], -p0[2])), θ), p0)

_coords_matrix(points::Path2D) = reinterpret(reshape, eltype(eltype(points)), points)

"""
    x_coords(points, [idxs])
    x_coords(polygon::Polygon)

The x co-ordinates of `points`, optionally only at the indices `idxs`.

For a `Polygon`, returns a vector of x co-ordinates for the exterior followed by each hole.
"""
x_coords(points::Path2D) = view(_coords_matrix(points), 1, :)
x_coords(points::Path2D, idxs::AbstractVector{Int}) = view(_coords_matrix(points), 1, idxs)
"""
    y_coords(points, [idxs])
    y_coords(polygon::Polygon)

The y co-ordinates of `points`, optionally only at the indices `idxs`.

For a `Polygon`, returns a vector of y co-ordinates for the exterior followed by each hole.
"""
y_coords(points::Path2D) = view(_coords_matrix(points), 2, :)
y_coords(points::Path2D, idxs::AbstractVector{Int}) = view(_coords_matrix(points), 2, idxs)
