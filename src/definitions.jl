const Point2D{T<:AbstractFloat} = NTuple{2,T}
const Segment2D{T<:AbstractFloat} = NTuple{2,Point2D{T}}
const Polygon2D{T<:AbstractFloat} = Vector{<:Point2D{T}}

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

function points_to_matrix(v::Vector{<:Point2D{T}}) where {T}
    X = Matrix{T}(undef, 2, length(v))
    for i in 1:length(v)
        X[1, i] = v[i][1]
        X[2, i] = v[i][2]
    end
    X
end

function matrix_to_points(m::AbstractMatrix{T}) where {T}
    N = size(m, 1)
    points = NTuple{N,T}[]
    for col in eachcol(m)
        push!(points, tuple(col...))
    end
    points
end

### basic operations on points
norm2(p::Point2D, q::Point2D) = (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2])
norm(p::Point2D, q::Point2D) = sqrt(norm2(p, q))

is_same_point(p::Point2D, q::Point2D; atol::AbstractFloat=1e-6) = norm(p, q) <= atol
translate(points::Vector{<:Point2D}, t::Point2D) = [(p[1] + t[1], p[2] + t[2]) for p in points]

"""
    rotate(points, θ, p0=(0, 0))

Rotate a set of points `θ` radians about `p0`.
"""
function rotate(points::Vector{<:Point2D}, θ::Float64)
    s = sin(θ)
    c = cos(θ)
    [(p[1] * c - p[2] * s, p[1] * s + p[2] * c) for p in points]
end

rotate(points::Vector{<:Point2D}, θ::Float64, p0::Point2D) = 
    translate(rotate(translate(points, (-p0[1], -p0[2])), θ), p0)

x_coords(points::Polygon2D) = [p[1] for p in points]
x_coords(points::Polygon2D, idxs::AbstractVector{Int}) = [p[1] for p in points[idxs]]
y_coords(points::Polygon2D) = [p[2] for p in points]
y_coords(points::Polygon2D, idxs::AbstractVector{Int}) = [p[2] for p in points[idxs]]

function separate(f, v::Vector)
    idxs = findall(f, v)
    other_idxs = setdiff(eachindex(v), idxs)
    (@view(v[idxs]), @view(v[other_idxs]))
end
  