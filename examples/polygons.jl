using PolygonAlgorithms: Point2D

function make_circle(
    n::Int=1000
    ;
    radius::AbstractFloat=1.0,
    origin::Point2D=(0.0, 0.0),
    )
    angles = range(0.0, 2π; length=n)
    x = radius * cos.(angles)
    y = radius * sin.(angles)
    collect(zip(x .+ origin[1], y .+ origin[2]))
end

function hilbert_curve(
    point::Tuple{T, T},
    horiz::Tuple{T, T},
    vert::Tuple{T, T},
    level::Int
    ) where T <: AbstractFloat
    # Algorithm from http://www.fundza.com/algorithmic/space_filling/hilbert/basics/
    if level <= 0
        next = ((point[1] + (horiz[1] + vert[1])/2), (point[2] + (horiz[2] + vert[2])/2))
        return [next]
    end
    points = [
        hilbert_curve(point, (horiz[2]/2, vert[2]/2), (horiz[1]/2, vert[1]/2), level-1)...,
        hilbert_curve((point[1] + horiz[1]/2, point[2] + horiz[2]/2), (horiz[1]/2, horiz[2]/2), (vert[1]/2, vert[2]/2), level-1)...,
        hilbert_curve((point[1] + horiz[1]/2+vert[1]/2, point[2] + horiz[2]/2+vert[2]/2), (horiz[1]/2, horiz[2]/2), (vert[1]/2, vert[2]/2), level-1)...,
        hilbert_curve((point[1] + horiz[1]/2+vert[1], point[2] + horiz[2]/2+vert[2]), (-vert[1]/2, -vert[2]/2), (-horiz[1]/2, -horiz[2]/2), level-1)...
    ]
    points
end

function make_hilbert_polygon(order::Int=6)
    points = hilbert_curve((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), order);
    tail = points[end]
    head = points[1]
    ## fill inner
    push!(points, (-tail[1], tail[2]))
    push!(points, (-head[1], head[2]))
    reverse(points)
end

function make_forks(n::Int=3
    ; origin::Point2D=(0.0, 0.0),
    step_size::AbstractFloat=1.0,
    height::AbstractFloat=4.0,
    base::AbstractFloat=1.0,
    )
    prev = origin
    points = [prev]
    for i in 1:(n-1)
        p1 = (prev[1], height)
        p2 = (prev[1] + step_size, height)
        p3 = (prev[1] + step_size, base)
        p4 = (prev[1] + 2 * step_size, base)
        push!(points, p1, p2, p3, p4)
        prev = p4
    end
    p1 = (prev[1], height)
    p2 = (prev[1] + step_size,  height)
    push!(points, p1, p2, (prev[1] + step_size, origin[2]))
    points
end

function random_spiked_balls(rng, n::Int=100, origin=(0.0, 0.0))
    points = [(origin[1] + rand(rng), origin[2] + rand(rng)) for i in 1:n]
    PolygonAlgorithms.sort_clockwise(points)
end

function random_concave_saw_from_points(points::Vector{<:Point2D})
    order = sort(points)
    left = order[1]
    right = order[end]
    divider = (left, right)
    indicators = [PolygonAlgorithms.is_above_or_on(pt, divider) for pt in order]
    above = order[indicators]
    below = reverse!(order[.!indicators])
    vcat(above, below)
end

function random_self_intersecting(rng, n::Int=100, origin=(0.0, 0.0))
     [(origin[1] + rand(rng), origin[2] + rand(rng)) for i in 1:n]
end

"""
    star_polygon(p, q)

Star polygon `{p/q}`: the `p` vertices of a regular polygon on the unit circle, joined by
connecting every `q`th vertex. For `q > 1` the paths are self-intersecting.

Returns a vector of paths. If `d = gcd(p, q) > 1` the star is a compound of `d` separate
polygons, each rotated by `2π/p` from the last; otherwise the vector has a single path.
"""
function star_polygon(p::Int, q::Int)
    d = gcd(p, q)
    n = p ÷ d # vertices per polygon
    if p < 3 || !(1 <= q < p) || n < 3
        throw(ArgumentError("star_polygon requires p >= 3, 1 <= q < p and p/gcd(p, q) >= 3; got p=$p, q=$q"))
    end
    # start at the top (offset by j for each polygon), step q vertices at a time
    [
        [sincospi(1/2 + 2 * mod(j + k * q, p) / p) |> reverse for k in 0:(n - 1)]
        for j in 0:(d - 1)
    ]
end
