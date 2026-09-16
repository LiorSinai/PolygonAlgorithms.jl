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

function separate(f, v::Vector)
    idxs = findall(f, v)
    other_idxs = setdiff(eachindex(v), idxs)
    (@view(v[idxs]), @view(v[other_idxs]))
end

function insert_in_order!(vec::Vector{T}, data::T; lt=isless, rev=false) where T
    idx = searchsortedfirst(vec, data; lt=lt, rev=rev)
    insert!(vec, idx, data)
end

function pop_key!(vec::Vector{T}, key::T) where T
    idx = findfirst(x->x===key, vec)
    if isnothing(idx)
        throw(KeyError(key))
    end
    popat!(vec, idx)
end

"""
    cyclic_equality(a, b)

A simple cyclic equality algorithm.

```
cyclic_equality([1, 2, 3], [2, 3, 1]) # true
cyclic_equality([1, 2, 3], [2, 1, 3]) # false
```

Warning: this function is used for testing only and it is not optimised.
This implementation is `O(n^2)`.
"""
function cyclic_equality(a::AbstractVector, b::AbstractVector)
    # TODO: KMP algorithm for O(n)
    if length(a) != length(b)
        return false
    elseif length(a) == 0
        return true
    end
    n = length(a)
    doubled_a = vcat(a, a)
    for i in 1:n
        if view(doubled_a, i:(i + n - 1)) == b
            return true
        end
    end
    false
end

"""
    minimal_rotation(array)

Find the lexicographically minimal string rotation (LMSR) in a cyclic array.

Warning: this function is used for testing only and it is not optimised.
This implementation is `O(n^2)`.

```julia
minimal_rotation([1, 4, 1, 2]) # [1, 2, 1, 4]
minimal_rotation([4, 6, 5, 1]) # [1, 4, 6, 5]
```

Reference
- https://en.wikipedia.org/wiki/Lexicographically_minimal_string_rotation
"""
function minimal_rotation(a::AbstractVector)
    # TODO: Booth's algorithm for O(n)
    n = length(a)
    n == 0 && return a
    best = a
    doubled_a = vcat(a, a)
    for i in 2:n
        candidate = view(doubled_a, i:(i+n-1))
        if candidate < best
            best = collect(candidate)
        end
    end
    best
end

function minimal_rotation_reverse(a::AbstractVector)
    forward = minimal_rotation(a)
    backward = minimal_rotation(reverse(a))
    min(forward, backward)
end

"""
    cyclic_set_equality(a, b; match_reverse=false)

Check that two sets of cyclic arrays are equivalent.

```
cyclic_set_equality(
    [[1, 2, 3], [4, 5, 6]],
    [[5, 6, 4], [3, 1, 2]]
) # true
```

Warning: this function is used for testing only and it is not optimised.
This implementation is `O(S*n^2 + S*n*log(S))` where `S` is the number of sets and `n` is the average length of the sets.
"""
function cyclic_set_equality(
    a::AbstractVector, b::AbstractVector;
    match_reverse::Bool=false
    )
    if length(a) != length(b)
        return false
    elseif length(a) == 0
        return true
    end
    fn = match_reverse ? minimal_rotation_reverse : minimal_rotation
    canonical_a = sort!(fn.(a))
    canonical_b = sort!(fn.(b))
    canonical_a == canonical_b
end

function compress_cyclic(a::AbstractVector)
    b = eltype(a)[]
    for x in a
        if isempty(b) || (x != b[end])
            push!(b, x)
        end
    end
    if length(b) > 0 && (b[end] == b[1])
        pop!(b)
    end
    b
end

function _normalise_points(pts::AbstractVector{<:Point2D{T}}; digits::Int=6) where T
    # Normalise points:
    # - round decimals
    # - convert -0.0 to 0.0
    # - remove repeat points
    map(pt -> round.(pt, digits=digits) .+ zero(T), pts) |> compress_cyclic
end

"""
    are_equivalent_polygons(polygon1, polygon2; digits=6, match_reverse=true)

Determine if the polygons are equivalent,
or alternatively, if a collections of polygons are equivalent.

Equivalence between two polygons is defined as a cyclic equality ignoring repeat points.
Equivalence between two collections is defined as a 1-1 equivalence for each polygon in `a` with a polygon from `b`.

Warning: this function is used only for testing and it is not optimised.
"""
function are_equivalent_polygons(
    a::AbstractVector{<:Point2D}, b::AbstractVector{<:Point2D}
    ; digits::Int=6, match_reverse::Bool=true
    )
    a = _normalise_points(a; digits=digits)
    b = _normalise_points(b; digits=digits)
    cyclic_equality(a, b) || (match_reverse && cyclic_equality(reverse(a), b))
end

function are_equivalent_polygons(
    a::AbstractVector{<:AbstractVector{<:Point2D}},
    b::AbstractVector{<:AbstractVector{<:Point2D}}
    ; digits::Int=6, match_reverse::Bool=true
    ) 
    a = map(pts -> _normalise_points(pts; digits=digits), a)
    b = map(pts -> _normalise_points(pts; digits=digits), b)
    cyclic_set_equality(a, b; match_reverse=match_reverse)
end
