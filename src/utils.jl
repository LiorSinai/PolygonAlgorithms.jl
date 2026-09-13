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
"""
function cyclic_equality(a::AbstractVector, b::AbstractVector)
    if length(a) != length(b)
        return false
    elseif length(a) == 0
        return true
    end
    n = length(a)
    doubled_a = vcat(a, a)
    for i in 1:n
        match = true
        for j in 1:n
            if doubled_a[i + j - 1] != b[j]
                match = false
                break
            end
        end
        match && return true
    end
    false
end

function are_equivalent_collections(
    a::AbstractVector{<:AbstractVector}, b::AbstractVector{<:AbstractVector};
    match_reverse::Bool=true
    )
    if length(a) != length(b)
        return false
    elseif length(a) == 0
        return true
    end
    matched = zeros(Int, length(a))
    for (i, a_i) in enumerate(a)
        for (j, b_j) in enumerate(b)
            if j in matched
                continue
            end
            if cyclic_equality(a_i, b_j) || (match_reverse && cyclic_equality(a_i, reverse(b_j)))
                matched[i] = j
                break
            end
        end
        if matched[i] == 0
            return false
        end
    end
    true
end

"""
    are_equivalent_polygons(a, b; digits=6, match_reverse=true)

Are equivalent collections of polygons.
That is, there is a 1-1 match (cyclic equality) for each polygon in `a` with a polygon from `b`.
"""
function are_equivalent_polygons(
    a::AbstractVector{<:AbstractVector{<:T}}, b::AbstractVector{<:AbstractVector{<:T}};
    digits::Int=6, match_reverse::Bool=true
    ) where T <: Point2D
    a = map(pts -> map(pt -> round.(pt, digits=digits) .+ 0.0, pts), a)
    b = map(pts -> map(pt -> round.(pt, digits=digits) .+ 0.0, pts), b)
    are_equivalent_collections(a, b; match_reverse=match_reverse)
end
