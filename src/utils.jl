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
