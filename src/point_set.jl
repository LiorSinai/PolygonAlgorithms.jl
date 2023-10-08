import Base: length, show, isempty, iterate, KeySet, union!, push!, show_vector, in, Set

"""
    PointSet([itr]; digits=6)

A specialised set based off `Base.Set`. It is used to compensate for minor numeric errors.
Each co-ordinate of the point is rounded to `digits` decimal places before being inserted in the set.
Likewise, the `in` operation rounds points before checking if they are in the set.
It also treats `-0.0` and `0.0` as equivalent unlike `Set([-0.0]) != Set([0.0])`.

For example:
```
s = PointSet([(1e-7, 0.0), (1.0, π)])
(0.0, 0.0) in s # true
(1.0, 3.141593) in s # true
(-0.0, -0.0) in s # true
```
"""
struct PointSet{T} <: AbstractSet{T}
    dict::Dict{Tuple{T, T},Nothing}
    digits::Int
    function PointSet(dict::Dict{Tuple{T, T},Nothing}, digits) where T
        point_set = new{T}(dict, digits)
        return point_set
    end
end

PointSet(;digits::Int=6) = PointSet(Dict{Tuple{Float64, Float64},Nothing}(), digits)
PointSet(itr::AbstractVector{<:Tuple}; digits::Int=6) = union!(PointSet(digits=digits), itr)
isempty(s::PointSet) = isempty(s.dict)
length(s::PointSet)  = length(s.dict)
iterate(s::PointSet, i...)  = iterate(KeySet(s.dict), i...)
in(x, s::PointSet) = haskey(s.dict, fudge(x, s.digits))

Set(s::PointSet) = Set(KeySet(s.dict))
PointSet(s::Set{T}) where T = PointSet(collect(KeySet(s.dict)))

function union!(s::PointSet, itr)
    for x in itr
        push!(s, x)
    end
    return s
end

function push!(s::PointSet, x::Tuple)
    s.dict[fudge(x, s.digits)] = nothing
    s
end

function fudge(point::Tuple, digits::Int)
    point = round.(point, digits=digits)
    if (point[1] == -0 || point[2] == -0)
        x, y = point
        point = (x == -0 ? zero(x) : x, y == -0 ? zero(y) : y)
    end
    point
end

function show(io::IO, s::PointSet)
    if isempty(s)
        show(io, typeof(s))
        print(io, "(digits=$(s.digits),)")
    else
        print(io, typeof(s), "(")
        print(io, "digits=$(s.digits),")
        show_vector(io, s)
        print(io, ')')
    end
end