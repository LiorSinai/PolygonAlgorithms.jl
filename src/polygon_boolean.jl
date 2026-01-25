include("boolean/weiler_atherton.jl")
include("boolean/chasing_edges.jl")
include("boolean/martinez_rueda.jl")

abstract type PolygonIntersectionAlgorithm end

"""
    PointSearchAlg()
- Time complexity: `O(nm)` where `n` and `m` vertices on polygon 1 and 2 respectively.
- Should only be used for convex polygons with a single region of intersection.
    Not guaranteed to give correct results.
- Description:
    1. Intersect all edge pairs.
    2. Selected all points in the other polygon.
    3. Sort results counter-clockwise. 
- For general non-intersecting polygons, the intersection points are valid but the order is not.
"""
struct PointSearchAlg <: PolygonIntersectionAlgorithm end

"""
    ChasingEdgesAlg():
- Time complexity: `O(n+m)` where `n` and `m` vertices on polygon 1 and 2 respectively.
- Only gives correct results for convex polygons with one region of intersection.
- Description: this algorithm rotates two pointers, one around each polygon. Each iteration it only
    advances one pointer based on a set of "advance rules". After two cycles around the polygons 
    it is guaranteed to have found all (zero, one or both) intersection points and 
    all the points in between.
- Algorithm is from "A New Linear Algorithm for Intersecting Convex Polygons" (1981) by Joseph O'Rourke et. al.
- Reference: https://www.cs.jhu.edu/~misha/Spring16/ORourke82.pdf
"""
struct ChasingEdgesAlg <: PolygonIntersectionAlgorithm end

"""
    WeilerAthertonAlg()
- Time complexity: `O(nm)` where `n` and `m` vertices on polygon 1 and 2 respectively.
- Returns regions, edges and single points of intersection. Only returns the larger type if one is within another e.g. an edge is also part of a region.
- Description: operates at a point level. Walks from point to point along `polygon2`. 
    It starts recording loops at "entry points" - crossings from `polygon2` to `polygon1` - and 
    stops recording when it gets back to the same entry point. So it continues walking along `polygon2`
    until it reaches the start point.
- Limitations
    1. This version does not cater for holes.
    2. It can fail completely for self-intersecting polygons.
"""
struct WeilerAthertonAlg <: PolygonIntersectionAlgorithm end

"""
    MartinezRuedaAlg()

- Time complexity: `O((n+m+k)log(n+m))` where `n` and `m` vertices on polygon 1 and 2 respectively
    and `k` is the total number of intersections between all segments.
- Returns regions and edges of interest. Does not return single points.
- Works for convex and concave polygons including with holes and self-intersections.
- Description: operates at a segment level and is an extension of the Bentley-Ottman line intersection algorithm.
    - Segments are scanned from left to right, bottom to top. 
    - The key assumption is that only the segments immediately above and below the current segment need to be inspected for intersections.
    This makes the algorithm fast but also sensitive to determining these segments correctly.
    - The segment that is immediately below (or empty space) is used to determine the fill annotations for the current segment.
    - Once all annotations are done, the desired segments can be selected that match a given criteria.
- Limitations:
    1. It can fail for improper polygons: polygons with lines sticking out.
    2. It is sensitive to numeric inaccuracies e.g. a line that is almost vertical or 
    tiny regions of intersection.
"""
struct MartinezRuedaAlg <: PolygonIntersectionAlgorithm end

# General

"""
    intersect_geometry([alg=MartinezRuedaAlg()], subject, others...; atol=default_atol, rtol=default_rtol)
    intersect_geometry(WeilerAthertonAlg(), subject, other; atol=default_atol) # only two arguments supported for this algorithm.

General case of polygon intersection.

The polygons can be of type:
- `Vector{<:Tuple{T, T}} where T<:AbstractFloat`. Output will be `Vector{Vector{Tuple{T, T}}}`. Possibly returns holes but does not classify them as holes.
- `Polygon{T<:AbstractFloat}`. Output will be `Vector{Polygon{T}}`. This will classify holes and match them to the smallest possible parent.
- The subject can also be the same type as the output. This enables the output to be an input. 
    Note that in this case the `subject` list is treated as a single polygon.
    If any of the polygons overlap, this is equivalent to passing a self-intersecting polygon and
    some areas might be classified as holes according to the even-odd rule.

Using multiple arguments is more efficient than looping through the polygons and applying `intersect_geometry` sequentially.

Use `intersect_convex` for convex polygons for an `O(n+m)` algorithm.
"""
intersect_geometry(polygon1::Path2D, others::Vararg{Path2D}; options...) = 
    intersect_geometry(MartinezRuedaAlg(), polygon1, others...; options...)

intersect_geometry(polygon1::Polygon, others::Vararg{Polygon}; options...) = 
    intersect_geometry(MartinezRuedaAlg(), polygon1, others...; options...)

intersect_geometry(polygon1::AbstractVector{<:Path2D}, others::Vararg{Path2D}; options...) = 
    intersect_geometry(MartinezRuedaAlg(), polygon1, others...; options...)

intersect_geometry(polygon1::AbstractVector{<:Polygon}, others::Vararg{Polygon}; options...) = 
    intersect_geometry(MartinezRuedaAlg(), polygon1, others...; options...)

function intersect_geometry(
    alg::WeilerAthertonAlg, polygon1::Path2D, polygon2::Path2D, 
    ; options...
    )
    weiler_atherton_algorithm(polygon1, polygon2; options...)
end

function intersect_geometry(
    alg::MartinezRuedaAlg, polygon1::Path2D, others::Vararg{Path2D}
    ; options...
    )
    martinez_rueda_algorithm(INTERSECTION_CRITERIA, polygon1, others...; options...)
end

function intersect_geometry(
    alg::MartinezRuedaAlg, polygon1::Polygon, others::Vararg{Polygon}
    ; options...
    )
    martinez_rueda_algorithm(INTERSECTION_CRITERIA, polygon1, others...; options...)
end

function intersect_geometry(
    alg::MartinezRuedaAlg, subjects::AbstractVector{<:Path2D}, others::Vararg{Path2D},
    ; options...
    )
    martinez_rueda_algorithm(INTERSECTION_CRITERIA, subjects, others...; options...)
end

function intersect_geometry(
    alg::MartinezRuedaAlg, subjects::AbstractVector{<:Polygon}, others::Vararg{Polygon}
    ; options...
    )
    martinez_rueda_algorithm(INTERSECTION_CRITERIA, subjects, others...; options...)
end


"""
    difference_geometry([alg=MartinezRuedaAlg()], subject, clips...; atol=default_atol, rtol=default_rtol)

General case of polygon difference: points in `subject` that are not in `clip(s)`.

The polygons can be of type:
- `Vector{<:Tuple{T, T}} where T<:AbstractFloat`. Output will be `Vector{Vector{Tuple{T, T}}}`. Possibly returns holes but does not classify them as holes.
- `Polygon{T<:AbstractFloat}`. Output will be `Vector{Polygon{T}}`. This will classify holes and match them to the smallest possible parent.
- The subject can also be the same type as the output. This enables the output to be an input. 
    Note that in this case the `subject` list is treated as a single polygon.
    If any of the polygons overlap, this is equivalent to passing a self-intersecting polygon and
    some areas might be classified as holes according to the even-odd rule.

Using multiple arguments is more efficient than looping through the polygons and applying `difference_geometry` sequentially.

`alg` can only be `MartinezRuedaAlg()`.
"""
difference_geometry(polygon1::Path2D, clips::Vararg{Path2D}; options...) = 
    difference_geometry(MartinezRuedaAlg(), polygon1, clips...; options...)

difference_geometry(polygon1::Polygon, clips::Vararg{Polygon}; options...) = 
    difference_geometry(MartinezRuedaAlg(), polygon1, clips...; options...)

difference_geometry(polygon1::AbstractVector{<:Path2D}, clips::Vararg{Path2D}; options...) = 
    difference_geometry(MartinezRuedaAlg(), polygon1, clips...; options...)

difference_geometry(polygon1::AbstractVector{<:Polygon}, clips::Vararg{Polygon}; options...) = 
    difference_geometry(MartinezRuedaAlg(), polygon1, clips...; options...)

function difference_geometry(
    alg::MartinezRuedaAlg, subject::Path2D, clips::Vararg{Path2D}
    ; options...
    )
    martinez_rueda_algorithm(DIFFERENCE_CRITERIA, subject, clips...; options...)
end

function difference_geometry(
    alg::MartinezRuedaAlg, subject::Polygon, clips::Vararg{Polygon}; options...
    )
    martinez_rueda_algorithm(DIFFERENCE_CRITERIA, subject, clips...; options...)
end

function difference_geometry(
    alg::MartinezRuedaAlg, subjects::AbstractVector{<:Path2D}, clips::Vararg{Path2D}; options...
    )
    martinez_rueda_algorithm(DIFFERENCE_CRITERIA, subjects, clips...; options...)
end

function difference_geometry(
    alg::MartinezRuedaAlg, subjects::AbstractVector{<:Polygon}, clips::Vararg{Polygon}; options...
    )
    martinez_rueda_algorithm(DIFFERENCE_CRITERIA, subjects, clips...; options...)
end

"""
    union_geometry([alg=MartinezRuedaAlg()], subject, others...; atol=default_atol, rtol=default_rtol)

General case of polygon union: in all polygons. 

The polygons can be of type:
- `Vector{<:Tuple{T, T}} where T<:AbstractFloat`. Output will be `Vector{Vector{Tuple{T, T}}}`. Possibly returns holes but does not classify them as holes.
- `Polygon{T<:AbstractFloat}`. Output will be `Vector{Polygon{T}}`. This will classify holes and match them to the smallest possible parent.
- The subject can also be the same type as the output. This enables the output to be an input. 
    Note that in this case the `subject` list is treated as a single polygon.
    If any of the polygons overlap, this is equivalent to passing a self-intersecting polygon and
    some areas might be classified as holes according to the even-odd rule.

Using multiple arguments is more efficient than looping through the polygons and applying `union_geometry` sequentially.

`alg` can only be `MartinezRuedaAlg()`.
"""
union_geometry(polygon1::Path2D, others::Vararg{Path2D}; options...) = 
    union_geometry(MartinezRuedaAlg(), polygon1, others...; options...)

union_geometry(polygon1::Polygon, others::Vararg{Polygon}; options...) = 
    union_geometry(MartinezRuedaAlg(), polygon1, others...; options...)

union_geometry(polygon1::AbstractVector{<:Path2D}, others::Vararg{Path2D}; options...) = 
    union_geometry(MartinezRuedaAlg(), polygon1, others...; options...)

union_geometry(polygon1::AbstractVector{<:Polygon}, others::Vararg{Polygon}; options...) = 
    union_geometry(MartinezRuedaAlg(), polygon1, others...; options...)

function union_geometry(
    alg::MartinezRuedaAlg, polygon1::Path2D, others::Vararg{Path2D}
    ; options...
    )
    martinez_rueda_algorithm(UNION_CRITERIA, polygon1, others...; options...)
end

function union_geometry(
    alg::MartinezRuedaAlg, polygon1::Polygon, others::Vararg{Polygon}
    ; options...
    )
    martinez_rueda_algorithm(UNION_CRITERIA, polygon1, others...; options...)
end

function union_geometry(
    alg::MartinezRuedaAlg, subjects::AbstractVector{<:Path2D}, others::Vararg{Path2D}
    ; options...
    )
    martinez_rueda_algorithm(UNION_CRITERIA, subjects, others...; options...)
end

function union_geometry(
    alg::MartinezRuedaAlg, subjects::AbstractVector{<:Polygon}, others::Vararg{Polygon}
    ; options...
    )
    martinez_rueda_algorithm(UNION_CRITERIA, subjects, others...; options...)
end

"""
    xor_geometry([alg=MartinezRuedaAlg()], subject, others...; atol=default_atol, rtol=default_rtol)

General case of polygon xor: in one polygon or the other but not a;;.
Possibly returns holes but does not classify them as holes.

The polygons can be of type:
- `Vector{<:Tuple{T, T}} where T<:AbstractFloat`. Output will be `Vector{Vector{Tuple{T, T}}}`. Possibly returns holes but does not classify them as holes.
- `Polygon{T<:AbstractFloat}`. Output will be `Vector{Polygon{T}}`. This will classify holes and match them to the smallest possible parent.
- The subject can also be the same type as the output. This enables the output to be an input. 
    Note that in this case the `subject` list is treated as a single polygon.
    If any of the polygons overlap, this is equivalent to passing a self-intersecting polygon and
    some areas might be classified as holes according to the even-odd rule.
Segment2D
Using multiple arguments is more efficient than looping through the polygons and applying `xor_geometry` sequentially.

`alg` can only be `MartinezRuedaAlg()`.
"""
xor_geometry(polygon1::Path2D, others::Vararg{Path2D}; options...) = 
    xor_geometry(MartinezRuedaAlg(), polygon1, others...; options...)

xor_geometry(polygon1::Polygon, others::Vararg{Polygon}; options...) = 
    xor_geometry(MartinezRuedaAlg(), polygon1, others...; options...)

xor_geometry(polygon1::AbstractVector{<:Path2D}, others::Vararg{Path2D}; options...) = 
    xor_geometry(MartinezRuedaAlg(), polygon1, others...; options...)

xor_geometry(polygon1::AbstractVector{<:Polygon}, others::Vararg{Polygon}; options...) = 
    xor_geometry(MartinezRuedaAlg(), polygon1, others...; options...)

function xor_geometry(
    alg::MartinezRuedaAlg, polygon1::Path2D, others::Vararg{Path2D}
    ; options...
    )
    martinez_rueda_algorithm(XOR_CRITERIA, polygon1, others...; options...)
end

function xor_geometry(
    alg::MartinezRuedaAlg, polygon1::Polygon, others::Vararg{Polygon}
    ; options...
    )
    martinez_rueda_algorithm(XOR_CRITERIA, polygon1, others...; options...)
end

function xor_geometry(
    alg::MartinezRuedaAlg, subjects::AbstractVector{<:Path2D}, others::Vararg{Path2D}
    ; options...
    )
    martinez_rueda_algorithm(XOR_CRITERIA, subjects, others...; options...)
end

function xor_geometry(
    alg::MartinezRuedaAlg, subjects::AbstractVector{<:Polygon}, others::Vararg{Polygon}
    ; options...
    )
    martinez_rueda_algorithm(XOR_CRITERIA, subjects, others...; options...)
end

#### Convex
"""
    intersect_convex([alg=ChasingEdgesAlg()], polygon1, polygon2; atol=default_atol)

Find the intersection points of convex polygons `polygon1` and `polygon2`.
They are not guaranteed to be unique.

The polygons can only be `Vector{<:Tuple{T, T}} where T<:AbstractFloat`. Holes are not supported.

A major assumption is that convex polygons only have one area of intersection.
This fails in the general case with non-convex polygons. 
See `intersect_geometry` for a more general `O(nm)` algorithm.

`alg` can either be `ChasingEdgesAlg()`, `MartinezRuedaAlg()`, `PointSearchAlg()` or `WeilerAthertonAlg()`.
The two general algorithms `MartinezRuedaAlg` and `WeilerAthertonAlg` will throw an error if more
than one area of intersection is found.
"""
intersect_convex(polygon1::Path2D, polygon2::Path2D; options...) = 
    intersect_convex(ChasingEdgesAlg(), polygon1, polygon2; options...)

function intersect_convex(alg::ChasingEdgesAlg, polygon1::Path2D, polygon2::Path2D,; options...)
    chasing_edges_algorithm(polygon1, polygon2; options...)
end

function intersect_convex(
    alg::WeilerAthertonAlg, polygon1::Path2D{T}, polygon2::Path2D{T}
    ; options...
    ) where T
    regions = weiler_atherton_algorithm(polygon1, polygon2; options...)
    if isempty(regions)
        return Point2D{T}[]
    end
    @assert length(regions) == 1 "Convex polygons can only have one region of intersection; $(length(regions)) found"
    regions[1]
end

function intersect_convex(
    alg::MartinezRuedaAlg, polygon1::Path2D{T}, polygon2::Path2D{T}
    ; options...
    ) where T
    regions = martinez_rueda_algorithm(INTERSECTION_CRITERIA, polygon1, polygon2; options...)
    if isempty(regions)
        return Point2D{T}[]
    end
    @assert length(regions) == 1 "Convex polygons can only have one region of intersection; $(length(regions)) found"
    regions[1]
end

function intersect_convex(
    alg::PointSearchAlg, polygon1::Path2D, polygon2::Path2D
    ; atol::AbstractFloat=default_atol
    )
    #https://www.swtestacademy.com/intersection-convex-polygons-algorithm/
    intersection_points = intersect_edges(polygon1, polygon2; atol=atol)
    
    i1_in_2 = [contains(polygon2, p) for p in polygon1]
    p1_in_2 = polygon1[i1_in_2]
    i2_in_1 = [contains(polygon1, p) for p in polygon2]
    p2_in_1 = polygon2[i2_in_1]

    if isempty(intersection_points) && isempty(p1_in_2) && isempty(p2_in_1)
        return intersection_points
    end

    points = vcat(intersection_points, p1_in_2, p2_in_1)
    sort_counter_clockwise!(points)
end
