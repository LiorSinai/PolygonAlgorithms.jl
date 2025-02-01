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
    intersect_geometry(polygon1::Vector{<:Point2D}, polygon2::Vector{<:Point2D}, alg=MartinezRuedaAlg())

General case of polygon intersection. 

`alg` can either be `MartinezRuedaAlg()` or `WeilerAthertonAlg()`.
There are slight differences in the result depending on the algorithm.
The `WeilerAthertonAlg` tends to be more robust to numeric inaccuracies.

Use `intersect_convex` for convex polygons for an `O(n+m)` algorithm.
"""
intersect_geometry(polygon1::Polygon2D, polygon2::Polygon2D; options...) = 
    intersect_geometry(polygon1, polygon2, MartinezRuedaAlg(); options...)

function intersect_geometry(
    polygon1::Polygon2D, polygon2::Polygon2D, alg::WeilerAthertonAlg
    ; options...
    )
    weiler_atherton_algorithm(polygon1, polygon2; options...)
end

function intersect_geometry(
    polygon1::Polygon2D, polygon2::Polygon2D, alg::MartinezRuedaAlg
    ; options...
    )
    martinez_rueda_algorithm(polygon1, polygon2, INTERSECTION_CRITERIA; options...)
end

"""
    difference_geometry(polygon1::Vector{<:Point2D}, polygon2::Vector{<:Point2D}, alg=MartinezRuedaAlg())

General case of polygon difference: points in `polygon1` that are not in `polygon2`.

`alg` can only be `MartinezRuedaAlg()`.
"""
function difference_geometry(
    polygon1::Polygon2D, polygon2::Polygon2D, alg::MartinezRuedaAlg=MartinezRuedaAlg()
    ; options...
    )
    martinez_rueda_algorithm(polygon1, polygon2, DIFFERENCE_CRITERIA)
end

"""
    union_geometry(polygon1::Vector{<:Point2D}, polygon2::Vector{<:Point2D}, alg=MartinezRuedaAlg())

General case of polygon union: in both polygons. 
Possibly returns holes but does not classify them as holes.

`alg` can only be `MartinezRuedaAlg()`.
"""
function union_geometry(
    polygon1::Polygon2D, polygon2::Polygon2D, alg::MartinezRuedaAlg=MartinezRuedaAlg()
    ; options...
    )
    martinez_rueda_algorithm(polygon1, polygon2, UNION_CRITERIA; options...)
end

"""
    xor_geometry(polygon1::Vector{<:Point2D}, polygon2::Vector{<:Point2D}, alg=MartinezRuedaAlg())

General case of polygon xor: in one polygon or the other but not both.
Possibly returns holes but does not classify them as holes.

`alg` can only be `MartinezRuedaAlg()`.
"""
function xor_geometry(
    polygon1::Polygon2D, polygon2::Polygon2D, alg::MartinezRuedaAlg=MartinezRuedaAlg()
    ; options...
    )
    martinez_rueda_algorithm(polygon1, polygon2, XOR_CRITERIA; options...)
end

#### Convex
"""
    intersect_convex(polygon1, polygon2, alg=ChasingEdgesAlg(); atol=PolygonAlgorithms.default_atol)

Find the intersection points of convex polygons `polygon1` and `polygon2`.
They are not guaranteed to be unique.

A major assumption is that convex polygons only have one area of intersection.
This fails in the general case with non-convex polygons. 
See `intersect_geometry` for a more general `O(nm)` algorithm.

`alg` can either be `ChasingEdgesAlg()`, `MartinezRuedaAlg()`, `PointSearchAlg()` or `WeilerAthertonAlg()`.
The two general algorithms `MartinezRuedaAlg` and `WeilerAthertonAlg` will throw an error if more
than one area of intersection is found.
"""
intersect_convex(polygon1::Polygon2D, polygon2::Polygon2D; options...) = 
    intersect_convex(polygon1, polygon2, ChasingEdgesAlg(); options...)

function intersect_convex(polygon1::Polygon2D, polygon2::Polygon2D, ::ChasingEdgesAlg; options...)
    chasing_edges_algorithm(polygon1, polygon2; options...)
end

function intersect_convex(
    polygon1::Polygon2D{T}, polygon2::Polygon2D{T}, ::WeilerAthertonAlg
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
    polygon1::Polygon2D{T}, polygon2::Polygon2D{T}, ::MartinezRuedaAlg
    ; options...
    ) where T
    regions = martinez_rueda_algorithm(polygon1, polygon2, INTERSECTION_CRITERIA; options...)
    if isempty(regions)
        return Point2D{T}[]
    end
    @assert length(regions) == 1 "Convex polygons can only have one region of intersection; $(length(regions)) found"
    regions[1]
end

function intersect_convex(
    polygon1::Polygon2D, polygon2::Polygon2D, ::PointSearchAlg
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
