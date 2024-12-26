include("boolean/weiler_atherton.jl")
include("boolean/chasing_edges.jl")
include("boolean/martinez_rueda.jl")

abstract type PolygonIntersectionAlgorithm end

struct PointSearchAlg <: PolygonIntersectionAlgorithm end
struct ChasingEdgesAlg <: PolygonIntersectionAlgorithm end
struct WeilerAthertonAlg <: PolygonIntersectionAlgorithm end
struct MartinezRuedaAlg <: PolygonIntersectionAlgorithm end

# General

"""
    intersect_geometry(polygon1::Vector{<:Point2D}, polygon2::Vector{<:Point2D}, alg=MartinezRuedaAlg())

General case of polygon intersection. 

`alg` can either be `MartinezRuedaAlg()` or `WeilerAthertonAlg()`.
There are slight differences in the result depending on the algorithm.
The `WeilerAthertonAlg` tends to be more robust to numeric inaccuracies.

For `n` and `m` vertices on polygon 1 and 2 respectively:
- `MartinezRuedaAlg`:
    - Time complexity: `O((n+m+k)log(n+m))` where `k` is the total number of intersections.
    - Returns regions and edges of intersection. Does not return single points of intersection.
    - Works for convex and concave polygons including with holes and self-intersections.
    - Limitations:
        1. It can fail for improper polygons: polygons with lines sticking out.
        2. It is sensitive to numeric inaccuracies e.g. a line that is almost vertical or tiny regions of intersection.
- `WeilerAthertonAlg`:
    - Time complexity: `O(nm)`. 
    - Returns regions, edges and single points of intersection. Only returns the larger type if one is within another e.g. an edge is also part of a region.
    - Limitations
        1. This version does not cater for holes.
        2. It can fail completely for self-intersecting polygons.

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
    difference_geometry(polygon1::Vector{<:Point2D}, polygon2::Vector{<:Point2D}, alg=WeilerAthertonAlg())

General case of polygon difference: points in `polygon1` that are not in `polygon2`.

`alg` can only be `MartinezRuedaAlg()`. Runs in `O((n+m+k)log(n+m))`.
"""
function difference_geometry(
    polygon1::Polygon2D, polygon2::Polygon2D, alg::MartinezRuedaAlg=MartinezRuedaAlg()
    ; options...
    )
    martinez_rueda_algorithm(polygon1, polygon2, DIFFERENCE_CRITERIA)
end

"""
    union_geometry(polygon1::Vector{<:Point2D}, polygon2::Vector{<:Point2D}, alg=WeilerAthertonAlg())

General case of polygon union: in both polygons.

`alg` can only be `MartinezRuedaAlg()`. Runs in `O((n+m+k)log(n+m))`. Possibly returns holes but does not classify them as holes.
"""
function union_geometry(
    polygon1::Polygon2D, polygon2::Polygon2D, alg::MartinezRuedaAlg=MartinezRuedaAlg()
    ; options...
    )
    martinez_rueda_algorithm(polygon1, polygon2, UNION_CRITERIA; options...)
end

"""
    xor_geometry(polygon1::Vector{<:Point2D}, polygon2::Vector{<:Point2D}, alg=WeilerAthertonAlg())

General case of polygon xor: in one polygon or the other but not both. 

`alg` can only be `MartinezRuedaAlg()`. Runs in `O((n+m+k)log(n+m))`. Possibly returns holes but does not classify them as holes.
"""
function xor_geometry(
    polygon1::Polygon2D, polygon2::Polygon2D, alg::MartinezRuedaAlg=MartinezRuedaAlg()
    ; options...
    )
    martinez_rueda_algorithm(polygon1, polygon2, XOR_CRITERIA; options...)
end

#### Convex
"""
    intersect_convex(polygon1::Vector{<:Point2D}, polygon2::Vector{<:Point2D}, alg=ChasingEdgesAlg())

Find the intersection points of convex polygons `polygon1` and `polygon2`.
They are not guaranteed to be unique.

A major assumption is that convex polygons only have one area of intersection.
This fails in the general case with non-convex polygons. 
See `intersect_geometry` for a more general `O(nm)` algorithm.

`alg` can either be `ChasingEdgesAlg()`, `MartinezRuedaAlg()`, `PointSearchAlg()` or `WeilerAthertonAlg()`.`

For `n` and `m` vertices on polygon 1 and 2 respectively:
- `ChasingEdgesAlg` (default):
    - Time complexity: `O(n+m)`. 
    - Algorithm is from "A New Linear Algorithm for Intersecting Convex Polygons" (1981) by Joseph O'Rourke et. al.
    - Reference: https://www.cs.jhu.edu/~misha/Spring16/ORourke82.pdf
- `MartinezRuedaAlg`:
    - Time complexity: `O((n+m+k)log(n+m))` where `k` is the total number of intersections.
    - Works for convex and concave polygons including with holes and self-intersections.
    - Will throw an error if there is more than one region of intersection.
- `PointSearchAlg`:
    - Time complexity: `O(nm)`. 
    - Algorithm: (1) Intersect all edge pairs. (2) Check all points in the other polygon. (3) Sort results counter-clockwise. 
    - For general non-intersecting polygons, the intersection points are valid but the order is not.
- `WeilerAthertonAlg`:
    - Time complexity: `O(nm)`. 
    - Designed for more complex concave polygons with multiple areas of intersection.
    - Will throw an error if there is more than one region of intersection.
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
    ; atol::AbstractFloat=1e-6
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
    points = sort_counter_clockwise(points)
    points
end
