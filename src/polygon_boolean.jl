include("boolean/weiler_atherton.jl")
include("boolean/chasing_edges.jl")
#include("martinez_rueda.jl")

abstract type PolygonIntersectionAlgorithm end

struct PointSearchAlg <: PolygonIntersectionAlgorithm end
struct ChasingEdgesAlg <: PolygonIntersectionAlgorithm end
struct WeilerAthertonAlg <: PolygonIntersectionAlgorithm end
#struct MartinezRueda <: PolygonIntersectionAlgorithm end

# General

intersect_geometry(polygon1::Polygon2D, polygon2::Polygon2D) = intersect_geometry(polygon1, polygon2, WeilerAthertonAlg())

"""
    intersect_geometry(polygon1::Vector{<:Point2D}, polygon2::Vector{<:Point2D}, alg)

Returns multiple regions, edges and single points of intersection. 
Only returns the larger type if one is within another e.g. an edge is also part of a region.

This uses the Weiler-Atherton algorithm.
It runs in `O(nm)` time where `n` and `m` are the number of vertices of polygon1 and polygon2 respectively.
Use `intersect_convex` for convex polygons for an `O(n+m)` algorithm.

Limitations
1. This version does not cater for holes.
2. It can fail completely for self-intersecting polygons.

For a more general algorithm see the Martinez-Rueda polygon clipping algorithm.
"""
function intersect_geometry(polygon1::Polygon2D, polygon2::Polygon2D, alg::WeilerAthertonAlg)
    weiler_atherton_algorithm(polygon1, polygon2)
end

# function intersect_geometry(polygon1::Polygon2D, polygon2::Polygon2D, alg::MartinezRueda)
#     martinez_rueda_algorithm(polygon1, polygon2)
# end

#### Convex
"""
    intersect_convex(polygon1::Vector{<:Point2D}, polygon2::Vector{<:Point2D}, alg=ChasingEdgesAlg())

Find the intersection points of convex polygons `polygon1` and `polygon2`.
They are not guaranteed to be unique.

A major assumption is that convex polygons only have one area of intersection.
This fails in the general case with non-convex polygons. 
See `intersect_geometry` for a more general `O(nm)` algorithm.

`alg` can either be `ChasingEdgesAlg()`, `PointSearchAlg()` or `WeilerAthertonAlg()`.

For `n` and `m` vertices on polygon 1 and 2 respectively:
- `ChasingEdgesAlg`:
    - Time complexity: `O(n+m)`. 
    - Algorithm is from "A New Linear Algorithm for Intersecting Convex Polygons" (1981) by Joseph O'Rourke et. al.
    - Reference: https://www.cs.jhu.edu/~misha/Spring16/ORourke82.pdf
- `PointSearchAlg`:
    - Time complexity: `O(nm)`. 
    - Algorithm: (1) Intersect all edge pairs. (2) Check all points in the other polygon. (3) Sort results counter-clockwise. 
    - For general non-intersecting polygons, the intersection points are valid but the order is not.
- `WeilerAthertonAlg`:
    - Time complexity: `O(nm)`. 
    - Designed for more complex concave polygons with multiple areas of intersection.
    - Will throw an error if there is more than one region of intersection.
"""
intersect_convex(polygon1::Polygon2D, polygon2::Polygon2D) = intersect_convex(polygon1, polygon2, ChasingEdgesAlg())

function intersect_convex(polygon1::Polygon2D, polygon2::Polygon2D, ::ChasingEdgesAlg)
    chasing_edges_algorithm(polygon1, polygon2)
end

function intersect_convex(polygon1::Polygon2D{T}, polygon2::Polygon2D{T}, ::WeilerAthertonAlg) where T
    regions = intersect_geometry(polygon1, polygon2)
    if isempty(regions)
        return Point2D{T}[]
    end
    @assert length(regions) == 1 "Convex polygons can only have one region of intersection; $(length(regions)) found"
    regions[1]
end

function intersect_convex(polygon1::Polygon2D, polygon2::Polygon2D, ::PointSearchAlg)
    #https://www.swtestacademy.com/intersection-convex-polygons-algorithm/
    intersection_points = intersect_edges(polygon1, polygon2)
    
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
