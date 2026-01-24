using Test
using PolygonAlgorithms
using PolygonAlgorithms: translate, rotate, Point2D, PointSet

function are_regions_equal(r1::Vector{<:Vector{<:Point2D}}, r2::Vector{<:Vector{<:Point2D}})
    if length(r1) != length(r2)
        return false
    end
    r1_sets = [PointSet(r) for r in r1]
    r2_sets = [PointSet(r) for r in r2]
    issetequal(r1_sets, r2_sets)
end

@testset verbose = true "PolygonAlgorithms" begin
    # data structures
    include("data_structures/linked_list.jl")
    include("data_structures/point_set.jl")
    # segments
    include("intersect.jl")
    include("line_sweep.jl")
    # polygon 
    include("convex_hull.jl")
    include("area-centroid.jl")
    include("bounds.jl")
    include("point_in_polygon.jl")
    include("data_structures/polygon.jl")
    # boolean
    include("boolean/martinez_rueda.jl")
    include("intersect_convex.jl")
    include("intersect_concave.jl")
    include("polygon_boolean.jl")
    include("polygon_boolean_holes.jl")
    include("polygon_boolean_multi.jl")
end