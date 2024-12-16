using Test
using PolygonAlgorithms
using PolygonAlgorithms: translate, rotate, PointSet

function are_regions_equal(r1::Vector{<:Vector{<:Point2D}}, r2::Vector{<:Vector{<:Point2D}})
    if length(r1) != length(r2)
        return false
    end
    r1_sets = [PointSet(r) for r in r1]
    r2_sets = [PointSet(r) for r in r2]
    issetequal(r1_sets, r2_sets)
end

@testset verbose = true "PolygonAlgorithms" begin
    include("linked_list.jl")
    include("point_set.jl")
    include("area-centroid.jl")
    include("bounds.jl")
    include("intersect.jl")
    include("point_in_polygon.jl")
    include("martinez_rueda.jl")
    include("intersect_convex.jl")
    include("intersect_concave.jl")
    include("convex_hull.jl")
end