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
    # Data structures
    include("data_structures/linked_list.jl")
    include("data_structures/point_set.jl")
    # Foundation algorithms
    include("bounds.jl")
    include("moments.jl")
    include("intersect.jl")
    include("point_in_polygon.jl")
    # Complex algorithms
    include("convex_hull.jl")
    include("line_sweep.jl")
    include("data_structures/polygon.jl")
    include("graphs/face_computation.jl")
    include("mappings.jl")
    include("boolean.jl")
end