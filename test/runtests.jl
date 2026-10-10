using Test
using PolygonAlgorithms
using PolygonAlgorithms: translate, rotate, Point2D, PointSet

@testset verbose = true "PolygonAlgorithms" begin
    include("utils.jl")
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