using Test
using PolygonAlgorithms
using PolygonAlgorithms: translate, rotate, PointSet

@testset verbose = true "PolygonAlgorithms" begin
    include("linked_list.jl")
    include("point_set.jl")
    include("area-centroid.jl")
    include("bounds.jl")
    include("intersect.jl")
    include("point_in_polygon.jl")
    include("intersect_convex.jl")
    include("intersect_concave.jl")
    include("intersect_numeric.jl")
    include("intersect_hilbert.jl")
    include("convex_hull.jl")
end