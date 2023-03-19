using Test
using PolygonAlgorithms
using PolygonAlgorithms: translate

@testset verbose = true "PolygonAlgorithms" begin
    include("linked_list.jl")
    include("area-centroid.jl")
    include("intersect.jl")
    include("point_in_polygon.jl")
    include("intersect_convex.jl")
    include("intersect_concave.jl")
    include("convex_hull.jl")
end