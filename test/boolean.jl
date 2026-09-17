using Test
using PolygonAlgorithms

@testset verbose = true "boolean" begin
    include("boolean/martinez_rueda.jl")
    include("boolean/intersect_convex.jl")
    include("boolean/intersect_concave.jl")
    include("boolean/pairs.jl")
    include("boolean/multi.jl")
    include("boolean/polygon_holes.jl")
end