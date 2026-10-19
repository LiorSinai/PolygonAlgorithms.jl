using PolygonAlgorithms: get_orientation, COLINEAR, CLOCKWISE, COUNTER_CLOCKWISE
using PolygonAlgorithms: on_segment, is_same_point
using PolygonAlgorithms: isless_orientation, isless_polar_angle
using PolygonAlgorithms: sort_counter_clockwise, sort_counter_clockwise!, sort_clockwise, sort_clockwise!
using PolygonAlgorithms: cross_product, dot_product, inner_angle, in_half_plane
using PolygonAlgorithms: is_vertical, is_above_or_on
using Test

@testset "orientation" begin

@testset "get_orientation" begin
    # left turn is counter-clockwise, right turn is clockwise
    @test get_orientation((0.0, 0.0), (1.0, 0.0), (1.0, 1.0)) == COUNTER_CLOCKWISE
    @test get_orientation((0.0, 0.0), (1.0, 0.0), (1.0, -1.0)) == CLOCKWISE
    @test get_orientation((0.0, 0.0), (1.0, 0.0), (2.0, 0.0)) == COLINEAR
    # r doubles back along the line
    @test get_orientation((0.0, 0.0), (1.0, 0.0), (-1.0, 0.0)) == COLINEAR
    # p and q are the same point
    @test get_orientation((1.0, 1.0), (1.0, 1.0), (5.0, 3.0)) == COLINEAR
end

@testset "get_orientation distance" begin
    # colinear if the distance of r from the line through p and q is at most atol=1e-6
    @test get_orientation((0.0, 0.0), (1.0, 1.0), (2.0, 2.0)) == COLINEAR
    @test get_orientation((0.0, 0.0), (1.0, 1.0), (2.0, 2.0 + 1e-6)) == COLINEAR # distance is 7.1e-7
    @test get_orientation((0.0, 0.0), (1.0, 1.0), (2.0, 2.0 + 2e-6)) == COUNTER_CLOCKWISE # distance is 1.4e-6
    @test get_orientation((0.0, 0.0), (500.0, 500.0), (1000.0, 1000.0)) == COLINEAR
    # the distance does not scale with the size of the coordinates
    @test get_orientation((0.0, 0.0), (500.0, 500.0), (1000.0, 1000.001)) == COUNTER_CLOCKWISE # distance is 7.1e-4
    @test get_orientation((0.0, 0.0), (500.0, 500.0), (1000.0, 1000.001); atol=1e-3) == COLINEAR
    @test get_orientation((0.0, 0.0), (1e9, 1e9), (2e9, 2e9 + 10)) == COUNTER_CLOCKWISE # distance is 7.1
    @test get_orientation((0.0, 0.0), (1e9, 1e9), (2e9, 2e9 + 10); atol=10.0) == COLINEAR
    # small coordinates: the cross product is tiny (1e-12) but the distance is 1e-6 * 100
    @test get_orientation((0.0, 0.0), (1e-6, 0.0), (2e-6, 1e-4)) == COUNTER_CLOCKWISE
end

@testset "get_orientation asymmetric" begin
    # p → q is the base line. Cyclic permutations have the same cross product but a different base.
    p, q, r = (0.0, 0.0), (10.0, 0.0), (10.5, 5e-7)
    @test get_orientation(p, q, r) == COLINEAR # r is 5e-7 from the line pq of length 10
    # same triangle, but p is 1e-5 from the line qr of length 0.5
    @test get_orientation(q, r, p) == COUNTER_CLOCKWISE
end

@testset "on_segment" begin
    segment = ((0.0, 0.0), (10.0, 0.0))
    @test on_segment((5.0, 0.0), segment)
    @test on_segment((0.0, 0.0), segment) # start point
    @test on_segment((10.0, 0.0), segment) # end point
    @test on_segment((5.0, 5e-7), segment) # within atol of the line
    @test !on_segment((5.0, 2e-6), segment) # beyond atol of the line
    @test on_segment((10.0 + 5e-7, 0.0), segment) # within atol past the end
    @test !on_segment((10.0 + 2e-6, 0.0), segment) # beyond atol past the end
    @test !on_segment((15.0, 0.0), segment) # colinear but not on the segment
    @test !on_segment((-5.0, 0.0), segment)

    segment = ((1.0, 1.0), (3.0, 5.0)) # diagonal
    @test on_segment((2.0, 3.0), segment)
    @test !on_segment((2.0, 3.1), segment)
    @test !on_segment((4.0, 7.0), segment) # colinear but not on the segment

    # with orientation already known
    segment = ((0.0, 0.0), (10.0, 0.0))
    @test on_segment((5.0, 0.0), segment, true)
    @test !on_segment((5.0, 0.0), segment, false)
    @test !on_segment((15.0, 0.0), segment, true) # outside the bounds
end

@testset "on_segment rectangular bounds" begin
    # The bounds are a rectangle expanded by atol, not a half-circle of radius atol at each end.
    # So a point just past a corner of the expanded rectangle can be on the segment while being
    # further than atol from the segment.
    segment = ((0.0, 0.0), (10.0, 0.0))
    atol = 1.0
    point = (10.9, 0.9)
    @test on_segment(point, segment; atol=atol) # more correct version would return false
    # the distance to the closest point on the segment is the distance to the end point
    @test !is_same_point(point, segment[2]; atol=atol) # distance is 1.27 > atol

    # For a diagonal segment this happens on the line itself, even for a small atol.
    segment = ((0.0, 0.0), (1.0, 1.0))
    point = (1.0 + 0.9e-6, 1.0 + 0.9e-6)
    @test on_segment(point, segment; atol=1e-6) # more correct version would return false
    @test !is_same_point(point, segment[2]; atol=1e-6) # distance is 1.27e-6 > atol
end

@testset "isless_orientation" begin
    p0 = (0.0, 0.0)
    # p is less than q if (p0, p, q) is counter-clockwise
    @test isless_orientation((1.0, 0.0), (0.0, 1.0), p0)
    @test !isless_orientation((0.0, 1.0), (1.0, 0.0), p0)
    # if colinear, the closer point is less
    @test isless_orientation((1.0, 1.0), (2.0, 2.0), p0)
    @test !isless_orientation((2.0, 2.0), (1.0, 1.0), p0)
    # sort points by angle around the bottom-left point, as for a convex hull
    points = [(1.0, 1.0), (0.0, 1.0), (1.0, 0.0), (0.5, 0.25)]
    sorted = sort(points; lt=(p, q) -> isless_orientation(p, q, p0))
    @test sorted == [(1.0, 0.0), (0.5, 0.25), (1.0, 1.0), (0.0, 1.0)]
end

@testset "isless_polar_angle" begin
    p0 = (0.0, 0.0)
    @test isless_polar_angle((1.0, 0.0), (0.0, 1.0), p0) # 0 < π/2
    @test !isless_polar_angle((0.0, 1.0), (1.0, 0.0), p0)
    @test isless_polar_angle((-1.0, -1.0), (1.0, -1.0), p0) # -3π/4 < -π/4
    # equal angles: the closer point is less
    @test isless_polar_angle((1.0, 1.0), (2.0, 2.0), p0)
    @test !isless_polar_angle((2.0, 2.0), (1.0, 1.0), p0)
end

@testset "sort clockwise and counter-clockwise" begin
    square = [(1.0, 1.0), (-1.0, -1.0), (-1.0, 1.0), (1.0, -1.0)]
    expected_ccw = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
    expected_cw = reverse(expected_ccw)

    @test sort_counter_clockwise(square) == expected_ccw
    @test sort_clockwise(square) == expected_cw
    @test square == [(1.0, 1.0), (-1.0, -1.0), (-1.0, 1.0), (1.0, -1.0)] # not modified

    points = copy(square)
    sort_counter_clockwise!(points)
    @test points == expected_ccw
    points = copy(square)
    sort_clockwise!(points)
    @test points == expected_cw

    # sorted around the mean, not the origin
    shifted = PolygonAlgorithms.translate(square, (10.0, 20.0))
    @test sort_counter_clockwise(shifted) == PolygonAlgorithms.translate(expected_ccw, (10.0, 20.0))
end

@testset "cross_product" begin
    edge1 = ((0.0, 0.0), (1.0, 0.0))
    @test cross_product(edge1, ((0.0, 0.0), (0.0, 1.0))) == 1.0 # counter-clockwise
    @test cross_product(edge1, ((0.0, 0.0), (0.0, -1.0))) == -1.0 # clockwise
    @test cross_product(edge1, ((5.0, 5.0), (7.0, 5.0))) == 0.0 # parallel
    @test cross_product(((0.0, 0.0), (2.0, 0.0)), ((0.0, 0.0), (0.0, 3.0))) == 6.0 # area of the parallelogram

    # the three point version is the same as the two edge version
    p, q, r = (1.0, 2.0), (4.0, 3.0), (2.0, 7.0) # left turn
    @test cross_product(p, q, r) > 0
    @test get_orientation(p, q, r) == COUNTER_CLOCKWISE
    @test cross_product(p, q, r) == cross_product((p, q), (q, r))
    p, q, r = (1.0, 2.0), (4.0, 3.0), (5.0, 1.0) # right turn
    @test cross_product(p, q, r) < 0
    @test get_orientation(p, q, r) == CLOCKWISE
    @test cross_product(p, q, r) == cross_product((p, q), (q, r))
end

@testset "dot_product and inner_angle" begin
    q = (0.0, 0.0)
    @test dot_product((1.0, 0.0), q, (0.0, 1.0)) == 0.0
    @test dot_product((2.0, 0.0), q, (3.0, 4.0)) == 6.0
    @test dot_product((1.0, 2.0), (1.0, 1.0), (1.0, 3.0)) == 2.0 # vectors start at q

    @test inner_angle((1.0, 0.0), q, (0.0, 1.0)) ≈ π/2
    @test inner_angle((1.0, 0.0), q, (1.0, 1.0)) ≈ π/4
    @test inner_angle((1.0, 0.0), q, (-1.0, 0.0)) ≈ π
    @test inner_angle((0.0, 1.0), q, (1.0, 0.0)) ≈ π/2 # not signed
    @test inner_angle((1.0, 1.0), (1.0, 0.0), (2.0, 0.0)) ≈ π/2
end

@testset "in_half_plane" begin
    edge = ((0.0, 0.0), (1.0, 0.0))
    # counter-clockwise polygon: inside is to the left
    @test in_half_plane(edge, (0.5, 1.0))
    @test !in_half_plane(edge, (0.5, -1.0))
    # clockwise polygon: inside is to the right
    @test !in_half_plane(edge, (0.5, 1.0), false)
    @test in_half_plane(edge, (0.5, -1.0), false)
    # the edge is extended to infinity
    @test in_half_plane(edge, (100.0, 1.0))
    # on the border
    @test in_half_plane(edge, (2.0, 0.0))
    @test !in_half_plane(edge, (2.0, 0.0); on_border_is_inside=false)
end

@testset "is_vertical" begin
    @test is_vertical(((1.0, 0.0), (1.0, 5.0)))
    @test is_vertical((1.0, 0.0), (1.0, 5.0))
    @test is_vertical(((1.0, 0.0), (1.0 + eps(1.0), 5.0)))
    @test !is_vertical(((1.0, 0.0), (1.0 + 1e-10, 5.0))) # uses eps, not atol
    @test !is_vertical(((0.0, 0.0), (1.0, 0.0)))
    @test is_vertical(((1.0f0, 0.0f0), (1.0f0, 5.0f0))) # Float32
end

@testset "is_above_or_on" begin
    segment = ((0.0, 0.0), (10.0, 10.0))
    @test is_above_or_on((5.0, 6.0), segment)
    @test !is_above_or_on((5.0, 4.0), segment)
    @test is_above_or_on((5.0, 5.0), segment) # on
    @test is_above_or_on((5.0, 5.0 - 1e-7), segment) # within atol
    @test !is_above_or_on((5.0, 5.0 - 1e-5), segment)
    @test is_above_or_on((20.0, 21.0), segment) # beyond the end: extended as a line
    # vertical segment: above the top
    segment = ((0.0, 0.0), (0.0, 10.0))
    @test is_above_or_on((0.0, 11.0), segment)
    @test is_above_or_on((3.0, 10.0), segment)
    @test !is_above_or_on((0.0, 5.0), segment)
end

end
