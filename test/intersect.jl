using PolygonAlgorithms: PointSet, translate, COLINEAR, COUNTER_CLOCKWISE

@testset "segment intersections" begin

@testset "intersect segments" begin
    seg1 = ((1.0, 1.0), (2.0, 3.0));

    seg2 = ((1.0, 3.0), (4.0, 1.0));
    @test do_intersect(seg1, seg2)
    @test intersect_geometry(seg1, seg2) == (1.75, 2.5)

    seg2 = ((1.0, 1.0), (4.0, 3.0)); # point is on vertex
    @test do_intersect(seg1, seg2)
    @test intersect_geometry(seg1, seg2) == (1.0, 1.0)

    seg2 = ((1.0, 0.0), (4.0, 3.0));
    @test !do_intersect(seg1, seg2)
    @test isnothing(intersect_geometry(seg1, seg2))

    seg2 = ((2.0, 2.0), (5.0, 0.0)); # point is on seg1 but not seg2
    @test !do_intersect(seg1, seg2)
    @test isnothing(intersect_geometry(seg1, seg2))

    seg2 = ((1.0, 0.0), (2.0, 2.0)); # parallel
    @test !do_intersect(seg1, seg2)
    @test isnothing(intersect_geometry(seg1, seg2))

    seg2 = ((1.5, 2.0), (2.5, 4.0)); # parallel & overlap
    @test do_intersect(seg1, seg2)
    @test isnothing(intersect_geometry(seg1, seg2))
end;

@testset "intersect 90 degree segments" begin
    seg1 = ((1.0, 1.0), (1.0, 4.0));

    seg2 = ((-1.0, 2.0), (4.0, 2.0)); # horizontal
    @test do_intersect(seg1, seg2)
    @test intersect_geometry(seg1, seg2) == (1.0, 2.0)

    seg2 = ((1.0, 2.0), (4.0, 2.0)); # point is on vertex
    @test do_intersect(seg1, seg2)
    @test intersect_geometry(seg1, seg2) == (1.0, 2.0)

    seg2 = ((2.0, 2.0), (4.0, 2.0)); # point is on seg1 but not seg2
    @test !do_intersect(seg1, seg2)
    @test isnothing(intersect_geometry(seg1, seg2))

    seg2 = ((3.0, 1.0), (3.0, 4.0)); # parallel
    @test !do_intersect(seg1, seg2)
    @test isnothing(intersect_geometry(seg1, seg2))
end;

@testset "intersect floating point error" begin
    seg1 = ((-1.0, 0.25), (1.0, 0.25))
    seg2 = ((0.4, 0.4), (0.4, -0.4))

    @test do_intersect(seg1, seg2)
    @test all(intersect_geometry(seg1, seg2) .≈ (0.4, 0.25))
end

@testset "intersect segments almost" begin
    # These don't intersect. The lines cross between 3e-7 and 6e-7 past the end of one segment.
    # So they intersect for atol=1e-6 but not for atol=1e-7.
    # do_intersect and intersect_geometry should agree.
    seg1 = ((0.2182133588430426, 0.42855583950407183), (0.6180465602541201, 0.8984103075263743))
    seg2 = ((0.830526104459048, 0.7943013174345828), (0.08723085309732159, 0.35031746510063355))
    @test do_intersect(seg1, seg2)
    @test intersect_geometry(seg1, seg2) == seg1[1] # 3.4e-7 before the start, so clamped to it
    @test !do_intersect(seg1, seg2; atol=1e-7)
    @test isnothing(intersect_geometry(seg1, seg2; atol=1e-7))

    seg1 = ((0.950520984763203, 0.9955447269810986), (0.8147603032908444, 0.0733985670511772))
    seg2 = ((0.9238819633881898, 0.8146036454030574), (0.7570281114911228, 0.8580955127719763))
    @test do_intersect(seg1, seg2)
    @test all(intersect_geometry(seg1, seg2) .≈ (0.9238823672633001, 0.8146035401295949)) # 4.2e-7 before the start of seg2
    @test !do_intersect(seg1, seg2; atol=1e-7)
    @test isnothing(intersect_geometry(seg1, seg2; atol=1e-7))

    seg1 = ((0.4856517528012627, 0.8793324264376292), (0.5066211310466722, 0.03827488541514157))
    seg2 = ((0.6514230424916896, 0.8877275468529608), (0.5042397387902877, 0.13379459881948708))
    @test do_intersect(seg1, seg2)
    @test all(intersect_geometry(seg1, seg2) .≈ (0.5042396322136217, 0.13379405289030877)) # 5.6e-7 past the end of seg2
    @test !do_intersect(seg1, seg2; atol=1e-7)
    @test isnothing(intersect_geometry(seg1, seg2; atol=1e-7))
end

@testset "intersect tolerances are distances" begin
    # small segments that clearly cross. The determinant is an area so must not be compared to atol.
    s = 5e-4
    seg1 = ((0.0, 0.0), (s, s))
    seg2 = ((0.0, s), (s, 0.0))
    @test do_intersect(seg1, seg2)
    @test all(intersect_geometry(seg1, seg2) .≈ (s/2, s/2))

    # long segment. The other segment stops 500*atol past its end.
    seg1 = ((0.0, 0.0), (1000.0, 0.0))
    seg2 = ((1000.0005, -1.0), (1000.0005, 1.0))
    @test !do_intersect(seg1, seg2)
    @test isnothing(intersect_geometry(seg1, seg2))

    # within atol past the end: clamped to the end point
    seg2 = ((1000.0000005, -1.0), (1000.0000005, 1.0))
    @test do_intersect(seg1, seg2)
    @test intersect_geometry(seg1, seg2) == (1000.0, 0.0)
end

@testset "orientation is rotation invariant" begin
    # r is 5e-4 = 500*atol off the line
    p, q, r = (0.0, 0.0), (1000.0, 0.0), (2000.0, 5e-4)
    @test get_orientation(p, q, r) == COUNTER_CLOCKWISE
    p_rot, q_rot, r_rot = PolygonAlgorithms.rotate([p, q, r], π/4)
    @test get_orientation(p_rot, q_rot, r_rot) == COUNTER_CLOCKWISE
    # within atol of the line
    @test get_orientation(p, q, (2000.0, 5e-7)) == COLINEAR
end

@testset "intersect tolerances disagree" begin
    # example from Clipper test 141
    segment1 = ((340.0, 710.0), (1450.0, 860.0))
    segment2 = ((830.0, 1150.0), (1540.0, 600.0))
    segment3 = ((1180.0, 170.0), (1260.0, 1040.0))
    p12 = intersect_geometry(segment1, segment2)
    p13 = intersect_geometry(segment1, segment3)
    p23 = intersect_geometry(segment2, segment3)
    s123 = (p12, p13) # norm is 5e-5
    s23 = (p23, segment3[2]) # norm is 209
    @test isnothing(intersect_geometry(s123, s23))
    @test isnothing(intersect_geometry(s123, s23; atol=1e-4))
    @test all(intersect_geometry(s123, s23; atol=4.3e-5) .≈ (1240.8493236866939, 831.7363950927964))
    @test isnothing(intersect_geometry(s123, s23; atol=1e-5))
    @test isnothing(intersect_geometry(s123, s23; atol=1e-6))
end

@testset "intersect rectangles" begin
    poly1 = [
        (0.0, 0.0), (0.0, 2.0), (2.0, 2.0), (2.0, 0.0)
    ]

    points = intersect_edges(poly1, poly1)
    @test PointSet(points) == PointSet(poly1)

    # overlap
    poly2 = translate(poly1, (1.0, 1.0))
    expected = [
        (1.0, 2.0), (2.0, 1.0)
        ]
    points = intersect_edges(poly1, poly2)
    @test PointSet(points) == PointSet(expected)

    # overlap along edges
    poly2 = translate(poly1, (1.0, 0.0))
    expected = [
        (1.0, 2.0), (2.0, 2.0), (2.0, 0.0), (1.0, 0.0)
        ]
    points = intersect_edges(poly1, poly2)
    @test PointSet(points) == PointSet(expected)

    # single point of itersection
    poly2 = translate(poly1, (2.0, 2.0))
    expected = [(2.0, 2.0)]
    points = intersect_edges(poly1, poly2)
    @test PointSet(points) == PointSet(expected)

    # intersect along edge
    poly2 = translate(poly1, (2.0, 1.0))
    expected = [
        (2.0, 2.0),
        (2.0, 1.0)
        ]
    points = intersect_edges(poly1, poly2)
    @test PointSet(points) == PointSet(expected)

    # no intersection
    poly2 = translate(poly1, (3.0, 3.0))
    expected = []
    points = intersect_edges(poly1, poly2)
    @test isempty(points)
end

@testset "cross" begin
    poly1 = [
        (0.0, 0.0), (0.0, 2.0), (1.0, 2.0), (1.0, 0.0)
    ]
    poly2 = [
        (-1.0, 1.5), (2.0, 1.5), (2.0, 0.5), (-1.0, 0.5)
    ]
    expected = [
        (0.0, 1.5), (0.0, 0.5), (1.0, 1.5), (1.0, 0.5)
    ]
    points = intersect_edges(poly1, poly2)
    @test PointSet(points) == PointSet(expected)
end

@testset "star of david" begin
    h = sqrt(3)
    poly1 = translate([
        (-1.0, 0.0), (0.0, h), (1.0, 0.0)
        ], (0.0, -h/3)
    )
    poly2 = translate([
        (-1.0, 0.0), (0.0, -h), (1.0, 0.0)
        ], (0.0, h/3)
    )
    
    expected = [
        (-2/3, 0.0), 
        (+2/3, 0.0),
        (-1/3, -h/3),
        (-1/3, +h/3),
        (+1/3, -h/3),
        (+1/3, +h/3),
        ]
    points = intersect_edges(poly1, poly2)
    @test PointSet(points) == PointSet(expected)
end

@testset "generic convex" begin
    poly1 = [
        (0.0, 0.4),
        (0.3, 0.2),
        (0.7, 0.3),
        (0.9, 0.7),
        (0.4, 0.7),
    ]
    poly2 = [
        (0.3, 0.4),
        (0.9, 0.1),
        (0.7, 0.4),
        (0.2, 0.9),
        (0.1, 0.6),
    ]

    expected = 
    [
        (0.4, 0.7),
        (0.566667, 0.266667),
        (0.171429, 0.528571),
        (0.728571, 0.357143),
    ]
    points = intersect_edges(poly1, poly2)
    @test PointSet(points) == PointSet(expected)
end

end