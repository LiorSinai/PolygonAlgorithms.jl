using PolygonAlgorithms: PointSet, translate

@testset "edge intersections" begin

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
    # these don't intersect, but because of error tolerances they look like they do
    seg1 = ((0.2182133588430426, 0.42855583950407183), (0.6180465602541201, 0.8984103075263743))
    seg2 = ((0.830526104459048, 0.7943013174345828), (0.08723085309732159, 0.35031746510063355))
    @test do_intersect(seg1, seg2)
    @test all(intersect_geometry(seg1, seg2) .≈ (0.21821313839369444, 0.42855558044826847))

    # these give conflicting answers (they don't intersect)
    seg1 = ((0.950520984763203, 0.9955447269810986), (0.8147603032908444, 0.0733985670511772))
    seg2 = ((0.9238819633881898, 0.8146036454030574), (0.7570281114911228, 0.8580955127719763))
    @test !do_intersect(seg1, seg2; atol=1e-7) # fails with bigger atol
    @test isnothing(intersect_geometry(seg1, seg2))

    seg1 = ((0.4856517528012627, 0.8793324264376292), (0.5066211310466722, 0.03827488541514157))
    seg2 = ((0.6514230424916896, 0.8877275468529608), (0.5042397387902877, 0.13379459881948708))
    @test !do_intersect(seg1, seg2; atol=1e-7) # fails with bigger atol
    @test isnothing(intersect_geometry(seg1, seg2))
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