using PolygonAlgorithms: translate
using PolygonAlgorithms: PointSearchAlg, ChasingEdgesAlg, WeilerAthertonAlg, MartinezRuedaAlg
using PolygonAlgorithms: are_equivalent_polygons

@testset "intersections convex - $alg" for alg in [
    PointSearchAlg(),
    ChasingEdgesAlg(),
    WeilerAthertonAlg(),
    MartinezRuedaAlg(), # Some results are different
]

@testset "one inside the other" begin
    poly1 = [
        (0.0, 0.0), (0.0, 3.0), (3.0, 3.0), (3.0, 0.0)
    ]
    poly2 = [
        (1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0)
    ]
    points = intersect_convex(alg, poly1, poly2)
    @test are_equivalent_polygons(points, poly2; match_reverse=true)
    points = intersect_convex(alg, poly2, poly1)
    @test are_equivalent_polygons(points, poly2; match_reverse=true)
end

@testset "rectangles" begin
    poly1 = [
        (0.0, 0.0), (0.0, 2.0), (2.0, 2.0), (2.0, 0.0)
    ]
    points = intersect_convex(alg, poly1, poly1)
    @test are_equivalent_polygons(points, poly1; match_reverse=true)

    # overlap
    poly2 = translate(poly1, (1.0, 1.0))
    expected = [
        (1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0), 
    ]
    points = intersect_convex(alg, poly1, poly2)
    @test are_equivalent_polygons(points, expected; match_reverse=true)

    # no intersection
    poly2 = translate(poly1, (4.0, 3.0))
    expected = []
    points = intersect_convex(alg, poly1, poly2)
    @test isempty(points)
end

@testset "rectangles - edge overlap" begin
    poly1 = [
        (0.0, 0.0), (0.0, 2.0), (2.0, 2.0), (2.0, 0.0)
    ]

    poly2 = translate(poly1, (2.0, 1.0))
    expected = (typeof(alg) == PolygonAlgorithms.MartinezRuedaAlg) ?
        Tuple{Float64, Float64}[] : [(2.0, 1.0), (2.0, 2.0)]
    points = intersect_convex(alg, poly1, poly2)
    @test are_equivalent_polygons(points, expected; match_reverse=true)
    points = intersect_convex(alg, poly2, poly1)
    @test are_equivalent_polygons(points, expected; match_reverse=true)

    poly2 = [
        (2.0, -1.0), (2.0, 3.0), (4.0, 3.0), (4.0, -1.0)
    ]
    expected = (typeof(alg) == PolygonAlgorithms.MartinezRuedaAlg) ?
        Tuple{Float64, Float64}[] : [(2.0, 0.0), (2.0, 2.0)]
    points = intersect_convex(alg, poly1, poly2)
    @test are_equivalent_polygons(points, expected; match_reverse=true)
    points = intersect_convex(alg, poly2, poly1)
    @test are_equivalent_polygons(points, expected; match_reverse=true)

    # edge + vertex overlap
    poly2 = translate(poly1, (1.0, 0.0))
    expected = [
        (1.0, 2.0), (2.0, 2.0), (2.0, 0.0), (1.0, 0.0)
        ]
    points = intersect_convex(alg, poly1, poly2)
    @test are_equivalent_polygons(points, expected; match_reverse=true)
    points = intersect_convex(alg, poly2, poly1)
    @test are_equivalent_polygons(points, expected; match_reverse=true)
end

@testset "rectangles - vertex overlap" begin
    poly1 = [
        (0.0, 0.0), (0.0, 2.0), (2.0, 2.0), (2.0, 0.0)
    ]
    # single point of itersection
    poly2 = translate(poly1, (2.0, 2.0))
    expected = (typeof(alg) == PolygonAlgorithms.MartinezRuedaAlg) ?
        Tuple{Float64, Float64}[] : [(2.0, 2.0)]
    points = intersect_convex(alg, poly1, poly2)
    @test are_equivalent_polygons(points, expected; match_reverse=true)

    points = intersect_convex(alg, poly2, poly1)
    @test are_equivalent_polygons(points, expected; match_reverse=true)
end

@testset "edge intersect vertex inner" begin 
    poly1 = [
        (0.0, 0.0), (0.0, 3.0), (2.0, 5.0), (5.0, 0.0)
    ]
    poly2 = [
        (0.0, 2.0), (2.0, 3.0), (6.0, 1.0)
    ];
    expected = [
        (4.222222, 1.296296),
        (3.714286, 2.142857),
        (2.0, 3.0),
        (0.0, 2.0),
    ]
    points = intersect_convex(alg, poly1, poly2)
    @test are_equivalent_polygons(points, expected; match_reverse=true)

    points = intersect_convex(alg, poly2, poly1)
    @test are_equivalent_polygons(points, expected; match_reverse=true)
end

@testset "edge-vertex pass through" begin 
    poly1 = [
        (0.0, 0.4),
        (0.3, 0.2),
        (0.7, 0.3),
        (0.9, 0.7),
        (0.4, 0.7),
    ]
    poly2 = [
        (0.3, 0.4),
        (0.7, 0.4),
        (0.2, 0.9),
        (0.1, 0.6),
    ];
    expected = [
        (0.3, 0.4),
        (0.4, 0.7),
        (0.171429, 0.528571),
        (0.7, 0.4),
    ]
    points = intersect_convex(alg, poly1, poly2)
    answer = Set([round.(p, digits=6) for p in points])
    @test issetequal(answer, expected)

    points = intersect_convex(alg, poly2, poly1)
    answer = Set([round.(p, digits=6) for p in points])
    @test issetequal(answer, expected)
end

@testset "vertex intersect vertex" begin 
    poly1 = [
        (0.0, 0.0), (0.0, 3.0), (2.0, 5.0), (5.0, 0.0)
    ]
    poly2 = [
        (0.0, 0.0), (2.0, 3.0), (6.0, 1.0)
    ];
    expected = [
        (0.0, 0.0),
        (4.545455, 0.757576),
        (3.714286, 2.142857),
        (2.0, 3.0),
    ]
    points = intersect_convex(alg, poly1, poly2)
    @test are_equivalent_polygons(points, expected; match_reverse=true)

    points = intersect_convex(alg, poly2, poly1)
    @test are_equivalent_polygons(points, expected; match_reverse=true)
end

@testset "vertex intersections" begin 
    poly1 = [
        (0.0, 0.0), (0.5, 1.0), (1.0, 0.0)
    ]
    poly2 = [
        (0.0, 1.0), (0.0, 3.0), (1.0, 3.0), (1.0, 1.0)
    ]

    # single point
    expected = (typeof(alg) == PolygonAlgorithms.MartinezRuedaAlg) ?
        Tuple{Float64, Float64}[] : [(0.5, 1.0)]
    points = intersect_convex(alg, poly1, poly2)
    @test are_equivalent_polygons(points, expected; match_reverse=true)
    points = intersect_convex(alg, poly2, poly1)
    @test are_equivalent_polygons(points, expected; match_reverse=true)

    # 2 points inside
    poly1_ = translate(poly1, ((0.0), (1.5)))
    expected = poly1_
    points = intersect_convex(alg, poly1_, poly2)
    @test are_equivalent_polygons(points, expected; match_reverse=true)
    points = intersect_convex(alg, poly2, poly1_)
    @test are_equivalent_polygons(points, expected; match_reverse=true)

    # 3 points inside
    poly1_ = translate(poly1, ((0.0), (2.0)))
    expected = poly1_
    points = intersect_convex(alg, poly1_, poly2)
    @test are_equivalent_polygons(points, expected; match_reverse=true)
    points = intersect_convex(alg, poly2, poly1_)
    @test are_equivalent_polygons(points, expected; match_reverse=true)
end

@testset "quads single vertex intersect" begin 
    poly1 = [
        (1.0, 1.0), (2.0, 4.0), (5.0, 5.0), (4.0, 2.0)
    ]
    poly2 = [
        (4.0, 2.0), (3.0, 4.0), (6.0, 4.0), (7.0, 1.0)
    ]

    expected = [(4.0, 2.0), (3.0, 4.0), (4.666666666666667, 4.0)]
    points = intersect_convex(alg, poly1, poly2)
    @test are_equivalent_polygons(points, expected; match_reverse=true)
    points = intersect_convex(alg, poly2, poly1)
    @test are_equivalent_polygons(points, expected; match_reverse=true)
end

@testset "share lines" begin 
    poly1 = [
        (5.0, 8.0), (9.0, 4.0), (1.0, 4.0),
    ]
    poly2 = [
        (5.0, 1.0), (1.0, 4.0), (9.0, 4.0)
    ]
    expected = (typeof(alg) == PolygonAlgorithms.MartinezRuedaAlg) ?
        Tuple{Float64, Float64}[] : [(1.0, 4.0), (9.0, 4.0)]
    points = intersect_convex(alg, poly1, poly2)
    @test are_equivalent_polygons(points, expected; match_reverse=true)
    points = intersect_convex(alg, poly2, poly1)
    @test are_equivalent_polygons(points, expected; match_reverse=true)

    # extra shared point on line
    poly1 = [
        (5.0, 8.0), (9.0, 4.0), (5.0, 4.0), (1.0, 4.0),
    ]
    poly2 = [
        (5.0, 1.0), (1.0, 4.0), (5.0, 4.0), (9.0, 4.0)
    ]
    expected = (typeof(alg) == PolygonAlgorithms.MartinezRuedaAlg) ?
        Tuple{Float64, Float64}[] : 
        (typeof(alg) == PolygonAlgorithms.PointSearchAlg) ? [(1.0, 4.0), (5.0, 4.0), (9.0, 4.0)] : 
        [(1.0, 4.0), (5.0, 4.0), (9.0, 4.0), (5.0, 4.0)]
    points = intersect_convex(alg, poly1, poly2)
    @test are_equivalent_polygons(points, expected; match_reverse=true)
    points = intersect_convex(alg, poly2, poly1)
    @test are_equivalent_polygons(points, expected; match_reverse=true)
end

@testset "cross" begin
    poly1 = [
        (0.0, 0.0), (0.0, 2.0), (1.0, 2.0), (1.0, 0.0)
    ]
    poly2 = [
        (-1.0, 1.5), (2.0, 1.5), (2.0, 0.5), (-1.0, 0.5)
    ]
    expected = [
        (1.0, 0.5), (1.0, 1.5), (0.0, 1.5), (0.0, 0.5)
    ]
    points = intersect_convex(alg, poly1, poly2)
    @test are_equivalent_polygons(points, expected; match_reverse=true)

    points = intersect_convex(alg, poly2, poly1)
    @test are_equivalent_polygons(points, expected; match_reverse=true)
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
    [(-0.33333333333333337, -0.5773502691896257), 
    (0.33333333333333326, -0.5773502691896257), 
    (0.6666666666666667, 0.0), 
    (0.33333333333333337, 0.5773502691896257),
     (-0.3333333333333333, 0.5773502691896257),
      (-0.6666666666666666, 0.0)]
    expected = [
        (-1/3, -h/3), (+1/3, -h/3),
        (2/3, 0.0), (+1/3, h/3),
        (-1/3, h/3), (-2/3, 0.0),  
    ]
    points = intersect_convex(alg, poly1, poly2)
    @test are_equivalent_polygons(points, expected; match_reverse=true)

    points = intersect_convex(alg, poly2, poly1)
    @test are_equivalent_polygons(points, expected; match_reverse=true)
end

@testset "multiple pockets + edge-vertex" begin 
    poly1 = [
        (1.0, 3.0),
        (1.0, 5.5),
        (2.0, 6.0),
        (5.0, 8.0),
        (6.0, 7.0),
        (8.0, 5.5),
        (8.0, 3.0),
        (6.0, 1.0),
    ]
    poly2 = [
        (0.0, 6.0),
        (5.2, 7.5),
        (7.0, 7.5),
        (9.0, 6.0),
        (8.0, 4.0),
        (4.0, 1.0),
    ];
    expected = [
        (3.058824, 2.176471), (4.695652, 1.521739), (8.0, 4.0), (8.0, 4.0),
        (8.0, 5.5), (6.0, 7.0), (5.5, 7.5), (5.2, 7.5), (3.525424, 7.016949),
        (2.0, 6.0), (1.0, 5.5), (1.0, 4.75)
    ]
    points = intersect_convex(alg, poly1, poly2)
    answer = Set([round.(p, digits=6) for p in points])
    @test are_equivalent_polygons(points, expected; match_reverse=true)

    points = intersect_convex(alg, poly2, poly1)
    answer = Set([round.(p, digits=6) for p in points])
    @test are_equivalent_polygons(points, expected; match_reverse=true)
end

end