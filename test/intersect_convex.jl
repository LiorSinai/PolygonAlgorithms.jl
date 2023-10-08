using PolygonAlgorithms: translate

@testset "convex intersections -$alg" for alg in [
    PolygonAlgorithms.PointSearchAlg(),
    PolygonAlgorithms.ChasingEdgesAlg(),
    PolygonAlgorithms.WeilerAthertonAlg(),
]

@testset "one inside the other" begin
    poly1 = [
        (0.0, 0.0), (0.0, 3.0), (3.0, 3.0), (3.0, 0.0)
    ]
    poly2 = [
        (1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0)
    ]
    points = intersect_convex(poly1, poly2, alg)
    @test issetequal(points, poly2)
    points = intersect_convex(poly2, poly1, alg)
    @test issetequal(points, poly2)
end

@testset "rectangles" begin
    poly1 = [
        (0.0, 0.0), (0.0, 2.0), (2.0, 2.0), (2.0, 0.0)
    ]
    points = intersect_convex(poly1, poly1, alg)
    @test issetequal(points, poly1)

    # overlap
    poly2 = translate(poly1, (1.0, 1.0))
    expected = [
        (1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0), 
    ]
    points = intersect_convex(poly1, poly2, alg)
    @test issetequal(points, expected)

    # no intersection
    poly2 = translate(poly1, (3.0, 3.0))
    expected = []
    points = intersect_convex(poly1, poly2, alg)
    @test issetequal(points, expected)
end

@testset "rectangles - edge overlap" begin
    poly1 = [
        (0.0, 0.0), (0.0, 2.0), (2.0, 2.0), (2.0, 0.0)
    ]

    poly2 = translate(poly1, (2.0, 1.0))
    expected = [
        (2.0, 1.0), (2.0, 2.0)
    ]
    points = intersect_convex(poly1, poly2, alg)
    @test issetequal(points, expected)
    points = intersect_convex(poly2, poly1, alg)
    @test issetequal(points, expected)

    poly2 = [
        (2.0, -1.0), (2.0, 3.0), (4.0, 3.0), (4.0, -1.0)
    ]
    expected = [
        (2.0, 0.0), (2.0, 2.0)
    ]
    points = intersect_convex(poly1, poly2, alg)
    @test issetequal(points, expected)
    points = intersect_convex(poly2, poly1, alg)
    @test issetequal(points, expected)

    # edge + vertix overlap
    poly2 = translate(poly1, (1.0, 0.0))
    expected = [
        (1.0, 2.0), (2.0, 2.0), (2.0, 0.0), (1.0, 0.0)
        ]
    points = intersect_convex(poly1, poly2, alg)
    @test issetequal(points, expected)
    points = intersect_convex(poly2, poly1, alg)
    @test issetequal(points, expected)
end

@testset "rectangles - vertix overlap" begin
    poly1 = [
        (0.0, 0.0), (0.0, 2.0), (2.0, 2.0), (2.0, 0.0)
    ]
    # single point of itersection
    poly2 = translate(poly1, (2.0, 2.0))
    expected = [(2.0, 2.0)]
    points = intersect_convex(poly1, poly2, alg)
    @test issetequal(points, expected)

    points = intersect_convex(poly2, poly1, alg)
    @test issetequal(points, expected)
end

@testset "edge intersect vertix inner" begin 
    poly1 = [
        (0.0, 0.0), (0.0, 3.0), (2.0, 5.0), (5.0, 0.0)
    ]
    poly2 = [
        (0.0, 2.0), (2.0, 3.0), (6.0, 1.0)
    ];
    expected = [
        (0.0, 2.0),
        (3.714286, 2.142857),
        (2.0, 3.0),
        (4.222222, 1.296296),
    ]
    points = intersect_convex(poly1, poly2, alg)
    answer = Set([round.(p, digits=6) for p in points])
    @test issetequal(answer, expected)

    points = intersect_convex(poly2, poly1, alg)
    answer = Set([round.(p, digits=6) for p in points])
    @test issetequal(answer, expected)
end

@testset "edge-vertix pass through" begin 
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
    points = intersect_convex(poly1, poly2, alg)
    answer = Set([round.(p, digits=6) for p in points])
    @test issetequal(answer, expected)

    points = intersect_convex(poly2, poly1, alg)
    answer = Set([round.(p, digits=6) for p in points])
    @test issetequal(answer, expected)
end

@testset "vertix intersect vertix" begin 
    poly1 = [
        (0.0, 0.0), (0.0, 3.0), (2.0, 5.0), (5.0, 0.0)
    ]
    poly2 = [
        (0.0, 0.0), (2.0, 3.0), (6.0, 1.0)
    ];
    expected = [
        (0.0, 0.0),
        (3.714286, 2.142857),
        (2.0, 3.0),
        (4.545455, 0.757576),
    ]
    points = intersect_convex(poly1, poly2, alg)
    answer = Set([round.(p, digits=6) for p in points])
    @test issetequal(answer, expected)

    points = intersect_convex(poly2, poly1, alg)
    answer = Set([round.(p, digits=6) for p in points])
    @test issetequal(answer, expected)
end

@testset "vertix intersections" begin 
    poly1 = [
        (0.0, 0.0), (0.5, 1.0), (1.0, 0.0)
    ]
    poly2 = [
        (0.0, 1.0), (0.0, 3.0), (1.0, 3.0), (1.0, 1.0)
    ]

    # single point
    expected = [(0.5, 1.0)]
    points = intersect_convex(poly1, poly2, alg)
    @test issetequal(points, expected)
    points = intersect_convex(poly2, poly1, alg)
    @test issetequal(points, expected)

    # 2 points inside
    poly1_ = translate(poly1, ((0.0), (1.5)))
    expected = poly1_
    points = intersect_convex(poly1_, poly2, alg)
    @test issetequal(points, expected)
    points = intersect_convex(poly2, poly1_, alg)
    @test issetequal(points, expected)

    # 3 points inside
    poly1_ = translate(poly1, ((0.0), (2.0)))
    expected = poly1_
    points = intersect_convex(poly1_, poly2, alg)
    @test issetequal(points, expected)
    points = intersect_convex(poly2, poly1_, alg)
    @test issetequal(points, expected)
end

@testset "quads single vertix intersect" begin 
    poly1 = [
        (1.0, 1.0), (2.0, 4.0), (5.0, 5.0), (4.0, 2.0)
    ]
    poly2 = [
        (4.0, 2.0), (3.0, 4.0), (6.0, 4.0), (7.0, 1.0)
    ]

    expected = [(4.0, 2.0), (3.0, 4.0), (4.666666666666667, 4.0)]
    points = intersect_convex(poly1, poly2, alg)
    @test issetequal(points, expected)
    points = intersect_convex(poly2, poly1, alg)
    @test issetequal(points, expected)
end

@testset "share lines" begin 
    poly1 = [
        (5.0, 8.0), (9.0, 4.0), (1.0, 4.0),
    ]
    poly2 = [
        (5.0, 1.0), (1.0, 4.0), (9.0, 4.0)
    ]
    expected = [(1.0, 4.0), (9.0, 4.0)]
    points = intersect_convex(poly1, poly2, alg)
    @test issetequal(points, expected)
    points = intersect_convex(poly2, poly1, alg)
    @test issetequal(points, expected)

    # shared point
    poly1 = [
        (5.0, 8.0), (9.0, 4.0), (5.0, 4.0), (1.0, 4.0),
    ]
    poly2 = [
        (5.0, 1.0), (1.0, 4.0), (5.0, 4.0), (9.0, 4.0)
    ]
    expected = [(1.0, 4.0), (5.0, 4.0), (9.0, 4.0), (5.0, 4.0),]
    points = intersect_convex(poly1, poly2, alg)
    @test issetequal(points, expected)
    points = intersect_convex(poly2, poly1, alg)
    @test issetequal(points, expected)
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
    points = intersect_convex(poly1, poly2, alg)
    @test issetequal(points, expected)

    points = intersect_convex(poly2, poly1, alg)
    @test issetequal(points, expected)
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
        (-2/3, 0.0), (2/3, 0.0), 
        (-1/3, h/3), (+1/3, h/3), 
        (-1/3, -h/3), (+1/3, -h/3), 
    ]
    points = intersect_convex(poly1, poly2, alg)
    @test issetequal(points, expected)

    points = intersect_convex(poly2, poly1, alg)
    @test issetequal(points, expected)
end

@testset "multiple pockets + edge-vertix" begin 
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
        (4.695652, 1.521739),
        (3.058824, 2.176471),
        (8.0, 4.0),
        (2.0, 6.0),
        (3.525424, 7.016949),
        (5.2, 7.5),
        (5.5, 7.5),
        (8.0, 5.5),
        (1.0, 5.5),
        (6.0, 7.0),
        (1.0, 4.75),
    ]
    points = intersect_convex(poly1, poly2, alg)
    answer = Set([round.(p, digits=6) for p in points])
    @test issetequal(answer, expected)

    points = intersect_convex(poly2, poly1, alg)
    answer = Set([round.(p, digits=6) for p in points])
    @test issetequal(answer, expected)
end

end