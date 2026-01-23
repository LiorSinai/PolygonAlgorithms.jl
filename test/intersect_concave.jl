using PolygonAlgorithms: translate, PointSet
using PolygonAlgorithms: WeilerAthertonAlg, MartinezRuedaAlg

@testset "intersections concave - $alg" for alg in [
    WeilerAthertonAlg(),
    MartinezRuedaAlg(),
]

@testset "rectangle jagged" begin 
    poly1 = [
        (0.0, 0.0), (0.0, 4.0), (3.0, 4.0), (3.0, 0.0)
    ]
    poly2 = [
        (2.0, 3.0), (5.0, 3.0), (5.0, 1.0), (2.0, 1.0), (4.0, 2.0)
    ];
    expected = [
        [(3.0, 1.0), (2.0, 1.0), (3.0, 1.5)],
        [(3.0, 2.5), (2.0, 3.0), (3.0, 3.0)]
    ]

    regions = intersect_geometry(alg, poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(alg, poly2, poly1)
    @test are_regions_equal(regions, expected)

    # now vertex intersects edge
    # creates cycle
    poly2_ = translate(poly2, (-1.0, 0.0))
    expected = (typeof(alg) == PolygonAlgorithms.MartinezRuedaAlg) ?
        [ [(1.0, 1.0), (3.0, 1.0), (3.0, 2.0)], [(1.0, 3.0), (3.0, 2.0), (3.0, 3.0)]] :
        [[(3.0, 1.0), (1.0, 1.0), (3.0, 2.0),(1.0, 3.0), (3.0, 3.0)]]

    regions = intersect_geometry(alg, poly1, poly2_,)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(alg, poly2_, poly1)
    @test are_regions_equal(regions, expected)

    # only points intercept
    poly2_ = translate(poly2, (1.0, 0.0))
    expected = (typeof(alg) == PolygonAlgorithms.MartinezRuedaAlg) ?
        Vector{Tuple{Float64, Float64}}[] :
        [[(3.0, 3.0)], [(3.0, 1.0)]]
    regions = intersect_geometry(alg, poly1, poly2_)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(alg, poly2_, poly1)
    @test are_regions_equal(regions, expected)
end


@testset "concave Cs" begin 
    poly1 = [
        (0.0, 0.0), (4.0, 0.0), (4.0, 6.0), (0.0, 6.0), (0.0, 4.0), (2.0, 4.0), (2.0, 2.0), (0.0, 2.0)
    ]
    poly2 = [
        (-2.0, 0.5), (5.0, 0.5), (5.0, 1.5), (-1.0, 1.5),
        (-1.0, 4.5), (3.0, 4.5), (3.0, 5.5), (-2.0, 5.5)
    ];
    expected = [
        [(0.0, 5.5), (3.0, 5.5), (3.0, 4.5), (0.0, 4.5)],
        [(0.0, 1.5), (4.0, 1.5), (4.0, 0.5), (0.0, 0.5)],
    ]

    regions = intersect_geometry(alg, poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(alg, poly2, poly1)
    @test are_regions_equal(regions, expected)
end

@testset "concave <>" begin 
    poly1 = [
        (0.0, 0.0), (-1.0, 1.0),  (4.0, 6.0), (-1.0, 11.0), (0.0, 12.0), (6.0, 6.0),  
    ]
    poly2 = [
        (2.0, 0.0), (3.0, 1.0),  (-2.0, 6.0), (3.0, 11.0), (2.0, 12.0), (-4.0, 6.0),  
    ]
    expected = [
        [(0.0, 10.0), (1.0, 11.0), (2.0, 10.0), (1.0, 9.0)],
        [(1.0, 3.0), (2.0, 2.0), (1.0, 1.0), (0.0, 2.0)],
    ]

    regions = intersect_geometry(alg, poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(alg, poly2, poly1)
    @test are_regions_equal(regions, expected)
end

@testset "concave arrows vertex intercepts" begin 
    poly1 = [
        (-2.0, 2.0), (2.0, 2.0), (-1.0, 0.0), (-0.5, 1.5)
    ]
    poly2 = [
        (-1.0, 3.0), (2.0, 0.0), (-2.0, 1.0), (-0.5, 1.5) 
    ]

    expected = [[
        (-0.666667, 2.0)
        (0.0, 2.0)
        (0.8, 1.2)
        (-0.181818, 0.545455)
        (-0.769231, 0.692308)
        (-0.5, 1.5)
    ]]
    regions = intersect_geometry(alg, poly1, poly2)
    regions = [[round.(p, digits=6) for p in r] for r in regions]
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(alg, poly2, poly1)
    regions = [[round.(p, digits=6) for p in r] for r in regions]
    @test are_regions_equal(regions, expected)
end

@testset "concave share edge inner" begin 
    poly1 = [
        (0.0, 0.0), (0.0, 1.0), (-1.0, 1.0), (-1.0, 2.0), (1.0, 2.0), (1.0, 0.0)
    ]
    poly2 = [
        (1.5, 0.5), (0.0, 0.5), (0.0, 1.0), (-0.5, 1.0), (0.5, 1.5) 
    ]

    expected = [[
        (0.0, 0.5),
        (0.0, 1.0),
        (-0.5, 1.0),
        (0.5, 1.5),
        (1.0, 1.0),
        (1.0, 0.5),
    ]]
    regions = intersect_geometry(alg, poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(alg, poly2, poly1)
    @test are_regions_equal(regions, expected)
end

@testset "concave share edge outer" begin 
    poly1 = [
        (0.0, 0.0), (0.0, 1.0), (-1.0, 1.0), (-1.0, 2.0), (1.0, 2.0), (1.0, 0.0)
    ]
    poly2 = [
        (-2.0, 1.0), (0.0, 1.0), (0.0, -1.0)
    ]

    expected = [[(0.0, 0.0), (0.0, 1.0), (-1.0, 1.0), (0.0, 1.0)]]
    regions = intersect_geometry(alg, poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(alg, poly2, poly1)
    @test are_regions_equal(regions, expected)
end

@testset "stars - numeric test" begin 
    poly1 = [
        (0.0, 18.0), (3.0, 5.0), (15.0, 5.0), (5.0, 0.0), (10.0, -12.0), (0.0, -2.0),
        (-10.0, -12.0), (-5.0, 0.0), (-15.0, 5.0), (-3.0, 5.0)
    ]
    ;
    poly2 = PolygonAlgorithms.rotate(poly1, π/1.0, (0.0, 0.0));
    expected = [[
        (-3.0, -5.0), (-7.083333333333333, -5.0), (-5.0, 0.0), (-7.083333333333333, 5.0), (-3.0, 5.0), 
        (0.0, 2.0), (3.0, 5.0), (7.083333333333333, 5.0), (5.0, 0.0), 
        (7.083333333333333, -5.0), (3.0, -5.0), (0.0, -2.0), (-3.0, -5.0)
    ]]

    regions = intersect_geometry(alg, poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(alg, poly2, poly1)
    @test are_regions_equal(regions, expected)
end

end

@testset "intersections concave - only WeilerAthertonAlg()" begin
#=
    The other algorithms can handle these cases, but may give slightly different results.
    E.g. exclude line/point intersections.
=#
alg = WeilerAthertonAlg()

@testset "concave outer share portion" begin 
    poly1 = [
        (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 2.0), (2.0, 2.0), (2.0, 0.0)
    ]
    poly2 = [
        (-1.0, 0.0), (1.0, 1.0), (0.0, -1.0)
    ]

    expected = [[(1.0, 1.0)]]
    regions = intersect_geometry(alg, poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(alg, poly2, poly1)
    @test are_regions_equal(regions, expected)

    # share edge to the left
    poly2 = [
        (-1.0, 0.0), (0.8, 1.0), (1.0, 1.0), (0.0, -1.0)
    ]
    expected = [[(0.8, 1.0), (1.0, 1.0)]]
    regions = intersect_geometry(alg, poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(alg, poly2, poly1)
    @test are_regions_equal(regions, expected)

    # share edge to the right
    poly2 = [
        (-1.0, 0.0), (1.0, 1.0), (1.0, 0.9), (0.0, -1.0)
    ]
    expected = [[(1.0, 1.0), (1.0, 0.9)]]
    regions = intersect_geometry(alg, poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(alg, poly2, poly1)
    @test are_regions_equal(regions, expected)

    # share 2 edges
    poly2 = [
        (-1.0, 0.0), (0.8, 1.0), (1.0, 1.0), (1.0, 0.9), (0.0, -1.0)
    ]
    expected = [[(0.8, 1.0), (1.0, 1.0), (1.0, 0.9)]]
    regions = intersect_geometry(alg, poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(alg, poly2, poly1)
    @test are_regions_equal(regions, expected)
end

@testset "concave convex share outer" begin 
    poly1 = [
        (0.0, 0.0), (0.0, 1.0), (1.0, 0.0)
    ]
    poly2 = [
        (2.0, 1.0), (2.0, -0.5), (1.0, 0.0),
    ]

    expected = [[(1.0, 0.0)]]
    regions = intersect_geometry(alg, poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(alg, poly2, poly1)
    @test are_regions_equal(regions, expected)

    # share edge below
    poly2 = [
        (2.0, 1.0), (2.0, -0.5), (0.5, 0.0), (1.0, 0.0)
    ]
    expected = [[(0.5, 0.0), (1.0, 0.0)]]
    regions = intersect_geometry(alg, poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(alg, poly2, poly1)
    @test are_regions_equal(regions, expected)

    # share edge above
    poly2 = [
        (2.0, 1.0), (2.0, -0.5), (-0.5, -0.5), (1.0, 0.0), (0.5, 0.5)
    ]
    expected = [[(1.0, 0.0), (0.5, 0.5)]]
    regions = intersect_geometry(alg, poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(alg, poly2, poly1)
    @test are_regions_equal(regions, expected)

    # share both
    poly2 = [
        (2.0, 1.0), (2.0, -0.5), (0.5, 0.0), (1.0, 0.0), (0.5, 0.5)
    ]
    expected = [[(0.5, 0.0), (1.0, 0.0), (0.5, 0.5)]]
    regions = intersect_geometry(alg, poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(alg, poly2, poly1)
    @test are_regions_equal(regions, expected)
end

@testset "concave outer share portion and inner" begin 
    # self-intersecting area
    poly1 = [
        (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 2.0), (2.0, 2.0), (2.0, 0.0)
    ]
    poly2 = [
        (-1.0, 0.0), (0.3, 1.5), (0.7, 1.0), (1.0, 1.0), (0.0, -1.0)
    ]
    expected = [[(1.0, 1.0), (0.7, 1.0), (0.0, 1.0), (0.0, 1.1538461538461537), (0.3, 1.5), (0.7, 1.0), (1.0, 1.0)]]
    regions = intersect_geometry(alg, poly1, poly2)
    @test_broken are_regions_equal(regions, expected)
    regions = intersect_geometry(alg, poly2, poly1)
    @test are_regions_equal(regions, expected)
end

@testset "concave saw + vertex intercepts" begin 
    poly1 = [
        (0.0, 1.0), (0.0, 2.0), (3.0, 2.0), (3.0, 1.0)
    ]
    poly2 = [
        (0.0, -0.5), (0.0, 0.0), (0.5, 1.0), (1.0, 0.0), (1.5, 1.0), (2.0, 0.0), (2.0, 2.0), (3.0, -0.5)
    ]

    expected = [
        [(2.0, 1.0), (2.0, 2.0), (2.4, 1.0)],
        [(0.5, 1.0)],
        [(1.5, 1.0)]
    ]
    
    regions = intersect_geometry(alg, poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(alg, poly2, poly1)
    @test are_regions_equal(regions, expected)
end

@testset "overlapping arches" begin
    poly1 = [
        (0.0, 0.0), (-1.0, 0.0), (-1.0, 2.0), (2.0, 2.0), (2.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)
    ]
    poly2 = [
        (0.0, 1.0), (0.0, 2.0), (2.0, 2.0), (2.0, -1.0), (0.0, -1.0), (0.0, 0.0),  (1.0, 0.0), (1.0, 1.0),
    ]
    expected = [
        [(1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 2.0), (2.0, 2.0), (2.0, 0.0)],
        [(0.0, 0.0)],
    ]
    regions = intersect_geometry(alg, poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(alg, poly2, poly1)
    @test are_regions_equal(regions, expected)

    # mirror image
    poly1 = [
        (0.0, 0.0), (1.0, 0.0), (1.0, 2.0), (-2.0, 2.0), (-2.0, 0.0), (-1.0, 0.0), (-1.0, 1.0), (0.0, 1.0)
    ]
    poly2 = [
        (0.0, 1.0), (0.0, 2.0), (-2.0, 2.0), (-2.0, -1.0), (0.0, -1.0), (0.0, 0.0),  (-1.0, 0.0), (-1.0, 1.0),
    ]
    expected = [
        [(-2.0, 0.0), (-2.0, 2.0), (0.0, 2.0), (0.0, 1.0), (-1.0, 1.0), (-1.0, 0.0),],
        [(0.0, 0.0)],
    ]
    regions = intersect_geometry(alg, poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(alg, poly2, poly1)
    @test are_regions_equal(regions, expected)
end


@testset "hilbert curve - order 2" begin 
    # This is an example of intersection regions which are self-intersecting
    # despite the original polygons not being self-intersecting.
    poly1 = [
        (-0.125, 0.125), (-0.125, 0.875), (0.125, 0.875), (0.125, 0.625),
        (0.375, 0.625), (0.375, 0.875), (0.625, 0.875), (0.875, 0.875),
        (0.875, 0.625), (0.625, 0.625), (0.625, 0.375), (0.875, 0.375),
        (0.875, 0.125), (0.625, 0.125), (0.375, 0.125), (0.375, 0.375),
        (0.125, 0.375), (0.125, 0.125),
    ]
    poly2 = PolygonAlgorithms.rotate(poly1, π/2.0, (0.5, 0.5))

    expected= [
        [
            (0.125, 0.875), (0.125, 0.625), 
            (0.375, 0.625), (0.375, 0.875), (0.375, 0.625), 
            (0.625, 0.625), (0.625, 0.875), (0.875, 0.875), (0.875, 0.625), (0.625, 0.625),
            (0.625, 0.375), (0.875, 0.375), (0.625, 0.375),
            (0.625, 0.125), (0.875, 0.125), (0.625, 0.125),
            (0.375, 0.125), (0.375, 0.375), (0.125, 0.375), (0.125, 0.625),
        ],
        [(0.125, 0.125)]
    ]
    regions = intersect_geometry(alg, poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(alg, poly2, poly1)
    @test are_regions_equal(regions, expected)
end

end

@testset "intersections concave - only MartinezRuedaAlg()" begin
#=
    Only the MartinezRuedaAlg can handle self intersecting polygons.
=#
alg = MartinezRuedaAlg()

@testset "self-intersect rectangle" begin
    self_intersect = [
        (0.0, 0.0), (2.0, 2.0), (6.0, -2.0), (11.0, 2.0), (11.0, 0.0)
    ]
    rectangle_horiz = [
        (-1.0, 0.0), (-1.0, 3.0), (12.0, 3.0), (12.0, 0.0)
    ];
    poly1 = self_intersect
    poly2 = rectangle_horiz
    expected = [
        [(4.0, 0.0), (0.0, 0.0), (2.0, 2.0)],
        [(8.5, -0.0), (11.0, 0.0), (11.0, 2.0)],
        [(4.0, 0.0), (8.5, -0.0)], # straight line
    ]
    regions = intersect_geometry(alg, poly1, poly2)
    @test are_regions_equal(regions, expected)
end

@testset "self-intersect star" begin
    self_intersect_star = [
        (-3.0, 2.0), (3.0, 2.0), (-2.0, -2.0), (0.0, 5.0), (2.0, -2.0)
    ]
    box = [
        (-2.0, 3.0), (2.0, 3.0), (2.0, -1.0), (-2.0, -1.0)
    ];
    poly1 = self_intersect_star
    poly2 = box
    expected= [
        [(-0.0, -0.4), (-0.75, -1.0), (-1.7142857142857144, -1.0), (-1.255813953488372, 0.6046511627906979), (-1.255813953488372, 0.6046511627906979), (-2.0, 1.2000000000000002), (-2.0, 2.0), (-0.857142857142857, 2.0), (-0.857142857142857, 2.0), (-1.255813953488372, 0.6046511627906979)],
        [(0.8571428571428572, 2.0), (-0.857142857142857, 2.0), (-0.5714285714285712, 3.0), (0.5714285714285715, 2.9999999999999996)],
        [(0.8571428571428572, 2.0), (1.255813953488372, 0.6046511627906976), (-0.0, -0.4), (0.7499999999999999, -1.0), (1.7142857142857142, -1.0), (1.255813953488372, 0.6046511627906976), (2.0, 1.2000000000000002), (2.0, 2.0)],
    ]
    regions = intersect_geometry(alg, poly1, poly2)
    @test are_regions_equal(regions, expected)
end

end