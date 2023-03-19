@testset verbose = true "intersect-concave" begin

function are_regions_equal(r1::Vector{Vector{T}}, r2::Vector{Vector{T}}) where T
    r1_sets = [Set(r) for r in r1]
    r2_sets = [Set(r) for r in r2]
    issetequal(r1_sets, r2_sets)
end

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

    regions = intersect_geometry(poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(poly2, poly1)
    @test are_regions_equal(regions, expected)

    # now vertix intersects edge
    # creates cycle
    poly2_ = translate(poly2, (-1.0, 0.0))
    expected = [
        [(3.0, 1.0), (1.0, 1.0), (3.0, 2.0),(1.0, 3.0), (3.0, 3.0)]
    ]

    regions = intersect_geometry(poly1, poly2_)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(poly2_, poly1)
    @test are_regions_equal(regions, expected)

    # only points intercept
    poly2_ = translate(poly2, (1.0, 0.0))
    expected = [
        [(3.0, 3.0)], [(3.0, 1.0)]
    ]
    regions = intersect_geometry(poly1, poly2_)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(poly2_, poly1)
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

    regions = intersect_geometry(poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(poly2, poly1)
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

    regions = intersect_geometry(poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(poly2, poly1)
    @test are_regions_equal(regions, expected)
end

@testset "concave arrows vertix intercepts" begin 
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
    regions = intersect_geometry(poly1, poly2)
    regions = [[round.(p, digits=6) for p in r] for r in regions]
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(poly2, poly1)
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
    regions = intersect_geometry(poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(poly2, poly1)
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
    regions = intersect_geometry(poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(poly2, poly1)
    @test are_regions_equal(regions, expected)
end

@testset "concave outer vertix intercept" begin 
    poly1 = [
        (0.0, 0.0), (0.0, 1.0), (-1.0, 1.0), (-1.0, 2.0), (1.0, 2.0), (1.0, 0.0)
    ]
    poly2 = [
        (-2.0, 0.0), (0.0, 1.0), (-1.0, -1.0)
    ]

    expected = [[(0.0, 1.0)]]
    regions = intersect_geometry(poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(poly2, poly1)
    @test are_regions_equal(regions, expected)
end

@testset "concave saw + vertix intercepts" begin 
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
    
    regions = intersect_geometry(poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(poly2, poly1)
    @test are_regions_equal(regions, expected)
end

end