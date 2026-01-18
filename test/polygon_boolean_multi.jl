using PolygonAlgorithms
using PolygonAlgorithms: MartinezRuedaAlg, PointSet

@testset "polygon boolean multi - $alg" for alg in [
    MartinezRuedaAlg(),
]

@testset verbose=true "single" begin
    rect1 = [
        (0.0, 0.0), (0.0, 2.0), (1.0, 2.0), (1.0, 0.0)
    ]
    points = union_geometry([rect1], Vector{Tuple{Float64, Float64}}[], alg)
    @test PointSet(points[1]) == PointSet(rect1)
end

@testset verbose=true "self-intersecting star" begin
    self_intersect_star = [
        (-3.0, 2.0), (3.0, 2.0), (-2.0, -2.0), (0.0, 5.0), (2.0, -2.0)
    ]
    expected = [
        [(0.0, -0.4), (-2.0, -2.0), (-1.255814, 0.604651), (-3.0, 2.0), (-0.857143, 2.0), (-0.857143, 2.0), (-1.255814, 0.604651)],
        [(0.857143, 2.0), (-0.857143, 2.0), (0.0, 5.0)],
        [(3.0, 2.0), (1.255814, 0.604651), (2.0, -2.0), (0.0, -0.4), (1.255814, 0.604651), (0.857143, 2.0)],
    ]
    regions = union_geometry([self_intersect_star], Vector{Tuple{Float64, Float64}}[], alg)
    @test are_regions_equal(regions, expected)
end

@testset verbose=true "cross + triangle" begin
    rect1 = [
        (0.0, 0.0), (0.0, 2.0), (1.0, 2.0), (1.0, 0.0)
    ]
    rect2 = [
        (-1.0, 1.5), (2.0, 1.5), (2.0, 0.5), (-1.0, 0.5)
    ]
    triangle = [
        (-0.5, 1.5), (1.5, 2.0), (0.5, 0.25)
    ]
    @testset "intersect" begin
        expected = [
            [(1.0, 1.5), (1.0, 1.125), (0.642857, 0.5), (0.3, 0.5), (0.0, 0.875), (0.0, 1.5)]
        ]
        points = intersect_geometry([rect1], [rect2, triangle], alg)
        @test PointSet(points[1]) == PointSet(expected[1])

        points2 = intersect_geometry(rect1, rect2, alg)
        point2 = intersect_geometry(points2[1], triangle, alg)
        @test points == point2
    end

    @testset "union" begin
        expected = [
            [
                (2.0, 1.5), (2.0, 0.5), (1.0, 0.5), (1.0, 0.0), (0.0, 0.0), (0.0, 0.5), (0.0, 0.5),
                (-1.0, 0.5), (-1.0, 1.5), (-0.5, 1.5), (0.0, 1.625), (0.0, 2.0), (1.0, 2.0),
                (1.0, 2.0), (1.0, 1.875), (1.5, 2.0), (1.5, 2.0), (1.2142857142857142, 1.5)
            ]
        ]
        points = union_geometry([rect1], [rect2, triangle], alg)
        @test PointSet(points[1]) == PointSet(expected[1])
        points2 = union_geometry(rect1, rect2, alg)
        point2 = union_geometry(points2[1], triangle, alg)
        @test points == point2
    end

    @testset "difference" begin
        expected = [
            [(1.0, 0.5), (1.0, 0.0), (0.0, 0.0), (0.0, 0.5), (0.3, 0.5), (0.5, 0.25), (0.642857, 0.5)], # top
            [(1.0, 2.0), (1.0, 1.875), (0.0, 1.625), (0.0, 2.0)], # bottom
        ]
        regions = difference_geometry([rect1], [rect2, triangle], alg)
        @test are_regions_equal(regions, expected)
        regions = difference_geometry(rect1, rect2, alg)
        regions2 = difference_geometry(regions, [triangle], alg)
        @test are_regions_equal(regions2, expected)
        regions2 = difference_geometry(regions[1:1], vcat([regions[2]], [triangle])) # top region is no longer primary
        @test are_regions_equal(regions2, expected[1:1]) # not a bug. The top region was excluded.

        # combine rect1 + triangle. They overlap so equivalent to self-intersecting with holes
        regions3 = difference_geometry([rect1, triangle], [rect2])
        expected = [
            [(-0.5, 1.5), (0.0, 1.5), (0.0, 1.625)],
            [(1.0, 0.5), (1.0, 0.0), (0.0, 0.0), (0.0, 0.5), (0.3, 0.5), (0.5, 0.25), (0.6428571428571428, 0.5)],
            [(1.5, 2.0), (1.2142857142857142, 1.5), (1.0, 1.5), (1.0, 1.875), (1.0, 1.875), (0.0, 1.625), (0.0, 2.0), (1.0, 2.0), (1.0, 2.0), (1.0, 1.875)],
        ]
        @test are_regions_equal(regions3, expected)
    end

    @testset "xor" begin
        expected = [
            [(0.3, 0.5), (0.5, 0.25), (0.642857, 0.5)], # hole
            [
                (2.0, 1.5), (2.0, 0.5), (1.0, 0.5), (1.0, 1.125), (1.0, 1.125),
                (0.642857, 0.5), (1.0, 0.5), (1.0, 0.0), (0.0, 0.0), (0.0, 0.5),
                (0.0, 0.5), (-1.0, 0.5), (-1.0, 1.5), (-0.5, 1.5), (0.0, 0.875), (0.0, 0.875),
                (0.0, 0.5), (0.3, 0.5), (0.0, 0.875), (0.0, 0.875), (0.0, 1.5), (-0.5, 1.5),
                (0.0, 1.625), (0.0, 1.625), (0.0, 1.5), (1.0, 1.5), (1.0, 1.5), (1.0, 1.125),
                (1.2142857142857142, 1.5), (1.0, 1.5), (1.0, 1.5), (1.0, 1.875), (1.0, 1.875),
                (0.0, 1.625), (0.0, 2.0), (1.0, 2.0), (1.0, 2.0), (1.0, 1.875), (1.5, 2.0),
                (1.5, 2.0), (1.2142857142857142, 1.5)
            ]
        ]
        regions = xor_geometry([rect1], [rect2, triangle], alg)
        @test are_regions_equal(regions, expected)
        regions = xor_geometry(rect1, rect2, alg)
        regions = xor_geometry(regions, [triangle], alg)
        @test are_regions_equal(regions, expected)
    end
end

@testset verbose=true "appearing holes" begin
    elbow = [
        (0.0, -1.0), (0.0, 2.0), (3.0, 2.0), (3.0, 1.0), (1.0, 1.0), (1.0, -1.0)
    ]
    triangle = [
        (0.8, -0.8), (2.5, 1.8), (2.5, -0.8)
    ]
    rect = [(0.5, 0.0), (0.5, 0.5), (2.0, 0.5), (2.0, 0.0)]

    @testset "intersect" begin
        regions = intersect_geometry([elbow], [triangle, rect], alg)
        @test isempty(regions)

        regions = intersect_geometry([elbow], [triangle], alg)
        regions = intersect_geometry(regions, [rect], alg)
        @test isempty(regions)
    end

    @testset "union" begin
        expected = [
            [(3.0, 2.0), (3.0, 1.0), (2.5, 1.0), (2.5, -0.8), (1.0, -0.8), (1.0, -1.0), (0.0, -1.0), (0.0, 2.0)],
            [(1.323077, 0.0), (1.0, -0.4941176), (1.0, -0.0)], # hole
            [(1.976923, 1.0), (1.65, 0.5), (1.0, 0.5), (1.0, 1.0)], # hole
        ]
        regions = union_geometry([elbow], [triangle, rect], alg)
        @test are_regions_equal(regions, expected)

        regions = union_geometry(elbow, triangle, alg)
        regions = union_geometry(regions, [rect], alg)
        @test are_regions_equal(regions, expected)
    end

    @testset "difference" begin
        expected = [
            [
                (3.0, 2.0), (3.0, 1.0), (2.5, 1.0), (2.5, 1.8), (2.5, 1.8),
                (1.976923, 1.0), (1.0, 1.0), (1.0, 0.5), (0.5, 0.5),
                (0.5, 0.0), (0.5, 0.0), (1.0, 0.0), (1.0, -0.4941176),
                (0.8, -0.8), (1.0, -0.8), (1.0, -1.0), (0.0, -1.0), (0.0, 2.0),
            ]
        ]
        regions = difference_geometry([elbow], [triangle, rect], alg)
        @test are_regions_equal(regions, expected)

        regions = difference_geometry(elbow, triangle, alg)
        regions = difference_geometry(regions, [rect], alg)
        @test are_regions_equal(regions, expected)
    end

    @testset "xor" begin
        expected = [
            [(3.0, 2.0), (3.0, 1.0), (2.5, 1.0), (2.5, 1.8), (2.5, 1.8), (1.976923, 1.0), (2.5, 1.0), (2.5, -0.8), (1.0, -0.8), (1.0, -0.4941176), (1.0, -0.4941176), (0.8, -0.8), (1.0, -0.8), (1.0, -1.0), (0.0, -1.0), (0.0, 2.0)],
            # holes:
            [(1.3230769, 0.0), (1.0, -0.4941176), (1.0, -0.0), (1.0, -0.0), (0.5, 0.0), (0.5, 0.5), (1.0, 0.5), (1.0, 0.5), (1.0, -0.0)],
            [(2.0, 0.5), (2.0, 0.0), (1.3230769, 0.0), (1.65, 0.5), (1.65, 0.5), (1.0, 0.5), (1.0, 1.0), (1.976923, 1.0), (1.976923, 1.0), (1.65, 0.5)],
        ]
        regions = xor_geometry([elbow], [triangle, rect], alg)
        @test are_regions_equal(regions, expected)

        regions = xor_geometry(elbow, triangle, alg)
        regions = xor_geometry(regions, [rect], alg)
        @test are_regions_equal(regions, expected)
    end
end

end