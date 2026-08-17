using PolygonAlgorithms
using PolygonAlgorithms: MartinezRuedaAlg, PointSet
using PolygonAlgorithms: MERGE_FACES, SPLIT_FACES

@testset "polygon boolean multi - $alg" for alg in [
    MartinezRuedaAlg(),
]

@testset verbose=true "single" begin
    rect1 = [
        (0.0, 0.0), (0.0, 2.0), (1.0, 2.0), (1.0, 0.0)
    ]
    points = union_geometry(alg, rect1)
    @test PointSet(points[1]) == PointSet(rect1)
end

@testset "overlap" begin
    poly1 = [(1.0, 2.0), (1.0, 4.0), (3.0, 4.0), (3.0, 2.0)];
    poly2 = [(1.0, 1.0), (2.0, 3.0), (3.0, 1.0)];
    # union two different polygons
    regions = union_geometry(alg, poly1, poly2)
    exterior = [[
        (3.0, 4.0), (3.0, 2.0), (2.5, 2.0), (3.0, 1.0), (1.0, 1.0), 
        (1.5, 2.0), (1.5, 2.0), (1.0, 2.0), (1.0, 4.0)
    ]]
    @test are_regions_equal(regions, exterior)
    # both are part of the subject
    regions = union_geometry(alg, [poly1, poly2])
    expected = [
        [(1.0, 1.0), (3.0, 1.0), (2.5, 2.0), (1.5, 2.0)],
        [(2.0, 3.0), (2.5, 2.0), (3.0, 2.0), (3.0, 4.0), (1.0, 4.0), (1.0, 2.0), (1.5, 2.0)]
    ]
    @test are_regions_equal(regions, expected)
end

@testset verbose=true "self-intersecting star" begin
    self_intersect_star = [
        (-3.0, 2.0), (3.0, 2.0), (-2.0, -2.0), (0.0, 5.0), (2.0, -2.0)
    ]
    regions = union_geometry(alg, self_intersect_star; face_selection=MERGE_FACES)
    expected = [
        # exterior
        [(0.0, -0.4), (2.0, -2.0), (1.255814, 0.604651), (3.0, 2.0), (0.857143, 2.0), (0.0, 5.0), (-0.857143, 2.0), (-3.0, 2.0), (-1.255814, 0.604651), (-2.0, -2.0)],
        # hole
        [(-0.857143, 2.0), (0.857143, 2.0), (1.255814, 0.604651), (0.0, -0.4), (-1.255814, 0.604651)],
    ]
    @test are_regions_equal(regions, expected)
    regions = union_geometry(alg, self_intersect_star; face_selection=SPLIT_FACES)
    expected = [
        [(-0.857143, 2.0), (0.857143, 2.0), (0.0, 5.0)],
        [(0.857143, 2.0), (1.255814, 0.604651), (3.0, 2.0)],
        [(-1.255814, 0.604651), (-2.0, -2.0), (0.0, -0.4)],
        [(-0.857143, 2.0), (-3.0, 2.0), (-1.255814, 0.604651)],
        [(0.0, -0.4), (2.0, -2.0), (1.255814, 0.604651)],
    ]
    @test are_regions_equal(regions, expected)
    inner = [(-1.0, 1.0), (0.0, 1.5), (1.0, 1.0)];
    regions = union_geometry(alg, self_intersect_star, inner; face_selection=SPLIT_FACES)
    expected2 = vcat(expected[1:3], [inner], expected[4:5])
    @test are_regions_equal(regions, expected2)
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
        points = intersect_geometry(alg, rect1, rect2, triangle)
        @test PointSet(points[1]) == PointSet(expected[1])

        points2 = intersect_geometry(alg, rect1, rect2)
        point2 = intersect_geometry(alg, points2[1], triangle)
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
        points = union_geometry(alg, rect1, rect2, triangle)
        @test PointSet(points[1]) == PointSet(expected[1])
        points2 = union_geometry(alg, rect1, rect2)
        point2 = union_geometry(alg, points2[1], triangle)
        @test points == point2
    end

    @testset "difference" begin
        expected = [
            [(1.0, 0.5), (1.0, 0.0), (0.0, 0.0), (0.0, 0.5), (0.3, 0.5), (0.5, 0.25), (0.642857, 0.5)], # top
            [(1.0, 2.0), (1.0, 1.875), (0.0, 1.625), (0.0, 2.0)], # bottom
        ]
        regions = difference_geometry(alg, rect1, rect2, triangle)
        @test are_regions_equal(regions, expected)
        regions = difference_geometry(alg, rect1, rect2)
        regions2 = difference_geometry(alg, regions, triangle)
        @test are_regions_equal(regions2, expected)
        regions2 = difference_geometry(alg, regions[2], regions[1], triangle) # top region is no longer primary
        @test are_regions_equal(regions2, [[(1.0, 2.0), (0.0, 2.0), (0.0, 1.625), (1.0, 1.875)]])

        # combine rect1 + triangle. They overlap so equivalent to self-intersecting with holes
        regions3 = difference_geometry(alg, [rect1, triangle], rect2; face_selection=MERGE_FACES)
        expected = [
             [(-0.5, 1.5), (0.0, 1.5), (0.0, 1.625), (1.0, 1.875), (1.0, 1.5), (1.214286, 1.5), (1.5, 2.0), (1.0, 1.875), (1.0, 2.0), (0.0, 2.0), (0.0, 1.625)],
             [(0.642857, 0.5), (0.5, 0.25), (0.3, 0.5), (0.0, 0.5), (0.0, 0.0), (1.0, 0.0), (1.0, 0.5)],
        ]
        @test are_regions_equal(regions3, expected)
        regions3 = difference_geometry(alg, [rect1, triangle], rect2)
        expected = [
            [(1.0, 1.875), (1.0, 2.0), (0.0, 2.0), (0.0, 1.625)],
            [(1.5, 2.0), (1.0, 1.875), (1.0, 1.5), (1.214286, 1.5)],
            [(-0.5, 1.5), (0.0, 1.5), (0.0, 1.625)],
            [(0.642857, 0.5), (0.5, 0.25), (0.3, 0.5), (0.0, 0.5), (0.0, 0.0), (1.0, 0.0), (1.0, 0.5)],
        ]
        @test are_regions_equal(regions3, expected)
    end

    @testset "xor" begin
        expected = [
            [(0.0, 1.625), (-0.5, 1.5), (0.0, 1.5)],
            [(1.0, 1.875), (1.0, 2.0), (0.0, 2.0), (0.0, 1.625)],
            [(0.0, 0.875), (-0.5, 1.5), (-1.0, 1.5), (-1.0, 0.5), (0.0, 0.5)],
            [(1.0, 0.5), (0.642857, 0.5), (0.5, 0.25), (0.3, 0.5), (0.0, 0.5), (0.0, 0.0), (1.0, 0.0)],
            [(1.0, 1.5), (1.214286, 1.5), (1.5, 2.0), (1.0, 1.875)],
            [(0.0, 0.875), (0.3, 0.5), (0.642857, 0.5), (1.0, 1.125), (1.0, 1.5), (0.0, 1.5)],
            [(1.214286, 1.5), (1.0, 1.125), (1.0, 0.5), (2.0, 0.5), (2.0, 1.5)],
        ]
        regions = xor_geometry(alg, rect1, rect2, triangle)
        @test are_regions_equal(regions, expected)
        regions = xor_geometry(alg, rect1, rect2)
        regions = xor_geometry(alg, regions, triangle)
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
        regions = intersect_geometry(alg, elbow, triangle, rect)
        @test isempty(regions)

        regions = intersect_geometry(alg, elbow, triangle)
        regions = intersect_geometry(alg, regions, rect)
        @test isempty(regions)
    end

    @testset "union" begin
        expected = [
            [(3.0, 2.0), (3.0, 1.0), (2.5, 1.0), (2.5, -0.8), (1.0, -0.8), (1.0, -1.0), (0.0, -1.0), (0.0, 2.0)],
            [(1.323077, 0.0), (1.0, -0.4941176), (1.0, -0.0)], # hole
            [(1.976923, 1.0), (1.65, 0.5), (1.0, 0.5), (1.0, 1.0)], # hole
        ]
        regions = union_geometry(alg, elbow, triangle, rect)
        @test are_regions_equal(regions, expected)

        regions = union_geometry(alg, elbow, triangle)
        regions = union_geometry(alg, regions, rect)
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
        regions = difference_geometry(alg, elbow, triangle, rect)
        @test are_regions_equal(regions, expected)

        regions = difference_geometry(alg, elbow, triangle)
        regions = difference_geometry(alg, regions, rect)
        @test are_regions_equal(regions, expected)
    end

    @testset "xor" begin
        expected = [
            [
                (0.5, 0.0), (0.5, 0.5), (1.0, 0.5), (1.0, 1.0), (1.976923, 1.0),
                (2.5, 1.8), (2.5, 1.0), (3.0, 1.0), (3.0, 2.0), (0.0, 2.0),
                (0.0, -1.0), (1.0, -1.0), (1.0, -0.8), (0.8, -0.8), (1.0, -0.494118),
                (1.0, 0.0)
            ],
            [(1.0, 0.0), (1.323077, 0.0), (1.65, 0.5), (1.0, 0.5)],
            [
                (1.323077, 0.0), (1.0, -0.494118), (1.0, -0.8), (2.5, -0.8),
                (2.5, 1.0), (1.976923, 1.0), (1.65, 0.5), (2.0, 0.5), (2.0, 0.0)
            ],
        ]
        regions = xor_geometry(alg, elbow, triangle, rect)
        @test are_regions_equal(regions, expected)

        regions = xor_geometry(alg, elbow, triangle)
        regions = xor_geometry(alg, regions, rect)
        @test are_regions_equal(regions, expected)
    end
end

@testset "expanding frontier" begin
    # Explicitly tests the expanding frontier in segments_to_paths
    # with face_selection=SPLIT_FACES
    poly1 = [(5.0, 1.0), (10.0, 6.0), (0.0, 6.0), (5.0, 1.0), (2.0, 5.0), (8.0, 5.0)]
    poly2 = [(5.0, 5.0), (6.5, 4.0), (3.5, 4.0)]
    poly3 = [(4.0, 4.0), (6.0, 4.0), (5.0, 3.0)]
    poly4 = [(5.0, 3.0), (5.5, 2.0), (4.5, 2.0)]
    regions = union_geometry(poly1, poly2, poly3, poly4,
        face_selection=PolygonAlgorithms.SPLIT_FACES);
    expected = [
        [(0.0, 6.0), (5.0, 1.0), (2.0, 5.0), (5.0, 5.0), (8.0, 5.0), (5.0, 1.0), (10.0, 6.0)],
        [(6.5, 4.0), (5.0, 5.0), (3.5, 4.0), (4.0, 4.0), (5.0, 3.0), (6.0, 4.0)],
        [(4.5, 2.0), (5.5, 2.0), (5.0, 3.0)],
    ]
    @test are_regions_equal(regions, expected)
    regions = union_geometry(poly1, poly2, poly3, poly4,
        face_selection=PolygonAlgorithms.MERGE_FACES);
    expected = [
        [(5.0, 1.0), (10.0, 6.0), (0.0, 6.0)],
        # hole:
        [
            (5.0, 1.0), (2.0, 5.0), (5.0, 5.0), (3.5, 4.0), (4.0, 4.0), 
            (5.0, 3.0), (4.5, 2.0), (5.5, 2.0), (5.0, 3.0), (6.0, 4.0), 
            (6.5, 4.0), (5.0, 5.0), (8.0, 5.0)
        ]
    ]
    @test are_regions_equal(regions, expected)
end

end