using PolygonAlgorithms
using PolygonAlgorithms: MartinezRuedaAlg, Polygon, PointSet, translate
using PolygonAlgorithms: are_equivalent_polygons

@testset "polygon boolean holes - $alg" for alg in [
    MartinezRuedaAlg(),
]
    
@testset "one inside the other" begin
    poly1 = Polygon([
        (3.0, 0.0), (3.0, 3.0), (0.0, 3.0), (0.0, 0.0)
    ])
    poly2 = Polygon([
        (2.0, 1.0), (2.0, 2.0), (1.0, 2.0), (1.0, 1.0)
    ])
    # difference
    regions = difference_geometry(alg, poly1, poly2)
    expected = [Polygon(poly1.exterior, holes=[reverse(poly2.exterior)])]
    @test are_equivalent_polygons(regions, expected)
    regions = difference_geometry(alg, poly2, poly1)
    @test isempty(regions)
    # intersection
    regions = intersect_geometry(alg, poly1, poly2)
    expected = [poly2]
    @test are_equivalent_polygons(regions, expected)
    # union
    regions = union_geometry(alg, poly1, poly2)
    expected = [poly1]
    @test are_equivalent_polygons(regions, expected)
    # in this case the polygon is treated as a hole
    regions = union_geometry(alg, [poly1, poly2])
    expected = [Polygon(poly1.exterior, holes=[reverse(poly2.exterior)])]
    @test are_equivalent_polygons(regions, expected)
    # XOR
    regions = xor_geometry(alg, poly1, poly2)
    expected = [Polygon(poly1.exterior; holes=[reverse(poly2.exterior)])]
    @test are_equivalent_polygons(regions, expected)
end

@testset "overlap" begin
    poly1 = Polygon([(3.0, 2.0), (3.0, 4.0), (1.0, 4.0), (1.0, 2.0)]);
    poly2 = Polygon([(3.0, 1.0),(2.0, 3.0), (1.0, 1.0)]);
    # union two different polygons
    regions = union_geometry(alg, poly1, poly2)
    expected = [Polygon([
        (1.0, 4.0), (1.0, 2.0), (1.5, 2.0), (1.0, 1.0),
        (3.0, 1.0), (2.5, 2.0), (3.0, 2.0), (3.0, 4.0)
    ])]
    @test are_equivalent_polygons(regions, expected)
    # both are part of the subject
    regions = union_geometry(alg, [poly1, poly2])
    expected = [
        Polygon([(1.0, 1.0), (3.0, 1.0), (2.5, 2.0), (1.5, 2.0)]),
        Polygon([(2.0, 3.0), (2.5, 2.0), (3.0, 2.0), (3.0, 4.0), (1.0, 4.0), (1.0, 2.0), (1.5, 2.0)]),
    ]
    @test are_equivalent_polygons(regions, expected)
end

@testset "rectangles" begin
    ## Same
    poly1 = Polygon([
        (3.0, 0.0), (3.0, 3.0), (0.0, 3.0), (0.0, 0.0)
    ])
    poly2 = poly1
    # difference
    regions = difference_geometry(alg, poly1, poly2)
    @test isempty(regions)
    # intersection
    regions = intersect_geometry(alg, poly1, poly2)
    expected = [poly1]
    @test are_equivalent_polygons(regions, expected)
    # union
    regions = union_geometry(alg, poly1, poly2)
    expected = [poly1]
    @test are_equivalent_polygons(regions, expected)
    # XOR
    regions = xor_geometry(alg, poly1, poly2)
    @test isempty(regions)

    ## Overlap
    poly2 = translate(poly1, (1.0, 1.0))
    # difference
    regions = difference_geometry(alg, poly1, poly2)
    expected = [
        Polygon([(1.0, 1.0), (1.0, 3.0), (0.0, 3.0), (0.0, 0.0), (3.0, 0.0), (3.0, 1.0)])
    ]
    @test are_equivalent_polygons(regions, expected)
    regions = difference_geometry(alg, poly2, poly1)
    expected = [
        Polygon([(1.0, 4.0), (1.0, 3.0), (3.0, 3.0),(3.0, 1.0), (4.0, 1.0), (4.0, 4.0)])
    ]
    @test are_equivalent_polygons(regions, expected)
    # intersection
    regions = intersect_geometry(alg, poly1, poly2)
    expected = [Polygon([(3.0, 3.0), (1.0, 3.0), (1.0, 1.0), (3.0, 1.0)])]
    @test are_equivalent_polygons(regions, expected)
    ## Union
    regions = union_geometry(alg, poly1, poly2)
    expected = [
        Polygon([(1.0, 4.0), (1.0, 3.0), (0.0, 3.0), (0.0, 0.0), (3.0, 0.0), (3.0, 1.0), (4.0, 1.0), (4.0, 4.0)])
    ]
    @test are_equivalent_polygons(regions, expected)
    ## XOR
    regions = xor_geometry(alg, poly1, poly2)
    expected = [
        Polygon([(0.0, 0.0), (3.0, 0.0), (3.0, 1.0), (1.0, 1.0), (1.0, 3.0), (0.0, 3.0)]),
        Polygon([(4.0, 1.0), (4.0, 4.0), (1.0, 4.0), (1.0, 3.0), (3.0, 3.0), (3.0, 1.0)]),
    ]
    @test are_equivalent_polygons(regions, expected)

    ## No overlap
    poly2 = translate(poly1, (4.0, 3.0))
    # difference
    regions = difference_geometry(alg, poly1, poly2)
    expected = [poly1]
    @test are_equivalent_polygons(regions, expected)
    regions = difference_geometry(alg, poly2, poly1)
    expected = [poly2]
    @test are_equivalent_polygons(regions, expected)
    # intersection
    regions = intersect_geometry(alg, poly1, poly2)
    @test isempty(regions)
    ## Union
    regions = union_geometry(alg, poly1, poly2)
    expected = [poly1, poly2]
    @test are_equivalent_polygons(regions, expected)
    ## XOR
    regions = xor_geometry(alg, poly1, poly2)
    expected = [poly1, poly2]
    @test are_equivalent_polygons(regions, expected)
end

@testset "triangle - holes" begin
    @testset "hole in centre" begin
        poly1 = Polygon([
            (9.0, 1.0), (5.0, 6.0), (1.0, 1.0)
        ])
        poly2 = Polygon([
            (4.0, 3.5), (6.0, 3.5), (5.0, 2.0)
        ])
        regions = difference_geometry(alg, poly1, poly2)
        expected = [Polygon(
            [(9.0, 1.0), (5.0, 6.0), (1.0, 1.0)],
            holes=[[(4.0, 3.5), (6.0, 3.5), (5.0, 2.0)]]
        )]
        @test are_equivalent_polygons(regions, expected)
    end

    @testset "hole touch centre" begin
        poly1 = Polygon([
            (9.0, 1.0), (5.0, 6.0), (1.0, 1.0)
        ])
        poly2 = Polygon([
            (3.0, 3.5), (6.0, 3.5), (5.0, 2.0)
        ])
        regions = difference_geometry(alg, poly1, poly2; face_selection=MERGE_FACES)
        expected = [Polygon(
            [(3.0, 3.5), (1.0, 1.0), (9.0, 1.0), (5.0, 6.0)],
            holes=[[(3.0, 3.5), (6.0, 3.5), (5.0, 2.0)]]
        )]
        @test are_equivalent_polygons(regions, expected)
        regions = difference_geometry(alg, poly1, poly2; face_selection=SPLIT_FACES)
        expected = [Polygon(
            [(3.0, 3.5), (6.0, 3.5), (5.0, 2.0), (3.0, 3.5), (1.0, 1.0), (9.0, 1.0), (5.0, 6.0)]
        )]
        @test are_equivalent_polygons(regions, expected)
    end

    @testset "hole only touching edges" begin
        poly1 = Polygon([
            (9.0, 1.0), (5.0, 6.0), (1.0, 1.0)
        ])
        poly2 = Polygon([
            (3.0, 3.5), (7.0, 3.5), (5.0, 1.0)
        ])
        regions = difference_geometry(alg, poly1, poly2; face_selection=MERGE_FACES)
        expected = [Polygon(
            [(3.0, 3.5), (1.0, 1.0), (5.0, 1.0), (9.0, 1.0), (7.0, 3.5), (5.0, 6.0)],
            holes=[[(3.0, 3.5), (7.0, 3.5), (5.0, 1.0)]]
        )]
        @test are_equivalent_polygons(regions, expected)
        regions = difference_geometry(alg, poly1, poly2; face_selection=SPLIT_FACES)
        expected = [
            Polygon([(5.0, 1.0), (3.0, 3.5), (1.0, 1.0)]),
            Polygon([(5.0, 1.0), (9.0, 1.0), (7.0, 3.5)]),
            Polygon([(7.0, 3.5), (5.0, 6.0), (3.0, 3.5)]),
        ]
        @test are_equivalent_polygons(regions, expected)
    end
end

@testset "concave <>" begin 
    poly1 = Polygon([
        (6.0, 6.0), (0.0, 12.0), (-1.0, 11.0), (4.0, 6.0), (-1.0, 1.0), (0.0, 0.0),  
    ])
    poly2 = Polygon([
        (-4.0, 6.0), (2.0, 12.0), (3.0, 11.0), (-2.0, 6.0), (3.0, 1.0), (2.0, 0.0),  
    ])
    # difference
    regions = difference_geometry(alg, poly1, poly2)
    expected = [
        Polygon([(0.0, 0.0), (1.0, 1.0), (0.0, 2.0), (-1.0, 1.0)]),
        Polygon([(2.0, 10.0), (1.0, 9.0), (4.0, 6.0), (1.0, 3.0), (2.0, 2.0), (6.0, 6.0)]),
        Polygon([(-1.0, 11.0), (0.0, 10.0), (1.0, 11.0), (0.0, 12.0)]),
    ]
    @test are_equivalent_polygons(regions, expected)
    regions = difference_geometry(alg, poly2, poly1)
    expected = [
        Polygon([(3.0, 11.0), (2.0, 12.0), (1.0, 11.0), (2.0, 10.0)]),
        Polygon([(3.0, 1.0), (2.0, 2.0), (1.0, 1.0), (2.0, 0.0)]),
        Polygon([(-2.0, 6.0), (1.0, 9.0), (0.0, 10.0), (-4.0, 6.0), (0.0, 2.0), (1.0, 3.0)]),
    ]
    @test are_equivalent_polygons(regions, expected)
    # intersection
    regions = intersect_geometry(alg, poly1, poly2)
    expected = [
        Polygon([(1.0, 11.0), (0.0, 10.0), (1.0, 9.0), (2.0, 10.0)]),
        Polygon([(1.0, 1.0), (2.0, 2.0), (1.0, 3.0), (0.0, 2.0)]),
    ]
    @test are_equivalent_polygons(regions, expected)
    ## Union
    regions = union_geometry(alg, poly1, poly2)
    expected = [
        Polygon(
            [(2.0, 10.0), (3.0, 11.0), (2.0, 12.0), (1.0, 11.0), (0.0, 12.0), (-1.0, 11.0), (0.0, 10.0), (-4.0, 6.0), (0.0, 2.0), (-1.0, 1.0), (0.0, 0.0), (1.0, 1.0), (2.0, 0.0), (3.0, 1.0), (2.0, 2.0), (6.0, 6.0)];
            holes=[[(4.0, 6.0), (1.0, 3.0), (-2.0, 6.0), (1.0, 9.0)]],
        )
    ]
    @test are_equivalent_polygons(regions, expected)
    ## XOR
    regions = xor_geometry(alg, poly1, poly2)
    expected = [
        Polygon([(1.0, 11.0), (2.0, 10.0), (3.0, 11.0), (2.0, 12.0)]),
        Polygon([(2.0, 2.0), (1.0, 1.0), (2.0, 0.0), (3.0, 1.0)]),
        Polygon([(-1.0, 11.0), (0.0, 10.0), (1.0, 11.0), (0.0, 12.0)]),
        Polygon([(1.0, 3.0), (-2.0, 6.0), (1.0, 9.0), (0.0, 10.0), (-4.0, 6.0), (0.0, 2.0)]),
        Polygon([(-1.0, 1.0), (0.0, 0.0), (1.0, 1.0), (0.0, 2.0)]),
        Polygon([(2.0, 10.0), (1.0, 9.0), (4.0, 6.0), (1.0, 3.0), (2.0, 2.0), (6.0, 6.0)]),
    ]
    @test are_equivalent_polygons(regions, expected)
end

@testset "self-intersect & rectangle" begin
    self_intersect = Polygon([
        (11.0, 0.0), (11.0, 2.0), (6.0, -2.0), (2.0, 2.0), (0.0, 0.0)
    ])
    rectangle_horiz = Polygon([
        (12.0, 0.0), (12.0, 3.0), (-1.0, 3.0), (-1.0, 0.0)
    ]);
    poly1 = self_intersect
    poly2 = rectangle_horiz
    # difference
    regions = difference_geometry(alg, poly1, poly2)
    expected = [Polygon([(4.0, 0.0), (6.0, -2.0), (8.5, 0.0)])]
    @test are_equivalent_polygons(regions, expected)
    regions = difference_geometry(alg, poly2, poly1)
    expected = [
        Polygon([
            (-1.0, 3.0), (-1.0, 0.0), (0.0, 0.0), (2.0, 2.0), (4.0, 0.0), (8.5, 0.0),
            (11.0, 2.0), (11.0, 0.0), (12.0, 0.0), (12.0, 3.0)
        ])
    ]
    @test are_equivalent_polygons(regions, expected)
    ## Intersection
    regions = intersect_geometry(alg, poly1, poly2)
    expected = [
        Polygon([(8.5, 0.0), (11.0, 0.0), (11.0, 2.0)]),
        Polygon([(4.0, 0.0), (2.0, 2.0), (0.0, 0.0)]),
    ]
    @test are_equivalent_polygons(regions, expected)
    ## Union
    regions = union_geometry(alg, poly1, poly2)
    expected = [
        Polygon([(12.0, 3.0), (-1.0, 3.0), (-1.0, 0.0), (0.0, 0.0), (4.0, 0.0), (6.0, -2.0), (8.5, 0.0), (11.0, 0.0), (12.0, 0.0)])
    ]
    @test are_equivalent_polygons(regions, expected)
    ## XOR
    regions = xor_geometry(alg, poly1, poly2)
    expected = [
        Polygon([(12.0, 3.0), (-1.0, 3.0), (-1.0, 0.0), (0.0, 0.0), (2.0, 2.0), (4.0, 0.0), (6.0, -2.0), (8.5, 0.0), (11.0, 2.0), (11.0, 0.0), (12.0, 0.0)])
    ]
    @test are_equivalent_polygons(regions, expected)
end

@testset "Ring & square" begin
    poly1 = Polygon(
        [(5.0, 5.0), (0.0, 10.0), (-5.0, 5.0), (0.0, 0.0)];
        holes=[[(0.0, 1.0), (-4.0, 5.0), (0.0, 9.0), (4.0, 5.0)]]
        )
    poly2 = Polygon([(8.0, 6.0), (3.0, 11.0), (-2.0, 6.0), (3.0, 1.0)])
    # difference
    regions = difference_geometry(alg, poly1, poly2)
    expected = [
        Polygon([
            (1.5, 2.5), (0.0, 1.0), (-4.0, 5.0), (0.0, 9.0), (0.5, 8.5), (1.0, 9.0),
            (0.0, 10.0), (-5.0, 5.0), (0.0, 0.0), (2.0, 2.0)
        ]),
    ]
    @test are_equivalent_polygons(regions, expected)
    regions = difference_geometry(alg, poly2, poly1)
    expected = [
        Polygon([(0.5, 8.5), (-2.0, 6.0), (1.5, 2.5), (4.0, 5.0)]),
        Polygon([(5.0, 5.0), (2.0, 2.0), (3.0, 1.0), (8.0, 6.0), (3.0, 11.0), (1.0, 9.0)]),
    ]
    @test are_equivalent_polygons(regions, expected)
    ## Intersection
    regions = intersect_geometry(alg, poly1, poly2)
    expected = [
        Polygon([(1.0, 9.0), (0.5, 8.5), (4.0, 5.0), (1.5, 2.5), (2.0, 2.0), (5.0, 5.0)])
    ]
    @test are_equivalent_polygons(regions, expected)
    ## Union
    regions = union_geometry(alg, poly1, poly2)
    expected = [
        Polygon(
            [(3.0, 11.0), (1.0, 9.0), (0.0, 10.0), (-5.0, 5.0), (0.0, 0.0), (2.0, 2.0), (3.0, 1.0), (8.0, 6.0)];
            holes=[[(1.5, 2.5), (0.0, 1.0), (-4.0, 5.0), (0.0, 9.0), (0.5, 8.5), (-2.0, 6.0)]]
        )
    ]
    @test are_equivalent_polygons(regions, expected)
    ## XOR
    regions = xor_geometry(alg, poly1, poly2)
    expected = [
        Polygon([(-4.0, 5.0), (0.0, 9.0), (0.5, 8.5), (1.0, 9.0), (0.0, 10.0), (-5.0, 5.0), (0.0, 0.0), (2.0, 2.0), (1.5, 2.5), (0.0, 1.0)]),
        Polygon([(5.0, 5.0), (2.0, 2.0), (3.0, 1.0), (8.0, 6.0), (3.0, 11.0), (1.0, 9.0)]),
        Polygon([(0.5, 8.5), (-2.0, 6.0), (1.5, 2.5), (4.0, 5.0)]),
    ]
    @test are_equivalent_polygons(regions, expected)
end

@testset "Overlapping rings" begin
    poly1 = Polygon(
        [(0.0, 0.0), (-5.0, 5.0), (0.0, 10.0), (5.0, 5.0)],
        [[(0.0, 1.0), (-4.0, 5.0), (0.0, 9.0), (4.0, 5.0)]]
        )
    poly2 = translate(poly1, (3.0, 1.0))
    # difference
    regions = difference_geometry(alg, poly1, poly2)
    expected = [
        Polygon([(2.5, 2.5), (5.0, 5.0), (1.5, 8.5), (1.0, 8.0), (4.0, 5.0), (2.0, 3.0)]),
        Polygon([(1.5, 2.5), (0.0, 1.0), (-4.0, 5.0), (0.0, 9.0), (0.5, 8.5), (1.0, 9.0), (0.0, 10.0), (-5.0, 5.0), (0.0, 0.0), (2.0, 2.0)]),
    ]
    @test are_equivalent_polygons(regions, expected)
    regions = difference_geometry(alg, poly2, poly1)
    expected = [
        Polygon([(2.0, 3.0), (-1.0, 6.0), (1.0, 8.0), (0.5, 8.5), (-2.0, 6.0), (1.5, 2.5)]),
        Polygon([(3.0, 2.0), (2.5, 2.5), (2.0, 2.0), (3.0, 1.0), (8.0, 6.0), (3.0, 11.0), (1.0, 9.0), (1.5, 8.5), (3.0, 10.0), (7.0, 6.0)]),
    ]
    @test are_equivalent_polygons(regions, expected)
    ## Intersection
    regions = intersect_geometry(alg, poly1, poly2)
    expected = [
        Polygon([(1.0, 8.0), (1.5, 8.5), (1.0, 9.0), (0.5, 8.5)]),
        Polygon([(2.0, 3.0), (1.5, 2.5), (2.0, 2.0), (2.5, 2.5)]),
    ]
    @test are_equivalent_polygons(regions, expected)
    ## Union
    regions = union_geometry(alg, poly1, poly2)
    expected = [
        Polygon(
            [(-5.0, 5.0), (0.0, 0.0), (2.0, 2.0), (3.0, 1.0), (8.0, 6.0), (3.0, 11.0), (1.0, 9.0), (0.0, 10.0)];
            holes=[
                [(3.0, 10.0), (7.0, 6.0), (3.0, 2.0), (2.5, 2.5), (5.0, 5.0), (1.5, 8.5)],
                [(0.5, 8.5), (-2.0, 6.0), (1.5, 2.5), (0.0, 1.0), (-4.0, 5.0), (0.0, 9.0)],
                [(2.0, 3.0), (-1.0, 6.0), (1.0, 8.0), (4.0, 5.0)]
            ]
        )
    ]
    @test are_equivalent_polygons(regions, expected)
    ## XOR
    regions = xor_geometry(alg, poly1, poly2)
    expected = [
        Polygon([(1.0, 9.0), (1.5, 8.5), (3.0, 10.0), (7.0, 6.0), (3.0, 2.0), (2.5, 2.5), (2.0, 2.0), (3.0, 1.0), (8.0, 6.0), (3.0, 11.0)]),
        Polygon([(0.0, 9.0), (0.5, 8.5), (1.0, 9.0), (0.0, 10.0), (-5.0, 5.0), (0.0, 0.0), (2.0, 2.0), (1.5, 2.5), (0.0, 1.0), (-4.0, 5.0)]),
        Polygon([(1.0, 8.0), (4.0, 5.0), (2.0, 3.0), (2.5, 2.5), (5.0, 5.0), (1.5, 8.5)]),
        Polygon([(2.0, 3.0), (-1.0, 6.0), (1.0, 8.0), (0.5, 8.5), (-2.0, 6.0), (1.5, 2.5)]),
    ]
    @test are_equivalent_polygons(regions, expected)
end

@testset verbose=true "multi - appearing holes" begin
    elbow = Polygon([
        (0.0, -1.0), (0.0, 2.0), (3.0, 2.0), (3.0, 1.0), (1.0, 1.0), (1.0, -1.0)
    ])
    triangle = Polygon([
        (0.8, -0.8), (2.5, 1.8), (2.5, -0.8)
    ])
    rect = Polygon([(0.5, 0.0), (0.5, 0.5), (2.0, 0.5), (2.0, 0.0)])

    ## Intersection
    regions = intersect_geometry(alg, elbow, triangle, rect)
    @test isempty(regions)
    # Union
    expected = [
        Polygon(
            [(0.0, 2.0), (0.0, -1.0), (1.0, -1.0), (1.0, -0.8), (2.5, -0.8), (2.5, 1.0), (3.0, 1.0), (3.0, 2.0)];
            holes=[
                [(1.0, 0.5), (1.0, 1.0), (1.976923, 1.0), (1.65, 0.5)],
                [(1.0, -0.494118), (1.0, 0.0), (1.323077, 0.0)],
            ]
        )
    ]
    regions = union_geometry(alg, elbow, triangle, rect)
    @test are_equivalent_polygons(regions, expected)
    ## Difference
    expected = [
        Polygon([
            (0.5, 0.0), (0.5, 0.5), (1.0, 0.5), (1.0, 1.0), (1.976923, 1.0),
            (2.5, 1.8), (2.5, 1.0), (3.0, 1.0), (3.0, 2.0), (0.0, 2.0), (0.0, -1.0),
            (1.0, -1.0), (1.0, -0.8), (0.8, -0.8), (1.0, -0.494118), (1.0, 0.0)
        ]),
        Polygon([
            (1.323077, 0.0), (1.0, -0.494118), (1.0, -0.8), (2.5, -0.8),
            (2.5, 1.0), (1.976923, 1.0), (1.65, 0.5), (2.0, 0.5), (2.0, 0.0)
        ]),
    ]
    regions = difference_geometry(alg, [elbow, triangle], rect)
    @test are_equivalent_polygons(regions, expected)

    # XOR
    expected = [
        Polygon([(0.5, 0.0), (0.5, 0.5), (1.0, 0.5), (1.0, 1.0), (1.976923, 1.0), (2.5, 1.8), (2.5, 1.0), (3.0, 1.0), (3.0, 2.0), (0.0, 2.0), (0.0, -1.0), (1.0, -1.0), (1.0, -0.8), (0.8, -0.8), (1.0, -0.494118), (1.0, 0.0)]),
        Polygon([(1.0, 0.0), (1.323077, 0.0), (1.65, 0.5), (1.0, 0.5)]),
        Polygon([(1.323077, 0.0), (1.0, -0.494118), (1.0, -0.8), (2.5, -0.8), (2.5, 1.0), (1.976923, 1.0), (1.65, 0.5), (2.0, 0.5), (2.0, 0.0)]),
    ]
    regions = xor_geometry(alg, elbow, triangle, rect)
    @test are_equivalent_polygons(regions, expected)
end

end