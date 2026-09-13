using PolygonAlgorithms
using PolygonAlgorithms: MartinezRuedaAlg, Polygon, translate
using PolygonAlgorithms: are_equivalent_polygons

@testset "polygon boolean pairs - $alg" for alg in [
    MartinezRuedaAlg(),
]

# See test/intersect_convex.jl and test/intersect_concave.jl for intersection tests.
    
@testset "one inside the other" begin
    poly1 = [
        (3.0, 0.0), (3.0, 3.0), (0.0, 3.0), (0.0, 0.0)
    ]
    poly2 = [
        (2.0, 1.0), (2.0, 2.0), (1.0, 2.0), (1.0, 1.0)
    ]
    # difference
    regions = difference_geometry(alg, poly1, poly2)
    expected = [poly1, reverse(poly2)] # second is a hole
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
    regions = difference_geometry(alg, poly2, poly1)
    @test isempty(regions)
    # union
    regions = union_geometry(alg, poly1, poly2)
    expected = [poly1]
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
    # XOR
    regions = xor_geometry(alg, poly1, poly2)
    expected = [reverse(poly2), poly1] # first is a hole
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
end

@testset "rectangles" begin
    ## Same
    poly1 = [
        (3.0, 0.0), (3.0, 3.0), (0.0, 3.0), (0.0, 0.0)
    ]
    poly2 = poly1
    # difference
    regions = difference_geometry(alg, poly1, poly2)
    @test isempty(regions)
    # union
    regions = union_geometry(alg, poly1, poly2)
    expected = [poly1]
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
    # XOR
    regions = xor_geometry(alg, poly1, poly2)
    @test isempty(regions)

    ## Overlap
    poly2 = translate(poly1, (1.0, 1.0))
    # difference
    regions = difference_geometry(alg, poly1, poly2)
    expected = [[(1.0, 1.0), (1.0, 3.0), (0.0, 3.0), (0.0, 0.0), (3.0, 0.0), (3.0, 1.0)]]
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
    regions = difference_geometry(alg, poly2, poly1)
    expected = [[(1.0, 4.0), (1.0, 3.0), (3.0, 3.0), (3.0, 1.0), (4.0, 1.0), (4.0, 4.0)]]
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
    ## Union
    regions = union_geometry(alg, poly1, poly2)
    expected = [[(1.0, 4.0), (1.0, 3.0), (0.0, 3.0), (0.0, 0.0), (3.0, 0.0), (3.0, 1.0), (4.0, 1.0), (4.0, 4.0)]]
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
    ## XOR
    regions = xor_geometry(alg, poly1, poly2; face_selection=MERGE_FACES)
    expected = [
        # exterior
        [(4.0, 4.0), (1.0, 4.0), (1.0, 3.0), (0.0, 3.0), (0.0, 0.0), (3.0, 0.0), (3.0, 1.0), (4.0, 1.0)],
        # holes
        [(3.0, 1.0), (1.0, 1.0), (1.0, 3.0), (3.0, 3.0)]
    ]
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
    regions = xor_geometry(alg, poly1, poly2)
    expected = [
        [(0.0, 0.0), (3.0, 0.0), (3.0, 1.0), (1.0, 1.0), (1.0, 3.0), (0.0, 3.0)],
        [(4.0, 1.0), (4.0, 4.0), (1.0, 4.0), (1.0, 3.0), (3.0, 3.0), (3.0, 1.0)],
    ]
    @test are_equivalent_polygons(regions, expected; match_reverse=false)

    ## No overlap
    poly2 = translate(poly1, (4.0, 3.0))
    # difference
    regions = difference_geometry(alg, poly1, poly2)
    expected = [poly1]
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
    regions = difference_geometry(alg, poly2, poly1)
    expected = [poly2]
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
    ## Union
    regions = union_geometry(alg, poly1, poly2)
    expected = [poly1, poly2]
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
    ## XOR
    regions = xor_geometry(alg, poly1, poly2)
    expected = [poly1, poly2]
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
end

@testset "cross" begin
    poly1 = [
        (1.0, 0.0), (1.0, 2.0), (0.0, 2.0), (0.0, 0.0)
    ]
    poly2 = [
        (-1.0, 0.5), (2.0, 0.5), (2.0, 1.5), (-1.0, 1.5)
    ]
    # difference
    regions = difference_geometry(alg, poly1, poly2)
    expected = [
        [(0.0, 0.5), (0.0, 0.0), (1.0, 0.0), (1.0, 0.5)],
        [(0.0, 2.0), (-0.0, 1.5), (1.0, 1.5), (1.0, 2.0)],
    ]
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
    regions = difference_geometry(alg, poly2, poly1)
    expected = [
        [(-1.0, 1.5), (-1.0, 0.5), (-0.0, 0.5), (0.0, 1.5)],
        [(1.0, 1.5), (1.0, 0.5), (2.0, 0.5), (2.0, 1.5)],
    ]
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
    ## Union
    regions = union_geometry(alg, poly1, poly2)
    expected = [
        [
            (1.0, 1.5), (1.0, 2.0), (0.0, 2.0), (-0.0, 1.5),
            (-1.0, 1.5), (-1.0, 0.5), (-0.0, 0.5), (0.0, 0.0),
            (1.0, 0.0), (1.0, 0.5), (2.0, 0.5), (2.0, 1.5)
        ]
    ]
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
    ## XOR
    regions = xor_geometry(alg, poly1, poly2)
    expected = [
        [(1.0, 0.5), (0.0, 0.5), (0.0, 0.0), (1.0, 0.0)],
        [(-1.0, 1.5), (-1.0, 0.5), (0.0, 0.5), (0.0, 1.5)],
        [(1.0, 1.5), (1.0, 0.5), (2.0, 0.5), (2.0, 1.5)],
        [(1.0, 1.5), (1.0, 2.0), (0.0, 2.0), (0.0, 1.5)],
    ]
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
end

@testset "concave <>" begin 
    poly1 = [
        (6.0, 6.0), (0.0, 12.0), (-1.0, 11.0), (4.0, 6.0), (-1.0, 1.0), (0.0, 0.0),  
    ]
    poly2 = [
        (2.0, 0.0), (3.0, 1.0),  (-2.0, 6.0), (3.0, 11.0), (2.0, 12.0), (-4.0, 6.0),  
    ]
    # difference
    regions = difference_geometry(alg, poly1, poly2)
    expected = [
        [(0.0, 2.0), (-1.0, 1.0), (0.0, 0.0), (1.0, 1.0)],
        [(0.0, 12.0), (-1.0, 11.0), (0.0, 10.0), (1.0, 11.0)],
        [(2.0, 10.0), (1.0, 9.0), (4.0, 6.0), (1.0, 3.0), (2.0, 2.0), (6.0, 6.0)],
    ]
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
    regions = difference_geometry(alg, poly2, poly1)
    expected = [
        [(0.0, 10.0), (-4.0, 6.0), (-0.0, 2.0), (1.0, 3.0), (-2.0, 6.0), (1.0, 9.0)],
        [(2.0, 2.0), (1.0, 1.0), (2.0, 0.0), (3.0, 1.0)],
        [(2.0, 12.0), (1.0, 11.0), (2.0, 10.0), (3.0, 11.0)],
    ]
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
    ## Union
    regions = union_geometry(alg, poly1, poly2)
    expected = [
        [(4.0, 6.0), (1.0, 3.0), (-2.0, 6.0), (1.0, 9.0)], # hole
        [
            (2.0, 10.0), (3.0, 11.0), (2.0, 12.0), (1.0, 11.0),
            (0.0, 12.0), (-1.0, 11.0), (0.0, 10.0), (-4.0, 6.0), (-0.0, 2.0),
            (-1.0, 1.0), (0.0, 0.0), (1.0, 1.0), (2.0, 0.0), (3.0, 1.0), (2.0, 2.0), (6.0, 6.0)]
    ]
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
    ## XOR
    regions = xor_geometry(alg, poly1, poly2)
    expected = [
        [(1.0, 11.0), (2.0, 10.0), (3.0, 11.0), (2.0, 12.0)],
        [(2.0, 2.0), (1.0, 1.0), (2.0, 0.0), (3.0, 1.0)],
        [(-1.0, 11.0), (0.0, 10.0), (1.0, 11.0), (0.0, 12.0)],
        [(1.0, 3.0), (-2.0, 6.0), (1.0, 9.0), (0.0, 10.0), (-4.0, 6.0), (0.0, 2.0)],
        [(-1.0, 1.0), (0.0, 0.0), (1.0, 1.0), (0.0, 2.0)],
        [(2.0, 10.0), (1.0, 9.0), (4.0, 6.0), (1.0, 3.0), (2.0, 2.0), (6.0, 6.0)],
    ]
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
end

@testset "self-intersect & rectangle" begin
    self_intersect = [
        (0.0, 0.0), (2.0, 2.0), (6.0, -2.0), (11.0, 2.0), (11.0, 0.0)
    ]
    rectangle_horiz = [
        (-1.0, 0.0), (-1.0, 3.0), (12.0, 3.0), (12.0, 0.0)
    ];
    poly1 = self_intersect
    poly2 = rectangle_horiz
    # difference
    regions = difference_geometry(alg, poly1, poly2)
    expected = [[(4.0, 0.0), (6.0, -2.0), (8.5, -0.0)]]
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
    regions = difference_geometry(alg, poly2, poly1)
    expected = [
        [
            (-1.0, 3.0), (-1.0, 0.0), (0.0, 0.0), (2.0, 2.0), (4.0, 0.0),
            (8.5, -0.0), (11.0, 2.0), (11.0, 0.0), (12.0, 0.0), (12.0, 3.0)
        ]
    ]
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
    ## Union
    regions = union_geometry(alg, poly1, poly2)
    expected = [
         [(12.0, 3.0), (-1.0, 3.0), (-1.0, 0.0), (0.0, 0.0), (4.0, 0.0), (6.0, -2.0), (8.5, 0.0), (11.0, 0.0), (12.0, 0.0)]
    ]
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
    ## XOR
    regions = xor_geometry(alg, poly1, poly2)
    expected = [
        [(12.0, 3.0), (-1.0, 3.0), (-1.0, 0.0), (0.0, 0.0), (2.0, 2.0), (4.0, 0.0), (6.0, -2.0), (8.5, 0.0), (11.0, 2.0), (11.0, 0.0), (12.0, 0.0)],
    ]
    @test are_equivalent_polygons(regions, expected; match_reverse=false)
end

end