using PolygonAlgorithms
using PolygonAlgorithms: Polygon, translate

@testset "polygon boolean holes - $alg" for alg in [
    MartinezRuedaAlg(),
]
    
@testset "one inside the other" begin
    poly1 = Polygon([
        (0.0, 0.0), (0.0, 3.0), (3.0, 3.0), (3.0, 0.0)
    ])
    poly2 = Polygon([
        (1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0)
    ])
    # difference
    regions = difference_geometry(poly1, poly2, alg)
    expected = [poly1.exterior, poly2.exterior] # second is a hole
    @test are_regions_equal(regions, expected)
    regions = difference_geometry(poly2, poly1, alg)
    @test isempty(regions)
    # intersection
    regions = intersect_geometry(poly1, poly2, alg)
    expected = [poly2.exterior]
    @test are_regions_equal(regions, expected)
    # union
    regions = union_geometry(poly1, poly2, alg)
    expected = [poly1.exterior]
    @test are_regions_equal(regions, expected)
    # XOR
    regions = xor_geometry(poly1, poly2, alg)
    expected = [poly2.exterior, poly1.exterior] # first is a hole
    @test are_regions_equal(regions, expected)
end

@testset "rectangles" begin
    ## Same
    poly1 = Polygon([
        (0.0, 0.0), (0.0, 3.0), (3.0, 3.0), (3.0, 0.0)
    ])
    poly2 = poly1
    # difference
    regions = difference_geometry(poly1, poly2, alg)
    @test isempty(regions)
    # intersection
    regions = intersect_geometry(poly1, poly2, alg)
    expected = [poly1.exterior]
    @test are_regions_equal(regions, expected)
    # union
    regions = union_geometry(poly1, poly2, alg)
    expected = [poly1.exterior]
    @test are_regions_equal(regions, expected)
    # XOR
    regions = xor_geometry(poly1, poly2, alg)
    @test isempty(regions)

    ## Overlap
    poly2 = translate(poly1, (1.0, 1.0))
    # difference
    regions = difference_geometry(poly1, poly2, alg)
    expected = [[(3.0, 1.0), (3.0, 0.0), (0.0, 0.0), (0.0, 3.0), (1.0, 3.0), (1.0, 3.0), (1.0, 1.0)]]
    @test are_regions_equal(regions, expected)
    regions = difference_geometry(poly2, poly1, alg)
    expected = [[(4.0, 4.0), (4.0, 1.0), (3.0, 1.0), (3.0, 3.0), (3.0, 3.0), (1.0, 3.0), (1.0, 4.0)]]
    @test are_regions_equal(regions, expected)
    # intersection
    regions = intersect_geometry(poly1, poly2, alg)
    expected = [[(3.0, 3.0), (3.0, 1.0), (1.0, 1.0), (1.0, 3.0)]]
    @test are_regions_equal(regions, expected)
    ## Union
    regions = union_geometry(poly1, poly2, alg)
    expected = [[(4.0, 4.0), (4.0, 1.0), (3.0, 1.0), (3.0, 0.0), (0.0, 0.0), (0.0, 3.0), (1.0, 3.0), (1.0, 4.0)]]
    @test are_regions_equal(regions, expected)
    ## XOR
    regions = xor_geometry(poly1, poly2, alg)
    expected = [
        [(3.0, 1.0), (3.0, 0.0), (0.0, 0.0), (0.0, 3.0), (1.0, 3.0), (1.0, 3.0), (1.0, 1.0)],
        [(4.0, 4.0), (4.0, 1.0), (3.0, 1.0), (3.0, 3.0), (3.0, 3.0), (1.0, 3.0), (1.0, 4.0)],
    ]
    @test are_regions_equal(regions, expected)

    ## No overlap
    poly2 = translate(poly1, (4.0, 3.0))
    # difference
    regions = difference_geometry(poly1, poly2, alg)
    expected = [poly1.exterior]
    @test are_regions_equal(regions, expected)
    regions = difference_geometry(poly2, poly1, alg)
    expected = [poly2.exterior]
    @test are_regions_equal(regions, expected)
    # intersection
    regions = intersect_geometry(poly1, poly2, alg)
    @test isempty(regions)
    ## Union
    regions = union_geometry(poly1, poly2, alg)
    expected = [poly1.exterior, poly2.exterior]
    @test are_regions_equal(regions, expected)
    ## XOR
    regions = xor_geometry(poly1, poly2, alg)
    expected = [poly1.exterior, poly2.exterior]
    @test are_regions_equal(regions, expected)
end

@testset "concave <>" begin 
    poly1 = Polygon([
        (0.0, 0.0), (-1.0, 1.0),  (4.0, 6.0), (-1.0, 11.0), (0.0, 12.0), (6.0, 6.0),  
    ])
    poly2 = Polygon([
        (2.0, 0.0), (3.0, 1.0),  (-2.0, 6.0), (3.0, 11.0), (2.0, 12.0), (-4.0, 6.0),  
    ])
    # difference
    regions = difference_geometry(poly1, poly2, alg)
    expected = [
        [(-0.0, 2.0), (-1.0, 1.0), (0.0, 0.0), (1.0, 1.0)],
        [(0.0, 12.0), (-1.0, 11.0), (0.0, 10.0), (1.0, 11.0)],
        [(2.0, 10.0), (1.0, 9.0), (1.0, 9.0), (4.0, 6.0), (1.0, 3.0), (2.0, 2.0), (6.0, 6.0)],
    ]
    @test are_regions_equal(regions, expected)
    regions = difference_geometry(poly2, poly1, alg)
    expected = [
        [(0.0, 10.0), (-4.0, 6.0), (-0.0, 2.0), (1.0, 3.0), (-2.0, 6.0), (1.0, 9.0)],
        [(2.0, 2.0), (1.0, 1.0), (2.0, 0.0), (3.0, 1.0)],
        [(2.0, 12.0), (1.0, 11.0), (2.0, 10.0), (3.0, 11.0)],
    ]
    @test are_regions_equal(regions, expected)
    # intersection
    regions = intersect_geometry(poly1, poly2, alg)
    expected = [
        [(0.0, 10.0), (1.0, 11.0), (2.0, 10.0), (1.0, 9.0)],
        [(1.0, 3.0), (2.0, 2.0), (1.0, 1.0), (0.0, 2.0)],
    ]
    @test are_regions_equal(regions, expected)
    ## Union
    regions = union_geometry(poly1, poly2, alg)
    expected = [
        [(1.0, 9.0), (-2.0, 6.0), (1.0, 3.0), (4.0, 6.0)], # hole
        [(2.0, 10.0), (3.0, 11.0), (3.0, 11.0), (2.0, 12.0), (1.0, 11.0), (0.0, 12.0), (-1.0, 11.0), (0.0, 10.0), (-4.0, 6.0), (-0.0, 2.0), (-1.0, 1.0), (0.0, 0.0), (1.0, 1.0), (2.0, 0.0), (3.0, 1.0), (2.0, 2.0), (6.0, 6.0)]
    ]
    @test are_regions_equal(regions, expected)
    ## XOR
    regions = xor_geometry(poly1, poly2, alg)
    expected = [
        [(0.0, 12.0), (-1.0, 11.0), (0.0, 10.0), (-4.0, 6.0), (-0.0, 2.0), (-1.0, 1.0), (0.0, 0.0), (1.0, 1.0), (-0.0, 2.0), (1.0, 3.0), (-2.0, 6.0), (1.0, 9.0), (0.0, 10.0), (1.0, 11.0)],
        [(2.0, 10.0), (3.0, 11.0), (3.0, 11.0), (2.0, 12.0), (1.0, 11.0), (2.0, 10.0), (1.0, 9.0), (1.0, 9.0), (4.0, 6.0), (1.0, 3.0), (2.0, 2.0), (1.0, 1.0), (2.0, 0.0), (3.0, 1.0), (2.0, 2.0), (6.0, 6.0)]
    ]
    @test are_regions_equal(regions, expected)
end

@testset "self-intersect & rectangle" begin
    self_intersect = Polygon([
        (0.0, 0.0), (2.0, 2.0), (6.0, -2.0), (11.0, 2.0), (11.0, 0.0)
    ])
    rectangle_horiz = Polygon([
        (-1.0, 0.0), (-1.0, 3.0), (12.0, 3.0), (12.0, 0.0)
    ]);
    poly1 = self_intersect
    poly2 = rectangle_horiz
    # difference
    regions = difference_geometry(poly1, poly2, alg)
    expected = [[(4.0, 0.0), (6.0, -2.0), (8.5, -0.0)]]
    @test are_regions_equal(regions, expected)
    regions = difference_geometry(poly2, poly1, alg)
    expected = [
        [(12.0, 3.0), (12.0, 0.0), (11.0, 0.0), (11.0, 2.0), (11.0, 2.0), (8.5, -0.0), (4.0, 0.0), (2.0, 2.0), (0.0, 0.0), (-1.0, 0.0), (-1.0, 3.0)]
    ]
    @test are_regions_equal(regions, expected)
    ## Intersection
    regions = intersect_geometry(poly1, poly2, alg)
    expected = [
        [(4.0, 0.0), (0.0, 0.0), (2.0, 2.0)],
        [(8.5, -0.0), (11.0, 0.0), (11.0, 2.0)],
        [(4.0, 0.0), (8.5, -0.0)], # straight line
    ]
    @test are_regions_equal(regions, expected)
    ## Union
    regions = union_geometry(poly1, poly2, alg)
    expected = [
        [(12.0, 3.0), (12.0, 0.0), (8.5, 0.0), (6.0, -2.0), (4.0, 0.0), (-1.0, 0.0), (-1.0, 3.0)]
    ]
    @test are_regions_equal(regions, expected)
    ## XOR
    regions = xor_geometry(poly1, poly2, alg)
    expected = [
        [(12.0, 3.0), (12.0, 0.0), (11.0, 0.0), (11.0, 2.0), (11.0, 2.0), (6.0, -2.0), (2.0, 2.0), (0.0, 0.0), (-1.0, 0.0), (-1.0, 3.0)]
    ]
    @test are_regions_equal(regions, expected)
end

@testset "Ring & square" begin
    poly1 = Polygon(
        [(0.0, 0.0), (-5.0, 5.0), (0.0, 10.0), (5.0, 5.0)],
        [[(0.0, 1.0), (-4.0, 5.0), (0.0, 9.0), (4.0, 5.0)]]
        )
    poly2 = Polygon([(3.0, 1.0), (-2.0, 6.0), (3.0, 11.0), (8.0, 6.0)])
    # difference
    regions = difference_geometry(poly1, poly2, alg)
    expected = [
        [(1.5, 2.5), (0.0, 1.0), (-4.0, 5.0), (0.0, 9.0), (0.5, 8.5), (1.0, 9.0), (1.0, 9.0), (0.0, 10.0), (-5.0, 5.0), (0.0, 0.0), (2.0, 2.0)],
    ]
    @test are_regions_equal(regions, expected)
    regions = difference_geometry(poly2, poly1, alg)
    expected = [
        [(4.0, 5.0), (1.5, 2.5), (-2.0, 6.0), (0.5, 8.5)],
        [(3.0, 11.0), (1.0, 9.0), (1.0, 9.0), (5.0, 5.0), (2.0, 2.0), (3.0, 1.0), (8.0, 6.0)]
    ]
    @test are_regions_equal(regions, expected)
    ## Intersection
    regions = intersect_geometry(poly1, poly2, alg)
    expected = [
        [(5.0, 5.0), (2.0, 2.0), (1.5, 2.5), (4.0, 5.0), (4.0, 5.0), (0.5, 8.5), (1.0, 9.0)]
    ]
    @test are_regions_equal(regions, expected)
    ## Union
    regions = union_geometry(poly1, poly2, alg)
    expected = [
        [(-2.0, 6.0), (0.5, 8.5), (0.5, 8.5), (0.0, 9.0), (-4.0, 5.0), (0.0, 1.0), (1.5, 2.5)], # hole
        [(3.0, 11.0), (1.0, 9.0), (0.0, 10.0), (-5.0, 5.0), (0.0, 0.0), (2.0, 2.0), (3.0, 1.0), (8.0, 6.0)]
    ]
    @test are_regions_equal(regions, expected)
    ## XOR
    regions = xor_geometry(poly1, poly2, alg)
    expected = [
        [(-2.0, 6.0), (0.5, 8.5), (0.5, 8.5), (0.0, 9.0), (-4.0, 5.0), (0.0, 1.0), (1.5, 2.5)], # hole
        [(0.5, 8.5), (1.0, 9.0), (1.0, 9.0), (0.0, 10.0), (-5.0, 5.0), (0.0, 0.0), (2.0, 2.0), (1.5, 2.5), (4.0, 5.0)],
        [(3.0, 11.0), (1.0, 9.0), (1.0, 9.0), (5.0, 5.0), (2.0, 2.0), (3.0, 1.0), (8.0, 6.0)],
    ]
    @test are_regions_equal(regions, expected)
end

@testset "Overlapping rings" begin
    poly1 = Polygon(
        [(0.0, 0.0), (-5.0, 5.0), (0.0, 10.0), (5.0, 5.0)],
        [[(0.0, 1.0), (-4.0, 5.0), (0.0, 9.0), (4.0, 5.0)]]
        )
    poly2 = translate(poly1, (3.0, 1.0))
    # difference
    regions = difference_geometry(poly1, poly2, alg)
    expected = [
        [(1.5, 2.5), (0.0, 1.0), (-4.0, 5.0), (0.0, 9.0), (0.5, 8.5), (1.0, 9.0), (1.0, 9.0), (0.0, 10.0), (-5.0, 5.0), (0.0, 0.0), (2.0, 2.0)],
        [(5.0, 5.0), (2.5, 2.5), (2.0, 3.0), (4.0, 5.0), (4.0, 5.0), (1.0, 8.0), (1.5, 8.5)],
    ]
    @test are_regions_equal(regions, expected)
    regions = difference_geometry(poly2, poly1, alg)
    expected = [
        [(2.0, 3.0), (1.5, 2.5), (-2.0, 6.0), (0.5, 8.5), (1.0, 8.0), (1.0, 8.0), (-1.0, 6.0)],
        [(3.0, 11.0), (1.0, 9.0), (1.5, 8.5), (3.0, 10.0), (7.0, 6.0), (7.0, 6.0), (3.0, 2.0), (2.5, 2.5), (2.0, 2.0), (3.0, 1.0), (8.0, 6.0)]
    ]
    @test are_regions_equal(regions, expected)
    ## Intersection
    regions = intersect_geometry(poly1, poly2, alg)
    expected = [
        [(1.0, 9.0), (0.5, 8.5), (1.0, 8.0), (1.5, 8.5)],
        [(2.0, 3.0), (1.5, 2.5), (2.0, 2.0), (2.5, 2.5)],
    ]
    @test are_regions_equal(regions, expected)
    ## Union
    regions = union_geometry(poly1, poly2, alg)
    expected = [
        [(-2.0, 6.0), (0.5, 8.5), (0.5, 8.5), (0.0, 9.0), (-4.0, 5.0), (0.0, 1.0), (1.5, 2.5)], # hole
        [(4.0, 5.0), (2.0, 3.0), (-1.0, 6.0), (1.0, 8.0)], # hole
        [(3.0, 10.0), (1.5, 8.5), (1.5, 8.5), (5.0, 5.0), (2.5, 2.5), (3.0, 2.0), (7.0, 6.0)], # hole
        [(3.0, 11.0), (1.0, 9.0), (0.0, 10.0), (-5.0, 5.0), (0.0, 0.0), (2.0, 2.0), (3.0, 1.0), (8.0, 6.0)],
    ]
    @test are_regions_equal(regions, expected)
    ## XOR
    regions = xor_geometry(poly1, poly2, alg)
    expected = [
        [(-2.0, 6.0), (0.5, 8.5), (0.5, 8.5), (0.0, 9.0), (-4.0, 5.0), (0.0, 1.0), (1.5, 2.5)], # hole
        [(-1.0, 6.0), (1.0, 8.0), (0.5, 8.5), (1.0, 9.0), (1.0, 9.0), (0.0, 10.0), (-5.0, 5.0), (0.0, 0.0), (2.0, 2.0), (1.5, 2.5), (2.0, 3.0)],
        [(3.0, 10.0), (1.5, 8.5), (1.5, 8.5), (5.0, 5.0), (2.5, 2.5), (3.0, 2.0), (7.0, 6.0)], # hole
        [(8.0, 6.0), (3.0, 1.0), (2.0, 2.0), (2.5, 2.5), (2.0, 3.0), (4.0, 5.0), (4.0, 5.0), (1.0, 8.0), (1.5, 8.5), (1.0, 9.0), (3.0, 11.0)],
    ]
    @test are_regions_equal(regions, expected)
end

end