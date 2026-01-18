using PolygonAlgorithms
using PolygonAlgorithms: MartinezRuedaAlg, Polygon, PointSet, translate

function are_equivalent(p1::Polygon, p2::Polygon)
    if length(p1.holes) != length(p2.holes)
        return false
    end
    if !issetequal(PointSet(p1.exterior), PointSet(p2.exterior))
        return false
    end
    for (hole1, hole2) in zip(p1.holes, p2.holes)
        if !issetequal(PointSet(hole1), PointSet(hole2))
            return false
        end
    end
    true
end

function are_equivalent(p1s::Vector{<:Polygon}, p2s::Vector{<:Polygon})
    if length(p1s) != length(p2s)
        return false
    end
    all(are_equivalent(p1, p2) for (p1, p2) in zip(p1s, p2s))
end

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
    expected = [Polygon(poly1.exterior, holes=[poly2.exterior])]
    @test are_equivalent(regions, expected)
    regions = difference_geometry(poly2, poly1, alg)
    @test isempty(regions)
    # intersection
    regions = intersect_geometry(poly1, poly2, alg)
    expected = [poly2]
    @test are_equivalent(regions, expected)
    # union
    regions = union_geometry(poly1, poly2, alg)
    expected = [poly1]
    @test are_equivalent(regions, expected)
    # in this case the polygon is treated as a hole
    regions = union_geometry([poly1, poly2], Polygon{Float64}[])
    expected = [Polygon(poly1.exterior, holes=[poly2.exterior])]
    @test are_equivalent(regions, expected)
    # XOR
    regions = xor_geometry(poly1, poly2, alg)
    expected = [Polygon(poly1.exterior; holes=[poly2.exterior])]
    @test are_equivalent(regions, expected)
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
    expected = [poly1]
    @test are_equivalent(regions, expected)
    # union
    regions = union_geometry(poly1, poly2, alg)
    expected = [poly1]
    @test are_equivalent(regions, expected)
    # XOR
    regions = xor_geometry(poly1, poly2, alg)
    @test isempty(regions)

    ## Overlap
    poly2 = translate(poly1, (1.0, 1.0))
    # difference
    regions = difference_geometry(poly1, poly2, alg)
    expected = [
        Polygon([(3.0, 1.0), (3.0, 0.0), (0.0, 0.0), (0.0, 3.0), (1.0, 3.0), (1.0, 3.0), (1.0, 1.0)])
    ]
    @test are_equivalent(regions, expected)
    regions = difference_geometry(poly2, poly1, alg)
    expected = [
        Polygon([(4.0, 4.0), (4.0, 1.0), (3.0, 1.0), (3.0, 3.0), (3.0, 3.0), (1.0, 3.0), (1.0, 4.0)])
    ]
    @test are_equivalent(regions, expected)
    # intersection
    regions = intersect_geometry(poly1, poly2, alg)
    expected = [Polygon([(3.0, 3.0), (3.0, 1.0), (1.0, 1.0), (1.0, 3.0)])]
    @test are_equivalent(regions, expected)
    ## Union
    regions = union_geometry(poly1, poly2, alg)
    expected = [
        Polygon([(4.0, 4.0), (4.0, 1.0), (3.0, 1.0), (3.0, 0.0), (0.0, 0.0), (0.0, 3.0), (1.0, 3.0), (1.0, 4.0)])
    ]
    @test are_equivalent(regions, expected)
    ## XOR
    regions = xor_geometry(poly1, poly2, alg)
    expected = [
        Polygon([(3.0, 1.0), (3.0, 0.0), (0.0, 0.0), (0.0, 3.0), (1.0, 3.0), (1.0, 3.0), (1.0, 1.0)]),
        Polygon([(4.0, 4.0), (4.0, 1.0), (3.0, 1.0), (3.0, 3.0), (3.0, 3.0), (1.0, 3.0), (1.0, 4.0)]),
    ]
    @test are_equivalent(regions, expected)

    ## No overlap
    poly2 = translate(poly1, (4.0, 3.0))
    # difference
    regions = difference_geometry(poly1, poly2, alg)
    expected = [poly1]
    @test are_equivalent(regions, expected)
    regions = difference_geometry(poly2, poly1, alg)
    expected = [poly2]
    @test are_equivalent(regions, expected)
    # intersection
    regions = intersect_geometry(poly1, poly2, alg)
    @test isempty(regions)
    ## Union
    regions = union_geometry(poly1, poly2, alg)
    expected = [poly1, poly2]
    @test are_equivalent(regions, expected)
    ## XOR
    regions = xor_geometry(poly1, poly2, alg)
    expected = [poly1, poly2]
    @test are_equivalent(regions, expected)
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
        Polygon([(-0.0, 2.0), (-1.0, 1.0), (0.0, 0.0), (1.0, 1.0)]),
        Polygon([(0.0, 12.0), (-1.0, 11.0), (0.0, 10.0), (1.0, 11.0)]),
        Polygon([(2.0, 10.0), (1.0, 9.0), (1.0, 9.0), (4.0, 6.0), (1.0, 3.0), (2.0, 2.0), (6.0, 6.0)]),
    ]
    @test are_equivalent(regions, expected)
    regions = difference_geometry(poly2, poly1, alg)
    expected = [
        Polygon([(0.0, 10.0), (-4.0, 6.0), (-0.0, 2.0), (1.0, 3.0), (-2.0, 6.0), (1.0, 9.0)]),
        Polygon([(2.0, 2.0), (1.0, 1.0), (2.0, 0.0), (3.0, 1.0)]),
        Polygon([(2.0, 12.0), (1.0, 11.0), (2.0, 10.0), (3.0, 11.0)]),
    ]
    @test are_equivalent(regions, expected)
    # intersection
    regions = intersect_geometry(poly1, poly2, alg)
    expected = [
        Polygon([(1.0, 3.0), (2.0, 2.0), (1.0, 1.0), (0.0, 2.0)]),
        Polygon([(0.0, 10.0), (1.0, 11.0), (2.0, 10.0), (1.0, 9.0)]),
    ]
    @test are_equivalent(regions, expected)
    ## Union
    regions = union_geometry(poly1, poly2, alg)
    expected = [
        Polygon(
            [(2.0, 10.0), (3.0, 11.0), (3.0, 11.0), (2.0, 12.0), (1.0, 11.0), (0.0, 12.0), (-1.0, 11.0), (0.0, 10.0), (-4.0, 6.0), (-0.0, 2.0), (-1.0, 1.0), (0.0, 0.0), (1.0, 1.0), (2.0, 0.0), (3.0, 1.0), (2.0, 2.0), (6.0, 6.0)];
            holes=[[(1.0, 9.0), (-2.0, 6.0), (1.0, 3.0), (4.0, 6.0)]],
        )
    ]
    @test are_equivalent(regions, expected)
    ## XOR
    regions = xor_geometry(poly1, poly2, alg)
    expected = [
        Polygon([(0.0, 12.0), (-1.0, 11.0), (0.0, 10.0), (-4.0, 6.0), (-0.0, 2.0), (-1.0, 1.0), (0.0, 0.0), (1.0, 1.0), (-0.0, 2.0), (1.0, 3.0), (-2.0, 6.0), (1.0, 9.0), (0.0, 10.0), (1.0, 11.0)]),
        Polygon([(2.0, 10.0), (3.0, 11.0), (3.0, 11.0), (2.0, 12.0), (1.0, 11.0), (2.0, 10.0), (1.0, 9.0), (1.0, 9.0), (4.0, 6.0), (1.0, 3.0), (2.0, 2.0), (1.0, 1.0), (2.0, 0.0), (3.0, 1.0), (2.0, 2.0), (6.0, 6.0)]),
    ]
    @test are_equivalent(regions, expected)
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
    expected = [Polygon([(4.0, 0.0), (6.0, -2.0), (8.5, -0.0)])]
    @test are_equivalent(regions, expected)
    regions = difference_geometry(poly2, poly1, alg)
    expected = [
        Polygon([(12.0, 3.0), (12.0, 0.0), (11.0, 0.0), (11.0, 2.0), (11.0, 2.0), (8.5, -0.0), (4.0, 0.0), (2.0, 2.0), (0.0, 0.0), (-1.0, 0.0), (-1.0, 3.0)])
    ]
    @test are_equivalent(regions, expected)
    ## Intersection
    regions = intersect_geometry(poly1, poly2, alg)
    expected = [
        Polygon([(4.0, 0.0), (0.0, 0.0), (2.0, 2.0)]),
        Polygon([(8.5, -0.0), (11.0, 0.0), (11.0, 2.0)]),
        Polygon([(4.0, 0.0), (8.5, -0.0)]), # straight line
    ]
    @test are_equivalent(regions, expected)
    ## Union
    regions = union_geometry(poly1, poly2, alg)
    expected = [
        Polygon([(12.0, 3.0), (12.0, 0.0), (8.5, 0.0), (6.0, -2.0), (4.0, 0.0), (-1.0, 0.0), (-1.0, 3.0)])
    ]
    @test are_equivalent(regions, expected)
    ## XOR
    regions = xor_geometry(poly1, poly2, alg)
    expected = [
        Polygon([(12.0, 3.0), (12.0, 0.0), (11.0, 0.0), (11.0, 2.0), (11.0, 2.0), (6.0, -2.0), (2.0, 2.0), (0.0, 0.0), (-1.0, 0.0), (-1.0, 3.0)])
    ]
    @test are_equivalent(regions, expected)
end

@testset "Ring & square" begin
    poly1 = Polygon(
        [(0.0, 0.0), (-5.0, 5.0), (0.0, 10.0), (5.0, 5.0)];
        holes=[[(0.0, 1.0), (-4.0, 5.0), (0.0, 9.0), (4.0, 5.0)]]
        )
    poly2 = Polygon([(3.0, 1.0), (-2.0, 6.0), (3.0, 11.0), (8.0, 6.0)])
    # difference
    regions = difference_geometry(poly1, poly2, alg)
    expected = [
        Polygon([(1.5, 2.5), (0.0, 1.0), (-4.0, 5.0), (0.0, 9.0), (0.5, 8.5), (1.0, 9.0), (1.0, 9.0), (0.0, 10.0), (-5.0, 5.0), (0.0, 0.0), (2.0, 2.0)]),
    ]
    @test are_equivalent(regions, expected)
    regions = difference_geometry(poly2, poly1, alg)
    expected = [
        Polygon([(4.0, 5.0), (1.5, 2.5), (-2.0, 6.0), (0.5, 8.5)]),
        Polygon([(3.0, 11.0), (1.0, 9.0), (1.0, 9.0), (5.0, 5.0), (2.0, 2.0), (3.0, 1.0), (8.0, 6.0)]),
    ]
    @test are_equivalent(regions, expected)
    ## Intersection
    regions = intersect_geometry(poly1, poly2, alg)
    expected = [
        Polygon([(5.0, 5.0), (2.0, 2.0), (1.5, 2.5), (4.0, 5.0), (4.0, 5.0), (0.5, 8.5), (1.0, 9.0)])
    ]
    @test are_equivalent(regions, expected)
    ## Union
    regions = union_geometry(poly1, poly2, alg)
    expected = [
        Polygon(
            [(3.0, 11.0), (1.0, 9.0), (0.0, 10.0), (-5.0, 5.0), (0.0, 0.0), (2.0, 2.0), (3.0, 1.0), (8.0, 6.0)];
            holes=[[(-2.0, 6.0), (0.5, 8.5), (0.5, 8.5), (0.0, 9.0), (-4.0, 5.0), (0.0, 1.0), (1.5, 2.5)]]
        )
    ]
    @test are_equivalent(regions, expected)
    ## XOR
    regions = xor_geometry(poly1, poly2, alg)
    expected = [
        Polygon(
            [(0.5, 8.5), (1.0, 9.0), (1.0, 9.0), (0.0, 10.0), (-5.0, 5.0), (0.0, 0.0), (2.0, 2.0), (1.5, 2.5), (4.0, 5.0)];
            holes=[[(-2.0, 6.0), (0.5, 8.5), (0.5, 8.5), (0.0, 9.0), (-4.0, 5.0), (0.0, 1.0), (1.5, 2.5)]]
        ),
        Polygon(
            [(3.0, 11.0), (1.0, 9.0), (1.0, 9.0), (5.0, 5.0), (2.0, 2.0), (3.0, 1.0), (8.0, 6.0)]
        ),
    ]
    @test are_equivalent(regions, expected)
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
        Polygon([(1.5, 2.5), (0.0, 1.0), (-4.0, 5.0), (0.0, 9.0), (0.5, 8.5), (1.0, 9.0), (1.0, 9.0), (0.0, 10.0), (-5.0, 5.0), (0.0, 0.0), (2.0, 2.0)]),
        Polygon([(5.0, 5.0), (2.5, 2.5), (2.0, 3.0), (4.0, 5.0), (4.0, 5.0), (1.0, 8.0), (1.5, 8.5)]),
    ]
    @test are_equivalent(regions, expected)
    regions = difference_geometry(poly2, poly1, alg)
    expected = [
        Polygon([(2.0, 3.0), (1.5, 2.5), (-2.0, 6.0), (0.5, 8.5), (1.0, 8.0), (1.0, 8.0), (-1.0, 6.0)]),
        Polygon([(3.0, 11.0), (1.0, 9.0), (1.5, 8.5), (3.0, 10.0), (7.0, 6.0), (7.0, 6.0), (3.0, 2.0), (2.5, 2.5), (2.0, 2.0), (3.0, 1.0), (8.0, 6.0)]),
    ]
    @test are_equivalent(regions, expected)
    ## Intersection
    regions = intersect_geometry(poly1, poly2, alg)
    expected = [
        Polygon([(1.0, 9.0), (0.5, 8.5), (1.0, 8.0), (1.5, 8.5)]),
        Polygon([(2.0, 3.0), (1.5, 2.5), (2.0, 2.0), (2.5, 2.5)]),
    ]
    @test are_equivalent(regions, expected)
    ## Union
    regions = union_geometry(poly1, poly2, alg)
    expected = [
        Polygon(
            [(3.0, 11.0), (1.0, 9.0), (0.0, 10.0), (-5.0, 5.0), (0.0, 0.0), (2.0, 2.0), (3.0, 1.0), (8.0, 6.0)];
            holes=[
                [(-2.0, 6.0), (0.5, 8.5), (0.5, 8.5), (0.0, 9.0), (-4.0, 5.0), (0.0, 1.0), (1.5, 2.5)],
                [(4.0, 5.0), (2.0, 3.0), (-1.0, 6.0), (1.0, 8.0)],
                [(3.0, 10.0), (1.5, 8.5), (1.5, 8.5), (5.0, 5.0), (2.5, 2.5), (3.0, 2.0), (7.0, 6.0)],
            ]
        )
    ]
    @test are_equivalent(regions, expected)
    ## XOR
    regions = xor_geometry(poly1, poly2, alg)
    expected = [
        Polygon(
            [(-1.0, 6.0), (1.0, 8.0), (0.5, 8.5), (1.0, 9.0), (1.0, 9.0), (0.0, 10.0), (-5.0, 5.0), (0.0, 0.0), (2.0, 2.0), (1.5, 2.5), (2.0, 3.0)];
            holes=[[(-2.0, 6.0), (0.5, 8.5), (0.5, 8.5), (0.0, 9.0), (-4.0, 5.0), (0.0, 1.0), (1.5, 2.5)]]
        ),
        Polygon(
            [(8.0, 6.0), (3.0, 1.0), (2.0, 2.0), (2.5, 2.5), (2.0, 3.0), (4.0, 5.0), (4.0, 5.0), (1.0, 8.0), (1.5, 8.5), (1.0, 9.0), (3.0, 11.0)];
            holes=[[(3.0, 10.0), (1.5, 8.5), (1.5, 8.5), (5.0, 5.0), (2.5, 2.5), (3.0, 2.0), (7.0, 6.0)]]
        ),
    ]
    @test are_equivalent(regions, expected)
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
    regions = intersect_geometry([elbow], [triangle, rect], alg)
    @test isempty(regions)
    # Union
    expected = [
        Polygon(
            [(3.0, 2.0), (3.0, 1.0), (2.5, 1.0), (2.5, -0.8), (1.0, -0.8), (1.0, -1.0), (0.0, -1.0), (0.0, 2.0)];
            holes=[
                [(1.323077, 0.0), (1.0, -0.4941176), (1.0, -0.0)],
            [(1.976923, 1.0), (1.65, 0.5), (1.0, 0.5), (1.0, 1.0)],
            ]
        )
    ]
    regions = union_geometry([elbow], [triangle, rect], alg)
    @test are_equivalent(regions, expected)
    ## Difference
    expected = [
        Polygon(
            [(3.0, 2.0), (3.0, 1.0), (2.5, 1.0), (2.5, 1.8), (2.5, 1.8), (1.976923, 1.0), (2.5, 1.0), (2.5, -0.8), (1.0, -0.8), (1.0, -0.4941176), (1.0, -0.4941176), (0.8, -0.8), (1.0, -0.8), (1.0, -1.0), (0.0, -1.0), (0.0, 2.0)];
            holes=[
                [(2.0, 0.5), (2.0, 0.0), (1.3230769, 0.0), (1.0, -0.4941176), (1.0, -0.0), (1.0, -0.0), (0.5, 0.0), (0.5, 0.5), (1.0, 0.5), (1.0, 1.0), (1.976923, 1.0), (1.976923, 1.0), (1.65, 0.5)]
            ]
        )
    ]
    regions = difference_geometry([elbow, triangle], [rect], alg)
    @test are_equivalent(regions, expected)

    # XOR
    expected = [
        Polygon(
            [(3.0, 2.0), (3.0, 1.0), (2.5, 1.0), (2.5, 1.8), (2.5, 1.8), (1.976923, 1.0), (2.5, 1.0), (2.5, -0.8), (1.0, -0.8), (1.0, -0.4941176), (1.0, -0.4941176), (0.8, -0.8), (1.0, -0.8), (1.0, -1.0), (0.0, -1.0), (0.0, 2.0)];
            holes=[
                [(1.3230769, 0.0), (1.0, -0.4941176), (1.0, -0.0), (1.0, -0.0), (0.5, 0.0), (0.5, 0.5), (1.0, 0.5), (1.0, 0.5), (1.0, -0.0)],
                [(2.0, 0.5), (2.0, 0.0), (1.3230769, 0.0), (1.65, 0.5), (1.65, 0.5), (1.0, 0.5), (1.0, 1.0), (1.976923, 1.0), (1.976923, 1.0), (1.65, 0.5)],
            ]
        ),
    ]
    regions = xor_geometry([elbow], [triangle, rect], alg)
    @test are_equivalent(regions, expected)
end

end