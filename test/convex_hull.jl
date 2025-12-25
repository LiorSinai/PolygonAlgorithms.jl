using PolygonAlgorithms: GiftWrappingAlg, GrahamScanAlg

@testset "convex hulls - $alg" for alg in [
    GiftWrappingAlg(),
    GrahamScanAlg(),
]

@testset "convex hull rectangle" begin
    h, w = 2.0, 3.0
    poly = [
        (0.0, 0.0), (0.0, h), (w, h), (w, 0.0)
    ]
    hull = convex_hull(poly, alg)
    expected = 1:4
    @test issetequal(hull, expected)

    poly = [
        (0.0, 0.0), (0.0, h/2), 
        (0.0, h), (w/2, h),
        (w, h), (w, h/2), 
        (w, 0.0), (w/2, 0.0),
    ]
    hull = convex_hull(poly, alg)
    expected = [1, 3, 5, 7]
    @test issetequal(hull, expected)

    interior = [(w * rand(), h * rand()) for i in 1:10]
    push!(poly, interior...)
    @test issetequal(hull, expected)
end;

@testset "convex hull H" begin
    poly = [
        (2.0, 10.0), (4.0, 10.0), (4.0, 6.0), (6.0, 6.0), (6.0, 10.0), (8.0, 10.0), 
        (8.0, 1.0), (6.0, 1.0), (6.0, 4.0), (4.0, 4.0), (4.0, 1.0), (2.0, 1.0)
    ]
    hull = convex_hull(poly, alg)
    expected = [1, 6, 7, 12]
    @test issetequal(hull, expected)
end;

@testset "convex hull grid" begin
    points = [(i, j) for i in 1.0:10.0 for j in 1.0:10.0]
    hull = convex_hull(points, alg)
    expected = [1, 10, 91, 100]
    @test issetequal(hull, expected)
end

@testset "generic" begin
    points = [
        (0.5, 0.7),
        (0.3, 0.4),
        (0.4, 0.4),
        (0.6, 1.0),
        (0.6, 0.4),
        (0.5, 0.9),
        (0.4, 0.3),
        (0.0, 0.3),
        (0.8, 0.1),
        (0.3, 0.9),
        (0.8, 0.8),
        (0.9, 0.2),
        (0.5, 0.2),
        (0.7, 0.6),
    ]
    ;
    hull = convex_hull(points, alg)
    expected = [4, 8, 9, 10, 11, 12]
    @test issetequal(hull, expected)
end

@testset "convex hull circle" begin
    # worst case for GiftWrappingAlg
    n = 10
    angles = 2π * rand(n)
    points = [(sin(t), cos(t)) for t in angles]
    hull = convex_hull(points, alg)
    expected = 1:n
    @test issetequal(hull, expected)
end

end