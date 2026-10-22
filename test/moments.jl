@testset "moments" begin

@testset "area rectangle" begin
    h, w = 2.0, 3.0
    poly = [
        (0.0, 0.0), (0.0, h), (w, h), (w, 0.0)
    ]
    area = area_polygon(poly)
    expected = h * w
    @test area == expected

    centroid = centroid_polygon(poly)
    expected = (w/2, h/2)
    @test centroid == expected
end;

@testset "area triangle" begin
    h, b1, b2 = 3.0, 2.0, 1.5
    b = b1 + b2
    poly = [
        (0.0, 0.0), (b1, h), (b, 0.0)
    ]
    area = area_polygon(poly)
    expected = 0.5 * h * b
    @test area == expected

    centroid = centroid_polygon(poly)
    expected = ((b1 + b)/3, h/3)
    @test centroid == expected
end;

@testset "area H" begin
    poly = [
        (2.0, 10.0), (4.0, 10.0), (4.0, 6.0), (6.0, 6.0), (6.0, 10.0), (8.0, 10.0), 
        (8.0, 1.0), (6.0, 1.0), (6.0, 4.0), (4.0, 4.0), (4.0, 1.0), (2.0, 1.0)
    ]
    area = area_polygon(poly)
    expected = 2 * 9 + 2 * 2 + 2 * 9
    @test area == expected

    centroid = centroid_polygon(poly)
    expected = (5.0, 5.45)
    @test centroid == expected
end;

@testset "direction" begin
    square = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    @test is_counter_clockwise(square)
    @test !is_clockwise(square)
    @test is_clockwise(reverse(square))
    # far from the origin
    square_far = PolygonAlgorithms.translate(square, (1e6, 1e6))
    @test is_counter_clockwise(square_far)
    @test is_clockwise(reverse(square_far))
    # Float32
    square32 = [(0.0f0, 0.0f0), (1.0f0, 0.0f0), (1.0f0, 1.0f0), (0.0f0, 1.0f0)]
    @test is_counter_clockwise(square32)
    @test is_clockwise(reverse(square32))
    # Polygon only checks the exterior
    polygon = PolygonAlgorithms.Polygon(square; holes=[[(0.25, 0.25), (0.25, 0.75), (0.75, 0.75), (0.75, 0.25)]])
    @test is_counter_clockwise(polygon)
    @test !is_clockwise(polygon)
end

@testset "zero area" begin
    # colinear points
    line = [(0.0, 0.0), (1.0, 1.0), (2.0, 2.0)]
    @test !is_clockwise(line)
    @test !is_counter_clockwise(line)
    # a line that loops back on itself
    # The first moment is not exactly zero (±1.8e-15) and its sign depends on the Julia version.
    loop = [(0.0, 0.0), (3.0, 3.0), (7.0, -1.0), (3.999, -3.999), (7.0, -1.0), (3.0, 3.0), (0.0, 0.0), (4.0, -4.0)]
    @test !is_clockwise(loop)
    @test !is_counter_clockwise(loop)
    # even with a very small atol, the error scales with the perimeter
    @test !is_clockwise(loop; atol=eps(Float64))
    @test !is_counter_clockwise(loop; atol=eps(Float64))
end

@testset "thin polygon" begin
    # zero area if thinner than atol: |area| ≤ atol * perimeter / 2
    # width=1e-7: area=1e-6, Manhattan perimeter≈20, so the error is 1e-5 for atol=1e-6
    thin = [(0.0, 0.0), (10.0, 0.0), (10.0, 1e-7), (0.0, 1e-7)]
    @test !is_counter_clockwise(thin)
    @test !is_clockwise(reverse(thin))
    # with a smaller atol it has a direction
    @test is_counter_clockwise(thin; atol=1e-8)
    @test is_clockwise(reverse(thin); atol=1e-8)
    @test is_counter_clockwise(thin; atol=eps(Float64))
    # Polygon passes atol to the exterior
    @test !is_counter_clockwise(PolygonAlgorithms.Polygon(thin))
    @test is_counter_clockwise(PolygonAlgorithms.Polygon(thin); atol=1e-8)
    # width=3e-6 is thicker than atol: area=3e-5 > 1e-5
    thicker = [(0.0, 0.0), (10.0, 0.0), (10.0, 3e-6), (0.0, 3e-6)]
    @test is_counter_clockwise(thicker)
    @test is_clockwise(reverse(thicker))
    # the Manhattan perimeter is up to √2 times larger for diagonal edges, so this
    # would fail if the Manhattan perimeter was used
    # area≈1.7e-5, perimeter≈28.3, Manhattan perimeter≈40, error≈2e-5
    d = 1.2e-6 / sqrt(2)
    thin_diagonal = [(0.0, 0.0), (10.0, 10.0), (10.0 - d, 10.0 + d), (-d, d)]
    @test is_counter_clockwise(thin_diagonal)
    # but the same rectangle aligned with the axes is thicker than atol
    # area≈1.7e-5, perimeter=Manhattan perimeter≈28.3, error≈1.4e-5
    thin_aligned = [(0.0, 0.0), (10.0 * sqrt(2), 0.0), (10.0 * sqrt(2), 1.2e-6), (0.0, 1.2e-6)]
    @test is_counter_clockwise(thin_aligned)
end

end