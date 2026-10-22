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
    # rotating colinear points introduces rounding errors
    line = PolygonAlgorithms.rotate([(0.0, 0.0), (1000.0, 0.0), (2000.0, 0.0)], π/4)
    @test !is_clockwise(line)
    @test !is_counter_clockwise(line)
end

@testset "thin polygon" begin
    # only zero area within floating point error: |area| ≤ n * rtol * Σ(|xᵢyᵢ₊₁| + |yᵢxᵢ₊₁|) / 2
    # a thin polygon near the origin has a direction
    thin = [(0.0, 0.0), (10.0, 0.0), (10.0, 1e-15), (0.0, 1e-15)]
    @test is_counter_clockwise(thin)
    @test is_clockwise(reverse(thin))
    @test is_counter_clockwise(PolygonAlgorithms.Polygon(thin))
    # the same polygon far from the origin is within floating point error
    thin_far = PolygonAlgorithms.translate(thin, (1000.0, 1000.0))
    @test !is_counter_clockwise(thin_far)
    @test !is_clockwise(reverse(thin_far))
    # a thicker polygon far from the origin has a direction, unless the relative tolerance is large
    thicker_far = PolygonAlgorithms.translate([(0.0, 0.0), (10.0, 0.0), (10.0, 1e-6), (0.0, 1e-6)], (1000.0, 1000.0))
    @test is_counter_clockwise(thicker_far)
    @test !is_counter_clockwise(thicker_far; rtol=1e-6)
end

end