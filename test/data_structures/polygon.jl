
using PolygonAlgorithms: Polygon, fully_contains, validate_polygon
using PolygonAlgorithms: x_coords, y_coords, translate, rotate

@testset "Polygon data structure" begin

@testset "area-centroid holes" begin
    polygon = Polygon(
        [(0.0, 0.0), (20.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
        [[(1.0, 8.0), (7.0, 8.0), (7.0, 3.0), (1.0, 3.0)]]
    )
    @test first_moment(polygon) == 120
    @test area_polygon(polygon) == 120

    @test centroid_polygon(polygon) == (8.722222222222223, 4.180555555555556)

    @test bounds(polygon) == (0.0, 0.0, 20.0, 10.0)

    @test !is_clockwise(polygon)
    @test is_counter_clockwise(polygon)
    @test is_clockwise(polygon.holes[1])

    # swap hole and exterior
    polygon = Polygon(
        [(1.0, 8.0), (7.0, 8.0), (7.0, 3.0), (1.0, 3.0)],
        [[(0.0, 0.0), (20.0, 0.0), (10.0, 10.0), (0.0, 10.0)]],
    ) # skips validation
    @test area_polygon(polygon) == -120
    @test centroid_polygon(polygon) == (8.722222222222223, 4.180555555555556)

    # exterior == hole
    polygon = Polygon(
        [(0.0, 0.0), (20.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
        [[(0.0, 0.0), (20.0, 0.0), (10.0, 10.0), (0.0, 10.0)]],
    ) # skips validation
    @test area_polygon(polygon) == 0.0
    @test all(isnan.(centroid_polygon(polygon)))
end

@testset "coordinates" begin
    polygon = Polygon(
        [(0.0, 0.0), (20.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
        [[(1.0, 8.0), (7.0, 8.0), (7.0, 3.0), (1.0, 3.0)]]
    )
    expected = [
        [0.0, 20.0, 10.0, 0.0], [1.0, 7.0, 7.0, 1.0]
    ]
    @test x_coords(polygon) == expected
    expected = [
        [0.0, 0.0, 10.0, 10.0], [8.0, 8.0, 3.0, 3.0]
    ]
    @test y_coords(polygon) == expected
end

@testset "translate and rotate" begin
    polygon = Polygon(
        [(0.0, 0.0), (20.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
        [[(1.0, 8.0), (7.0, 8.0), (7.0, 3.0), (1.0, 3.0)]]
    )
    
    poly2 = translate(polygon, (1.0, 2.5))
    expected = Polygon(
        [(1.0, 2.5), (21.0, 2.5), (11.0, 12.5), (1.0, 12.5)],
        [[(2.0, 10.5), (8.0, 10.5), (8.0, 5.5), (2.0, 5.5)]]
    )
    @test poly2 == expected

    poly2 = rotate(polygon, pi/2)
    expected = [
        (0.0, 0.0), (0.0, 20.0), (-10.0, 10.0), (-10.0, 0.0),
        (-8.0, 1.0), (-8.0, 7.0), (-3.0, 7.0), (-3.0, 1.0),
    ]
    points = vcat(poly2.exterior, poly2.holes...)
    @test all([all(isapprox.(p1, p2, atol=1e-6)) for (p1, p2) in zip(points, expected)])

    poly2 = rotate(polygon, pi/2)
    expected = [
        (0.0, 0.0), (0.0, 20.0), (-10.0, 10.0), (-10.0, 0.0),
        (-8.0, 1.0), (-8.0, 7.0), (-3.0, 7.0), (-3.0, 1.0),
    ]
    points = vcat(poly2.exterior, poly2.holes...)
    @test all([all(isapprox.(p1, p2, atol=1e-6)) for (p1, p2) in zip(points, expected)])

    poly2 = rotate(polygon, pi/3, (20.0, 0.0))
    expected = [
        (10.0, -17.320508), (20.0, 0.0), (6.339746, -3.660254), (1.339745, -12.320508),
        (3.571797, -12.454483), (6.571797, -7.2583302), (10.901923, -9.758330), (7.901924, -14.954483),
    ]
    points = vcat(poly2.exterior, poly2.holes...)
    @test all([all(isapprox.(p1, p2, atol=1e-6)) for (p1, p2) in zip(points, expected)])
end


@testset "point in polygon - holes" begin
    polygon = Polygon(
        [(0.0, 0.0), (20.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
        [[(1.0, 3.0), (7.0, 3.0), (7.0, 8.0), (1.0, 8.0)]]
    )
    @test contains(polygon, (10.0, 5.0)) # inside
    @test !contains(polygon, (15.0, 10.0)) # outside
    @test !contains(polygon, (5.0, 5.0)) # in hole
    @test contains(polygon, (5.0, 3.0)) # on hole border
    @test !contains(polygon, (5.0, 3.0); on_border_is_inside=false) # on hole border
end

@testset "polygon in polygon" begin
    exterior = [(0.0, 0.0), (20.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
    hole = [(1.0, 3.0), (7.0, 3.0), (7.0, 8.0), (1.0, 8.0)]
    @test fully_contains(exterior, hole) # inside
    # inside on boundary
    polygon2 = translate(hole, (-1.0, 0.0))
    @test fully_contains(exterior, polygon2)
    # intersects
    polygon2 = translate(hole, (-2.0, 0.0))
    @test !fully_contains(exterior, polygon2)
    # outside on boundary
    polygon2 = translate(hole, (-7.0, 0.0))
    @test !fully_contains(exterior, polygon2)
    polygon2 = [(0.0, 3.0), (0.0, 8.0), (-6.0, 8.0)]
    @test !fully_contains(exterior, polygon2)
    # outside
    polygon2 = translate(hole, (-9.0, 0.0))
    @test !fully_contains(exterior, polygon2)
    # inside near boundary
    polygon2 = translate(hole, (4.999, 0.0))
    @test fully_contains(exterior, polygon2)
    # inside touch boundary
    polygon2 = translate(hole, (5.0, 0.0))
    @test fully_contains(exterior, polygon2)
    # outside intersects
    polygon2 = translate(hole, (9.0, 0.0))
    @test !fully_contains(exterior, polygon2)
    # outside
    polygon2 = translate(hole, (10.0, 10.0))
    @test !fully_contains(exterior, polygon2)
end

@testset "polygon in polygon - edge case" begin
    hour_glass = [(0.0, 0.0), (4.0, 5.0), (0.0, 10.0), (10.0, 10.0), (6.0, 5.0), (10.0, 0.0)]
    rectangle = [(2.0, 1.0),(2.0, 8.0), (8.0, 8.0), (8.0, 1.0)]
    @test all(p -> contains(hour_glass, p), rectangle)
    @test !fully_contains(hour_glass, rectangle)
end


@testset "polygon in polygon - self-intersect" begin
    pentagon = [
        (-0.8, 0.0), (0.0, 0.6), (0.8, 0.0), (0.5, -1.0), (-0.5, -1.0)
    ]
    pentagram = pentagon[[1, 3, 5, 2, 4]]
    # touches vertices
    @test_broken fully_contains(pentagon, pentagram) # cannot handle self-intersecting polygons
    pieces = [ # but after dividing into simple polygons
        [(0.0, -0.615384), (-0.5, -1.0), (-0.306201, -0.379844), (-0.8, 0.0), (-0.1875, 0.0), (-0.1875, 0.0), (-0.306201, -0.379844)],
        [(0.0, 0.6), (0.1875, 0.0), (-0.1875, 0.0)],
        [(0.306201, -0.379845), (0.5, -1.0), (0.0, -0.615384), (0.306201, -0.379845), (0.1875, 0.0), (0.8, 0.0)],
    ]
    @test all(p -> fully_contains(pentagon, p), pieces) 
    pentagram2 = [p .* 0.9 for p in pentagram]
    @test_broken fully_contains(pentagon, pentagram2) # cannot handle self-intersecting polygons
    pieces2 = [[p .* 0.9 for p in piece] for piece in pieces]
    @test all(p -> fully_contains(pentagon, p), pieces2) 
end

@testset "validate polygon" begin
    # single point
    @test_throws "exterior requires at least 3 points" validate_polygon([(1.0, 2.0)])
    @test Polygon([(1.0, 2.0)]) isa Polygon # skips validation checks
    # Validation checks
    exterior = [(0.0, 0.0), (20.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
    holes = [[(1.0, 8.0), (7.0, 8.0), (7.0, 3.0), (1.0, 3.0)]]
    @test validate_polygon(exterior, holes=holes)
    @test_throws "Invalid hole 1: requires at least 3 points." validate_polygon(exterior, holes=[[(1.0, 2.0)]])
    polygon2 = translate(holes[1], (-2.0, 0.0))
    @test_throws "Hole 1 intersects with the exterior." validate_polygon(exterior; holes=[polygon2])
    polygon2 = translate(holes[1], (-10.0, 0.0))
    @test_throws "Hole 1 is outside the polygon." validate_polygon(exterior; holes=[polygon2])
    # worst case: no segments intersect and all points on polygon
    @test_throws "Hole 1 is outside the polygon." validate_polygon(exterior; holes=[exterior])
    # self intersect
    self_intersect = [
        (2.0, 3.0), (4.0, 5.0), (8.0, 1.0), (13.0, 5.0), (13.0, 3.0)
    ]
    @test_throws "Invalid exterior: edges self-intersect." validate_polygon(self_intersect; holes=holes)
    @test_throws "Invalid hole 1: edges self-intersect." validate_polygon(exterior; holes=[self_intersect])
end

end