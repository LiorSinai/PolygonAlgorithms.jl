@testset "area-centroids" begin

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

end