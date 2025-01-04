@testset "point in polygon" begin

H_polygon = (
    [
        (2.0, 10.0), (4.0, 10.0), (4.0, 6.0), (6.0, 6.0), (6.0, 10.0), (8.0, 10.0), 
        (8.0, 1.0), (6.0, 1.0), (6.0, 4.0), (4.0, 4.0), (4.0, 1.0), (2.0, 1.0)
    ],
    BitArray([
        0  0  0  0  0  0  0  0  0  0  0
        0  0  1  1  1  0  1  1  1  0  0;
        0  0  1  1  1  0  1  1  1  0  0;
        0  0  1  1  1  0  1  1  1  0  0;
        0  0  1  1  1  1  1  1  1  0  0;
        0  0  1  1  1  1  1  1  1  0  0;
        0  0  1  1  1  1  1  1  1  0  0;
        0  0  1  1  1  0  1  1  1  0  0;
        0  0  1  1  1  0  1  1  1  0  0;
        0  0  1  1  1  0  1  1  1  0  0;
        0  0  1  1  1  0  1  1  1  0  0;
        0  0  0  0  0  0  0  0  0  0  0
    ])
    )
L_polygon = (
    [(2.0, 10.0), (4.0, 10.0), (4.0, 3.0), (7.0, 3.0), (7.0, 1.0), (2.0, 1.0)],
    BitArray([
        0  0  0  0  0  0  0  0  0  0  0;
        0  0  1  1  1  1  1  1  0  0  0;
        0  0  1  1  1  1  1  1  0  0  0;
        0  0  1  1  1  1  1  1  0  0  0;
        0  0  1  1  1  0  0  0  0  0  0;
        0  0  1  1  1  0  0  0  0  0  0;
        0  0  1  1  1  0  0  0  0  0  0;
        0  0  1  1  1  0  0  0  0  0  0;
        0  0  1  1  1  0  0  0  0  0  0;
        0  0  1  1  1  0  0  0  0  0  0;
        0  0  1  1  1  0  0  0  0  0  0;
        0  0  0  0  0  0  0  0  0  0  0
    ])
    )
rectangle = (
    [(2.0, 9.0), (8.0, 9.0), (8.0, 2.0), (2.0, 2.0)],
    BitArray([
        0  0  0  0  0  0  0  0  0  0  0;
        0  0  0  0  0  0  0  0  0  0  0;
        0  0  1  1  1  1  1  1  1  0  0;
        0  0  1  1  1  1  1  1  1  0  0;
        0  0  1  1  1  1  1  1  1  0  0;
        0  0  1  1  1  1  1  1  1  0  0;
        0  0  1  1  1  1  1  1  1  0  0;
        0  0  1  1  1  1  1  1  1  0  0;
        0  0  1  1  1  1  1  1  1  0  0;
        0  0  1  1  1  1  1  1  1  0  0;
        0  0  0  0  0  0  0  0  0  0  0;
        0  0  0  0  0  0  0  0  0  0  0;
    ])
    )
pentagon = (
    [(1.0, 7.0), (5.0, 10.0), (9.0, 7.0), (7.0, 1.0), (3.0, 1.0)],
    BitArray([
        0  0  0  0  0  0  0  0  0  0  0;
        0  0  0  1  1  1  1  1  0  0  0;
        0  0  0  1  1  1  1  1  0  0  0;
        0  0  0  1  1  1  1  1  0  0  0;
        0  0  1  1  1  1  1  1  1  0  0;
        0  0  1  1  1  1  1  1  1  0  0;
        0  0  1  1  1  1  1  1  1  0  0;
        0  1  1  1  1  1  1  1  1  1  0;
        0  0  0  1  1  1  1  1  0  0  0;
        0  0  0  0  1  1  1  0  0  0  0;
        0  0  0  0  0  1  0  0  0  0  0;
        0  0  0  0  0  0  0  0  0  0  0
    ])
    )
skew_H = (
    [
        (2.0, 10.0), (4.0, 10.0), (4.0, 6.0), (6.0, 6.0), (6.0, 10.0), 
        (8.0, 10.0), (8.0, 1.0), (6.0, 1.0), (6.0, 4.0), (4.0, 4.0), (4.0, 2.0), (2.0, 2.0)
    ],
    BitArray([
        0  0  0  0  0  0  0  0  0  0  0
        0  0  0  0  0  0  1  1  1  0  0
        0  0  1  1  1  0  1  1  1  0  0
        0  0  1  1  1  0  1  1  1  0  0
        0  0  1  1  1  1  1  1  1  0  0
        0  0  1  1  1  1  1  1  1  0  0
        0  0  1  1  1  1  1  1  1  0  0
        0  0  1  1  1  0  1  1  1  0  0
        0  0  1  1  1  0  1  1  1  0  0
        0  0  1  1  1  0  1  1  1  0  0
        0  0  1  1  1  0  1  1  1  0  0
        0  0  0  0  0  0  0  0  0  0  0
    ])
    )

@testset "centre points" begin
    @test contains(rectangle[1], (5.0, 5.0),)

    @test contains(H_polygon[1], (3.0, 8.0))
    @test contains(H_polygon[1], (7.0, 8.0))
    @test contains(H_polygon[1], (5.0, 5.0))

    @test contains(L_polygon[1], (3.0, 7.0))
    @test contains(L_polygon[1], (6.0, 2.0))
    @test contains(L_polygon[1], (3.0, 2.0))

    @test contains(pentagon[1], (5.0, 5.0))
    @test contains(pentagon[1], (5.0, 8.0))
end;

@testset "outside points" begin
    @test !contains(H_polygon[1], (5.0, 3.0))
    @test !contains(H_polygon[1], (5.0, 8.0))
    @test !contains(H_polygon[1], (1.0, 8.0))

    @test !contains(L_polygon[1], (8.0, 2.0))
    @test !contains(L_polygon[1], (6.0, 8.0))

    @test !contains(rectangle[1], (5.0, 10.0))
    @test !contains(rectangle[1], (1.0, 5.0))
    @test !contains(rectangle[1], (10.0, 5.0))

    @test !contains(pentagon[1], (2.0, 9.0))
    @test !contains(pentagon[1], (8.0, 9.0))
end;

@testset "horiztonal edges" begin
    @test contains(H_polygon[1], (5.0, 6.0))
    @test !contains(H_polygon[1], (5.0, 10.0))
    @test contains(H_polygon[1], (3.0, 6.0))

    @test !contains(L_polygon[1], (8.0, 3.0))
    @test contains(L_polygon[1], (3.0, 3.0))
    @test contains(L_polygon[1], (6.0, 3.0))

    @test !contains(rectangle[1], (5.0, 10.0))
    @test !contains(rectangle[1], (1.0, 5.0))
    @test !contains(rectangle[1], (10.0, 5.0))

    @test !contains(pentagon[1], (1.0, 1.0))
    @test contains(pentagon[1], (5.0, 1.0))

    @test !contains(skew_H[1], (1.0, 2.0))
    @test !contains(skew_H[1], (5.0, 2.0))
    @test contains(skew_H[1], (7.0, 2.0))
end;

@testset "vertex intersections" begin
    @test !contains(pentagon[1], (2.0, 10.0))
    @test contains(pentagon[1], (5.0, 10.0))
    @test !contains(pentagon[1], (8.0, 10.0))

    @test !contains(pentagon[1], (0.0, 7.0))
    @test contains(pentagon[1], (1.0, 7.0))
    @test contains(pentagon[1], (5.0, 7.0))
    @test contains(pentagon[1], (9.0, 7.0))
    @test !contains(pentagon[1], (10.0, 7.0))
end;

@testset "full grids" begin
    polygons = [
        H_polygon, 
        L_polygon,
        rectangle,
        pentagon,
        skew_H
    ]
    for (polygon, expected) in polygons
        grid_ =  BitArray(undef, 12, 11)
        for x in 0:10
            for y in 0:11
                grid_[y + 1, x + 1] = contains(polygon, float.((x, y)))
            end
        end
        @test expected == grid_
    end
end;

end