@testset verbose = true "point in polygon" begin

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
    @test point_in_polygon((5.0, 5.0), rectangle[1])

    @test point_in_polygon((3.0, 8.0), H_polygon[1])
    @test point_in_polygon((7.0, 8.0), H_polygon[1])
    @test point_in_polygon((5.0, 5.0), H_polygon[1])

    @test point_in_polygon((3.0, 7.0), L_polygon[1])
    @test point_in_polygon((6.0, 2.0), L_polygon[1])
    @test point_in_polygon((3.0, 2.0), L_polygon[1])

    @test point_in_polygon((5.0, 5.0), pentagon[1])
    @test point_in_polygon((5.0, 8.0), pentagon[1])
end;

@testset "outside points" begin
    @test !point_in_polygon((5.0, 3.0), H_polygon[1])
    @test !point_in_polygon((5.0, 8.0), H_polygon[1])
    @test !point_in_polygon((1.0, 8.0), H_polygon[1])

    @test !point_in_polygon((8.0, 2.0), L_polygon[1])
    @test !point_in_polygon((6.0, 8.0), L_polygon[1])

    @test !point_in_polygon((5.0, 10.0), rectangle[1])
    @test !point_in_polygon((1.0, 5.0), rectangle[1])
    @test !point_in_polygon((10.0, 5.0), rectangle[1])

    @test !point_in_polygon((2.0, 9.0), pentagon[1])
    @test !point_in_polygon((8.0, 9.0), pentagon[1])
end;

@testset "horiztonal edges" begin
    @test point_in_polygon((5.0, 6.0), H_polygon[1])
    @test !point_in_polygon((5.0, 10.0), H_polygon[1])
    @test point_in_polygon((3.0, 6.0), H_polygon[1])

    @test !point_in_polygon((8.0, 3.0), L_polygon[1])
    @test point_in_polygon((3.0, 3.0), L_polygon[1])
    @test point_in_polygon((6.0, 3.0), L_polygon[1])

    @test !point_in_polygon((5.0, 10.0), rectangle[1])
    @test !point_in_polygon((1.0, 5.0), rectangle[1])
    @test !point_in_polygon((10.0, 5.0), rectangle[1])

    @test !point_in_polygon((1.0, 1.0), pentagon[1])
    @test point_in_polygon((5.0, 1.0), pentagon[1])

    @test !point_in_polygon((1.0, 2.0), skew_H[1])
    @test !point_in_polygon((5.0, 2.0), skew_H[1])
    @test point_in_polygon((7.0, 2.0), skew_H[1])
end;

@testset "vertix intersections" begin
    @test !point_in_polygon((2.0, 10.0), pentagon[1])
    @test point_in_polygon((5.0, 10.0), pentagon[1])
    @test !point_in_polygon((8.0, 10.0), pentagon[1])

    @test !point_in_polygon((0.0, 7.0), pentagon[1])
    @test point_in_polygon((1.0, 7.0), pentagon[1])
    @test point_in_polygon((5.0, 7.0), pentagon[1])
    @test point_in_polygon((9.0, 7.0), pentagon[1])
    @test !point_in_polygon((10.0, 7.0), pentagon[1])
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
                grid_[y + 1, x + 1] = point_in_polygon(float.((x, y)), polygon)
            end
        end
        @test expected == grid_
    end
end;

end