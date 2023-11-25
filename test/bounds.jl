@testset "bounds" begin
    @testset "random polygon" begin
        polygon = [
            (2.0, 1.0), (4.0, 6.0), (8.0, -3.0), (3.0, -6.0)
        ]
        rect = bounds(polygon)
        @test rect == (2.0, -6.0, 8.0, 6.0)
    end

    @testset "box" begin
        polygon = [
            (1.5, 0.5), (1.5, 9.5), (10.5, 9.5), (10.5, 0.5)
        ]
        rect = bounds(polygon)
        @test rect == (1.5, 0.5, 10.5, 9.5)
    end

    @testset "line" begin
        polygon = [
            (1.5, 0.5), (10.5, 7.5)
        ]
        rect = bounds(polygon)
        @test rect == (1.5, 0.5, 10.5, 7.5)
    end

    @testset "point" begin
        polygon = [
            (1.5, 0.5)
        ]
        rect = bounds(polygon)
        @test rect == (1.5, 0.5, 1.5, 0.5)
    end

    @testset "empty" begin
        polygon = Tuple{Float64, Float64}[]
        @test_throws BoundsError bounds(polygon)
    end

    @testset "multi polygon" begin
        polygons = [
            [(2.0, 1.0), (4.0, 6.0), (8.0, -3.0), (3.0, -6.0)],
            [(5.0, -10.0)],
            [(-2.0, 2.0), (2.0, 2.0), (-1.0, 0.0), (-0.5, 1.5)]
        ]
        rect = bounds(polygons)
        @test rect == (-2.0, -10.0, 8.0, 6.0)
    end
end
