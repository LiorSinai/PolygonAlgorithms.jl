using Test
using PolygonAlgorithms: points_to_matrix, matrix_to_points
using PolygonAlgorithms: compress_cyclic, cyclic_equality

@testset "utils" begin

@testset "points to matrix" begin
    points = [(1.0, 3.0), (2.0, 4.0)]
    m = points_to_matrix(points)
    @test m == [1.0 2.0 ; 3.0 4.0]
end

@testset "matrix to points" begin
    m = [1.0 2.0 ; 3.0 4.0]
    points = matrix_to_points(m)
    @test points == [(1.0, 3.0), (2.0, 4.0)]
end

@testset "cyclic equality" begin
    @test cyclic_equality([1, 2, 3], [1, 2, 3])
    @test cyclic_equality([1, 2, 3], [2, 3, 1])
    @test cyclic_equality([1, 2, 3], [3, 1, 2])
    @test !cyclic_equality([1, 2, 3], [1, 3, 2])
    @test !cyclic_equality([1, 2, 3], [3, 2, 1])
    @test cyclic_equality([1, 2, 1, 3], [2, 1, 3, 1])
    @test cyclic_equality([], [])
end

@testset "compress cyclic" begin
    a = [1, 2, 3, 2]
    @test compress_cyclic(a) == a
    @test compress_cyclic(Int[]) == Int[]
    a = [1, 2, 2, 3, 1, 1, 2, 2, 1, 1]
    @test compress_cyclic(a) == [1, 2, 3, 1, 2]
    a = [1, 2, 2, 3, 1, 1, 2, 2, 1, 1, 4]
    @test compress_cyclic(a) == [1, 2, 3, 1, 2, 1, 4]
end

end