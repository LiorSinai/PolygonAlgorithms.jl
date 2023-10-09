
@testset "PointSet" begin
    s = PointSet()
    @test isempty(s)
    
    # raw inputs. Not recommended
    s = PointSet(Dict((1.0, 2.0)=>nothing, (1.0, 3.0)=>nothing), 6)
    @test !isempty(s)
    @test (1.0, 2.0) in s
    @test !((2.0, 2.0) in s)

    s = PointSet([(1.123, 2.124), (1.54, 1.55), (-0.0, 0.0)]; digits=1)
    @test !isempty(s)
    @test (1.1, 2.1) in s
    @test (1.5, 1.6) in s
    @test !((1.5, 1.5) in s)
    @test (0.0, -0.0) in s
    @test (0.0, 0.0) in s
    @test (-0.0, -0.0) in s
    @test (-1e-2, 1e-5) in s

    s = PointSet([(1e-7, 0.0), (1.0, π)])
    @test (0.0, 0.0) in s
    @test (1.0, 3.141593) in s
end
