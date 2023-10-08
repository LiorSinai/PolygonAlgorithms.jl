using PolygonAlgorithms: rotate

@testset "intersect-numeric" begin

function are_regions_equal(r1::Vector{Vector{T}}, r2::Vector{Vector{T}}) where T
    if length(r1) != length(r2)
        return false
    end
    r1_sets = [PointSet(r) for r in r1]
    r2_sets = [PointSet(r) for r in r2]
    issetequal(r1_sets, r2_sets)
end

@testset "stars" begin 
    poly1 = [
        (0.0, 18.0), (3.0, 5.0), (15.0, 5.0), (5.0, 0.0), (10.0, -12.0), (0.0, -2.0),
        (-10.0, -12.0), (-5.0, 0.0), (-15.0, 5.0), (-3.0, 5.0)
    ]
    ;
    poly2 = PolygonAlgorithms.rotate(poly1, π/1.0, (0.0, 0.0));
    expected = [[
        (-3.0, -5.0), (-7.083333333333333, -5.0), (-5.0, 0.0), (-7.083333333333333, 5.0), (-3.0, 5.0), 
        (0.0, 2.0), (3.0, 5.0), (7.083333333333333, 5.0), (5.0, 0.0), 
        (7.083333333333333, -5.0), (3.0, -5.0), (0.0, -2.0), (-3.0, -5.0)
    ]]

    regions = intersect_geometry(poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(poly2, poly1)
    @test are_regions_equal(regions, expected)
end

@testset "Hilbert curve order 4 - top corner" begin 
    poly1 = [
        (0.09375, 0.65625),(0.09375, 0.71875),(0.03125, 0.71875),(0.03125, 0.78125),(0.03125, 0.84375),
        (0.09375, 0.84375),(0.09375, 0.78125),(0.15625, 0.78125),(0.21875, 0.78125),(0.21875, 0.84375),
        (0.15625, 0.84375),(0.15625, 0.90625),(0.21875, 0.90625),(0.21875, 0.96875),(0.15625, 0.96875),
        (0.09375, 0.96875),(0.09375, 0.90625),(0.03125, 0.90625),(0.03125, 0.96875),(-0.03125, 0.96875),
        (-0.03125, 0.65625),
    ]
    poly2 = [
        (0.09375, 0.65625),(0.03125, 0.65625),(0.03125, 0.71875),(0.03125, 0.78125),(0.09375, 0.78125),
        (0.09375, 0.84375),(0.03125, 0.84375),(0.03125, 0.90625),(0.031250000000000056, 0.96875),
        (0.09375000000000006, 0.96875),(0.09375, 0.90625),(0.15625, 0.90625),
        (0.15625000000000006, 0.96875),(0.21875000000000006, 0.96875),(0.21875, 0.90625),(0.21875, 0.84375),
        (0.15625, 0.84375),(0.15625, 0.78125),(0.21875, 0.78125),(0.28125, 0.78125),(0.34375, 0.78125),
        (0.34375, 0.65625), 
    ]
    expected = [
        [
            (0.09375, 0.65625), (0.03125, 0.65625), (0.03125, 0.71875), (0.03125, 0.78125), (0.03125, 0.71875), (0.09375, 0.71875)
        ],
        [
            (0.09375, 0.78125), (0.09375, 0.84375), (0.03125, 0.84375), (0.03125, 0.90625), (0.03125, 0.96875), (0.03125, 0.90625), 
            (0.09375, 0.90625), (0.09375, 0.96875), (0.09375, 0.90625), (0.15625, 0.90625), (0.15625, 0.96875), (0.21875, 0.96875), 
            (0.21875, 0.90625), (0.15625, 0.90625), (0.15625, 0.84375), (0.21875, 0.84375), (0.15625, 0.84375), (0.15625, 0.78125), 
            (0.21875, 0.78125), (0.15625, 0.78125)
        ],
    ]
    regions = intersect_geometry(poly1, poly2)
    @test are_regions_equal(regions, expected)
    regions = intersect_geometry(poly2, poly1)
    @test are_regions_equal(regions, expected)
end

end