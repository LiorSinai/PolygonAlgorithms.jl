using Test
using PolygonAlgorithms: AnnotatedSegment, SegmentAnnotations, Polygon
using PolygonAlgorithms: is_hole, match_holes_polygons, events_to_paths

@testset "mappings" begin
    reverse_chain(vec::Vector{<:AnnotatedSegment}) = reverse!([reverse(x) for x in vec])

    @testset "classify holes" begin
        @testset "square hole" begin
            segments = [
                AnnotatedSegment((2.0, 2.0), (2.0, 1.0), SegmentAnnotations(true, false)),
                AnnotatedSegment((2.0, 1.0), (1.0, 1.0), SegmentAnnotations(false, true)),
                AnnotatedSegment((1.0, 1.0), (1.0, 2.0), SegmentAnnotations(false, true)),
                AnnotatedSegment((1.0, 2.0), (2.0, 2.0), SegmentAnnotations(true, false)),
            ]
            @test is_hole(segments)
            @test is_hole(segments, false)
            @test is_hole(reverse_chain(segments))
            # throws error if direction cannot be calculated
            @test_throws AssertionError is_hole(segments[[3, 1, 2, 4]])
            # but passes if is given
            @test is_hole(segments[[3, 1, 2, 4]], false)
        end

        @testset "square interior" begin
            segments = [
                AnnotatedSegment((2.0, 2.0), (1.0, 2.0), SegmentAnnotations(false, true)),
                AnnotatedSegment((1.0, 2.0), (1.0, 1.0), SegmentAnnotations(true, false)),
                AnnotatedSegment((1.0, 1.0), (2.0, 1.0), SegmentAnnotations(true, false)),
                AnnotatedSegment((2.0, 1.0), (2.0, 2.0), SegmentAnnotations(false, true)),
            ]
            @test !is_hole(segments)
            @test !is_hole(segments, true)
            @test !is_hole(reverse_chain(segments))
        end

        @testset "concave star" begin
            segments = [
                AnnotatedSegment(((1.0, 1.0), (2.0, 2.0)), SegmentAnnotations(true, false)),
                AnnotatedSegment(((2.0, 2.0), (3.0, 1.0)), SegmentAnnotations(true, false)),
                AnnotatedSegment(((3.0, 1.0), (2.5, 2.0)), SegmentAnnotations(false, true)),
                AnnotatedSegment(((2.5, 2.0), (3.5, 3.0)), SegmentAnnotations(true, false)),
                AnnotatedSegment(((3.5, 3.0), (2.5, 2.5)), SegmentAnnotations(false, true)),
                AnnotatedSegment(((2.5, 2.5), (1.5, 3.0)), SegmentAnnotations(false, true)),
                AnnotatedSegment(((1.5, 3.0), (2.0, 2.5)), SegmentAnnotations(true, false)),
                AnnotatedSegment(((2.0, 2.5), (1.0, 1.0)), SegmentAnnotations(false, true)),
            ]
            @test !is_hole(segments)
            @test !is_hole(segments, true)
            @test !is_hole(reverse_chain(segments))
        end

        @testset "concave star hole" begin
            segments = [
                AnnotatedSegment(((1.0, 1.0), (2.0, 2.5)), SegmentAnnotations(true, false)),
                AnnotatedSegment(((2.0, 2.5), (1.5, 3.0)), SegmentAnnotations(false, true)),
                AnnotatedSegment(((1.5, 3.0), (2.5, 2.5)), SegmentAnnotations(true, false)),
                AnnotatedSegment(((2.5, 2.5), (3.5, 3.0)), SegmentAnnotations(true, false)),
                AnnotatedSegment(((3.5, 3.0), (2.5, 2.0)), SegmentAnnotations(false, true)),
                AnnotatedSegment(((2.5, 2.0), (3.0, 1.0)), SegmentAnnotations(true, false)),
                AnnotatedSegment(((3.0, 1.0), (2.0, 2.0)), SegmentAnnotations(false, true)),
                AnnotatedSegment(((2.0, 2.0), (1.0, 1.0)), SegmentAnnotations(false, true)),
            ]
            @test is_hole(segments)
            @test is_hole(segments, false)
            @test is_hole(reverse_chain(segments))
        end
    end

    @testset "match holes" begin
        @testset "two squares with holes" begin
            polygons = [
                Polygon([(3.0, 3.0), (3.0, 0.0), (0.0, 0.0), (0.0, 3.0)]),
                Polygon([(7.0, 3.0), (7.0, 0.0), (4.0, 0.0), (4.0, 3.0)]),
            ]
            holes = [
                [(2.0, 2.0), (2.0, 1.0), (1.0, 1.0), (1.0, 2.0)],
                [(6.0, 2.0), (6.0, 1.0), (5.0, 1.0), (5.0, 2.0)],
            ]
            parents = match_holes_polygons(polygons, holes)
            @test parents == [1, 2]
        end

        @testset "nested squares with holes" begin
            poly1 = Polygon(
                [(0.0, 0.0), (0.0, 7.0), (7.0, 7.0), (7.0, 0.0)]
            )
            poly2 = Polygon(
                [(2.0, 2.0), (2.0, 5.0), (5.0, 5.0), (5.0, 2.0)];
            )
            holes = [
                [(1.0, 1.0), (1.0, 6.0), (6.0, 6.0), (6.0, 1.0)],
                [(3.0, 3.0), (3.0, 4.0), (4.0, 4.0), (4.0, 3.0)],
            ]
            # order does not matter
            parents = match_holes_polygons([poly1, poly2], holes)
            @test parents == [1, 2]
            parents = match_holes_polygons([poly2, poly1], holes)
            @test parents == [2, 1]
        end
    end

    @testset "events to path" begin
        @testset "diamond" begin
            segments = [
                SegmentEvent(((0.0, 0.0), (4.0, -4.0)), true),
                SegmentEvent(((0.0, 0.0), (3.0, 3.0)), true),
                SegmentEvent(((3.0, 3.0), (7.0, -1.0)), true),
                SegmentEvent(((4.0, -4.0), (7.0, -1.0)), true),
            ]
            exteriors, holes = events_to_paths(segments)
            expected = [
                [(4.0, -4.0), (7.0, -1.0), (3.0, 3.0), (0.0, 0.0)]
            ]
            @test exteriors == expected
            @test isempty(holes)
        end

        @testset "small gap" begin
            segments = [
                SegmentEvent(((0.0, 0.0), (4.0, -4.0)), true),
                SegmentEvent(((0.0, 0.0), (3.0, 3.0)), true),
                SegmentEvent(((3.0, 3.0), (7.0, -1.0)), true),
                SegmentEvent(((3.999, -3.999), (7.0, -1.0)), true),
            ]
            # open chain, gap not closed
            exteriors, holes = events_to_paths(segments)
            expected = [
                [
                    (0.0, 0.0), (3.0, 3.0), (7.0, -1.0), (3.999, -3.999),
                    (7.0, -1.0), (3.0, 3.0), (0.0, 0.0), (4.0, -4.0)
                ]
            ]
            @test exteriors == expected
            @test isempty(holes)
            # gap closed
            exteriors, holes = events_to_paths(segments; digits=2)
            expected = [
                [(4.0, -4.0), (7.0, -1.0), (3.0, 3.0), (0.0, 0.0)]
            ]
            @test exteriors == expected
            @test isempty(holes)
        end

        @testset "improper" begin
            segments = [
                SegmentEvent(((0.0, 0.0), (2.0, 2.0)), true),
                SegmentEvent(((2.0, 2.0), (5.0, 1.0)), true),
                SegmentEvent(((2.0, 2.0), (3.0, 5.0)), true),
                SegmentEvent(((3.0, 5.0), (5.0, 1.0)), true),
            ]
            exteriors, holes = events_to_paths(segments)
            expected = [
                [(2.0, 2.0), (0.0, 0.0), (2.0, 2.0), (5.0, 1.0), (3.0, 5.0)]
            ]
            @test exteriors == expected
            @test isempty(holes)
        end

        @testset "improper with hole" begin
            segments = [
                SegmentEvent(((0.0, 0.0), (0.0, 4.0)), true, true, SegmentAnnotations(false, true)),
                SegmentEvent(((2.0, 1.0), (2.0, 2.0)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((0.0, 4.0), (2.0, 2.0)), true, true, SegmentAnnotations(true, true) ),
                SegmentEvent(((2.0, 1.0), (3.0, 3.0)), true, true, SegmentAnnotations(false, true)),
                SegmentEvent(((2.0, 2.0), (3.0, 3.0)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((0.0, 0.0), (5.0, 1.0)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((5.0, 1.0), (5.0, 4.0)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((0.0, 4.0), (5.0, 4.0)), true, true, SegmentAnnotations(false, true)),
            ]
            exteriors, holes = events_to_paths(segments)
            expected = [
                [(0.0, 4.0), (0.0, 0.0), (5.0, 1.0), (5.0, 4.0)]
            ]
            @test exteriors == expected
            expected = [
                [(2.0, 1.0), (2.0, 2.0), (3.0, 3.0)]
            ]
            @test holes == expected
        end
    end
end
