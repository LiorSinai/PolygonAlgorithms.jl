@testset "chain segments" begin
    @testset "chain segments - rectangle" begin
        segments = [
            SegmentEvent(((0.0, 0.0), (4.0, -4.0)), true),
            SegmentEvent(((0.0, 0.0), (3.0, 3.0)), true),
            SegmentEvent(((3.0, 3.0), (7.0, -1.0)), true),
            SegmentEvent(((4.0, -4.0), (7.0, -1.0)), true),
        ]
        regions = chain_segments(segments)
        expected = [[
                SegmentEvent(((3.0, 3.0), (7.0, -1.0)), false),
                SegmentEvent(((0.0, 0.0), (3.0, 3.0)),  false),
                SegmentEvent(((0.0, 0.0), (4.0, -4.0)), true),
                SegmentEvent(((0.0, 0.0), (4.0, -4.0)), false),
        ]]
        @test regions == expected
    end

    @testset "close fuzzy" begin
        segments = [
            SegmentEvent(((0.0, 0.0), (4.0, -4.0)), true),
            SegmentEvent(((0.0, 0.0), (3.0, 3.0)), true),
            SegmentEvent(((3.0, 3.0), (7.0, -1.0)), true),
            SegmentEvent(((3.999, -3.999), (7.0, -1.0)), true),
        ]
        regions = chain_segments(segments)
        expected = [[
            SegmentEvent(((3.999, -3.999), (7.0, -1.0)), true),
            SegmentEvent(((3.0, 3.0), (7.0, -1.0)), false),
            SegmentEvent(((0.0, 0.0), (3.0, 3.0)), false),
            SegmentEvent(((0.0, 0.0), (4.0, -4.0)), true),
            SegmentEvent(((0.0, 0.0), (4.0, -4.0)), false),
        ]]
        @test regions == expected
    end

    @testset "chain segments - improper" begin
        segments = [
            SegmentEvent(((0.0, 0.0), (2.0, 2.0)), true),
            SegmentEvent(((2.0, 2.0), (5.0, 1.0)), true),
            SegmentEvent(((2.0, 2.0), (3.0, 5.0)), true),
            SegmentEvent(((3.0, 5.0), (5.0, 1.0)), true),
        ]
        @test_throws AssertionError chain_segments(segments)
    end

    @testset "chain segments - improper hole" begin
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
        @test_throws AssertionError chain_segments(segments)
    end
end