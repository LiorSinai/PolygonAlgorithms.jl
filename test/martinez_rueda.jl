using Test
using PolygonAlgorithms: SegmentEvent, SegmentAnnotations
using PolygonAlgorithms: add_annotated_segment!, compare_events, convert_to_event_queue
using PolygonAlgorithms: check_and_divide_intersection!, event_loop!
using PolygonAlgorithms: chain_segments
using PolygonAlgorithms: is_hole, match_holes_polygons, Polygon
using PolygonAlgorithms: BLANK

@testset "Martinez-Rueda algorithm" verbose=false begin

    @testset "almost vertical" begin
        # the original segment was ((0.5345889348336309, 0.9626445530693036), (0.5346819812512368, 0.7153822809806809))
        # it was not vertical, but the two intersecting lines make a tiny vertical segment
        event_queue = [
            SegmentEvent(((0.5346362515873916, 0.8369046444971744), (0.5346819812512368, 0.7153822809806809)), true),
            SegmentEvent(((0.5346362515873916, 0.8369046444971744), (0.702057739662388, 0.63591526606886130)), true),
            SegmentEvent(((0.5346358170779485, 0.8380593132624886), (0.5346362515873916, 0.8369046444971744)), true),
            SegmentEvent(((0.5346358170779485, 0.8380593132624886), (0.6336845198791616, 0.7251516629897194)), true),
        ]
        # the event_loop! modified event_queue[3] with update_end!
        # now, the tail should be inserted after it.
        tail = SegmentEvent(event_queue[3].segment, false)
        idx = searchsortedfirst(event_queue, tail; lt=compare_events)
        @test_broken idx == 5
    end

    @testset "divide" begin
        @testset "divide no intersection" begin
            # setup
            ev1_start = SegmentEvent(((3.0, 3.0), (7.0, -1.0)), true)
            ev1_end = SegmentEvent(((3.0, 3.0), (7.0, -1.0)), false)
            ev1_start.other = ev1_end
            ev1_end.other = ev1_start
            ev2_start = SegmentEvent(((4.0, 5.0), (8.0, 3.0)), true)
            ev2_end = SegmentEvent(((4.0, 5.0), (8.0, 3.0)), false)
            ev2_start.other = ev2_end
            ev2_end.other = ev2_start
            queue = [ev1_start, ev2_start, ev1_end, ev2_end]
            # test
            check_and_divide_intersection!(queue, ev1_start, ev2_start, false)
            @test length(queue) == 4
        end

        @testset "divide intersection" begin
            # setup
            ev1_start = SegmentEvent(((3.0, 3.0), (7.0, -1.0)), true)
            ev1_end = SegmentEvent(((3.0, 3.0), (7.0, -1.0)), false)
            ev1_start.other = ev1_end
            ev1_end.other = ev1_start
            ev2_start = SegmentEvent(((4.0, 0.0), (7.0, 3.0)), true)
            ev2_end = SegmentEvent(((4.0, 0.0), (7.0, 3.0)), false)
            ev2_start.other = ev2_end
            ev2_end.other = ev2_start
            queue = [ev1_start, ev2_start, ev1_end, ev2_end]
            # test
            check_and_divide_intersection!(queue, ev1_start, ev2_start, false)
            expected = [
                ev1_start,
                ev2_start,
                ev2_end,
                ev1_end,
                SegmentEvent(((5.0, 1.0), (7.0, -1.0)), true),
                SegmentEvent(((5.0, 1.0), (7.0, 3.0)), true),
                SegmentEvent(((5.0, 1.0), (7.0, -1.0)), false),
                SegmentEvent(((5.0, 1.0), (7.0, 3.0)), false),
            ]
            @test queue == expected
        end

        @testset "divide along1" begin
            # setup
            ev1_start = SegmentEvent(((3.0, 3.0), (7.0, -1.0)), true)
            ev1_end = SegmentEvent(((3.0, 3.0), (7.0, -1.0)), false)
            ev1_start.other = ev1_end
            ev1_end.other = ev1_start
            ev2_start = SegmentEvent(((4.0, 2.0), (7.0, 5.0)), true)
            ev2_end = SegmentEvent(((4.0, 2.0), (7.0, 5.0)), false)
            ev2_start.other = ev2_end
            ev2_end.other = ev2_start
            queue = [ev1_start, ev2_start, ev1_end, ev2_end]
            # test
            check_and_divide_intersection!(queue, ev1_start, ev2_start, false)
            expected = [
                ev1_start,
                ev1_end,
                SegmentEvent(((4.0, 2.0), (7.0, -1.0)), true), 
                ev2_start,
                SegmentEvent(((4.0, 2.0), (7.0, -1.0)), false), 
                ev2_end
            ]
            @test queue == expected
        end

        @testset "divide coincident" begin
            # setup
            ev1_start = SegmentEvent(((3.0, 3.0), (5.0, 1.0)), true)
            ev1_end = SegmentEvent(((3.0, 3.0), (5.0, 1.0)), false)
            ev1_start.other = ev1_end
            ev1_end.other = ev1_start
            ev2_start = SegmentEvent(((4.0, 2.0), (7.0, -1.0)), true)
            ev2_end = SegmentEvent(((4.0, 2.0), (7.0, -1.0)), false)
            ev2_start.other = ev2_end
            ev2_end.other = ev2_start
            queue = [ev1_start, ev2_start, ev1_end, ev2_end]
            # test
            check_and_divide_intersection!(queue, ev1_start, ev2_start, false)
            @test length(queue) == 4 # no change because ev1 is processed before ev2
            check_and_divide_intersection!(queue, ev2_start, ev1_start, false)
            expected = [
                ev1_start,
                ev1_end,
                SegmentEvent(((4.0, 2.0), (5.0, 1.0)), true),
                ev2_start,
                SegmentEvent(((4.0, 2.0), (5.0, 1.0)), false),
                ev2_end,
                SegmentEvent(((5.0, 1.0), (7.0, -1.0)), true),
                SegmentEvent(((5.0, 1.0), (7.0, -1.0)), false),
            ]
            @test queue == expected
        end

        @testset "divide coincident same start" begin
            # setup
            ev1_start = SegmentEvent(((3.0, 3.0), (5.0, 1.0)), true)
            ev1_end = SegmentEvent(((3.0, 3.0), (5.0, 1.0)), false)
            ev1_start.other = ev1_end
            ev1_end.other = ev1_start
            primary = true
            self_annotations = SegmentAnnotations(true, false)
            ev2_start = SegmentEvent(((3.0, 3.0), (7.0, -1.0)), true, primary, self_annotations)
            ev2_end = SegmentEvent(((3.0, 3.0), (7.0, -1.0)), false, primary, self_annotations)
            ev2_start.other = ev2_end
            ev2_end.other = ev2_start
            queue = [ev1_start, ev2_start, ev1_end, ev2_end]
            # test
            check_and_divide_intersection!(queue, ev1_start, ev2_start, true)
            expected = [
                ev2_start,
                ev2_end,
                SegmentEvent(((5.0, 1.0), (7.0, -1.0)), true, true, self_annotations),
                SegmentEvent(((5.0, 1.0), (7.0, -1.0)), false, true, self_annotations),
            ]
            @test queue == expected
        end

        @testset "divide same intersection" begin
            # setup
            ev1_start = SegmentEvent(((3.0, 3.0), (5.0, 1.0)), true)
            ev1_end = SegmentEvent(((3.0, 3.0), (5.0, 1.0)), false)
            ev1_start.other = ev1_end
            ev1_end.other = ev1_start
            primary = true
            self_annotations = SegmentAnnotations(true, false)
            ev2_start = SegmentEvent(((3.0, 3.0), (5.0, 1.0)), true, primary, self_annotations)
            ev2_end = SegmentEvent(((3.0, 3.0), (5.0, 1.0)), false, primary, self_annotations)
            ev2_start.other = ev2_end
            ev2_end.other = ev2_start
            queue = [ev1_start, ev2_start, ev1_end, ev2_end]
            # test
            check_and_divide_intersection!(queue, ev1_start, ev2_start, true)
            @test queue == [ev2_start, ev2_end]
        end
    end

    @testset "annotations" begin
        @testset "rectangle self annotations" begin
            rectangle = [
                (3.0, 3.0), (7.0, -1.0), (4.0, -4.0), (0.0, 0.0)
            ];
            event_queue = convert_to_event_queue(rectangle)
            annotated_segments = event_loop!(event_queue; self_intersection=true)
            expected = [
                SegmentEvent(((0.0, 0.0), (3.0, 3.0)), true, true, SegmentAnnotations(false, true)),
                SegmentEvent(((0.0, 0.0), (4.0, -4.0)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((4.0, -4.0), (7.0, -1.0)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((3.0, 3.0), (7.0, -1.0)), true, true, SegmentAnnotations(false, true)),
            ]
            @test annotated_segments == expected
        end

        @testset "rectangle horiz self annotations" begin
            rectangle_horiz = [
                (-1.0, 0.0), (-1.0, 3.0), (12.0, 3.0), (12.0, 0.0)
            ];
            event_queue = convert_to_event_queue(rectangle_horiz)
            annotated_segments = event_loop!(event_queue; self_intersection=true)
            expected = [
                SegmentEvent(((-1.0, 0.0), (-1.0, 3.0)), true, true, SegmentAnnotations(false, true)),
                SegmentEvent(((-1.0, 0.0), (12.0, 0.0)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((12.0, 0.0), (12.0, 3.0)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((-1.0, 3.0), (12.0, 3.0)), true, true, SegmentAnnotations(false, true)),
            
            ]
            @test annotated_segments == expected
        end

        @testset "self-intersect self annotations" begin
            self_intersect = [
                (0.0, 0.0), (2.0, 2.0), (6.0, -2.0), (11.0, 2.0), (11.0, 0.0)
            ]
            event_queue = convert_to_event_queue(self_intersect)
            annotated_segments = event_loop!(event_queue; self_intersection=true)
            expected = [
                SegmentEvent(((0.0, 0.0), (2.0, 2.0)), true, true, SegmentAnnotations(false, true)),
                SegmentEvent(((0.0, 0.0), (4.0, 0.0)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((2.0, 2.0), (4.0, 0.0)), true, true, SegmentAnnotations(false, true)),
                SegmentEvent(((4.0, 0.0), (6.0, -2.0)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((6.0, -2.0), (8.5, -0.0)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((4.0, 0.0), (8.5, -0.0)), true, true, SegmentAnnotations(false, true)),
                SegmentEvent(((8.5, -0.0), (11.0, 0.0)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((11.0, 0.0), (11.0, 2.0)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((8.5, -0.0), (11.0, 2.0)), true, true, SegmentAnnotations(false, true)),
            ]
            @test annotated_segments == expected
        end

        @testset "pentagram self annotations" begin
            pentagon = [
                (-0.8, 0.0), (0.0, 0.6), (0.8, 0.0), (0.5, -1.0), (-0.5, -1.0)
            ]
            pentagram = pentagon[[1, 3, 5, 2, 4]]
            event_queue = convert_to_event_queue(pentagram)
            annotated_segments = event_loop!(event_queue; self_intersection=true)
            expected = [
                SegmentEvent(((-0.5, -1.0), (-0.3062015503875969, -0.3798449612403101)), true, true, SegmentAnnotations(false, true)),
                SegmentEvent(((-0.8, 0.0), (-0.3062015503875969, -0.3798449612403101)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((-0.3062015503875969, -0.3798449612403101), (-0.18749999999999997, 5.551115123125783e-17)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((-0.8, 0.0), (-0.18749999999999997, 5.551115123125783e-17)), true, true, SegmentAnnotations(false, true)),
                SegmentEvent(((-0.5, -1.0), (0.0, -0.6153846153846154)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((-0.3062015503875969, -0.3798449612403101), (0.0, -0.6153846153846154)), true, true, SegmentAnnotations(false, true)),
                SegmentEvent(((-0.18749999999999997, 5.551115123125783e-17), (0.0, 0.6)), true, true, SegmentAnnotations(false, true)),
                SegmentEvent(((-0.18749999999999997, 5.551115123125783e-17), (0.18749999999999997, 0.0)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((0.0, 0.6), (0.18749999999999997, 0.0)), true, true, SegmentAnnotations(false, true)),
                SegmentEvent(((0.0, -0.6153846153846154), (0.3062015503875969, -0.37984496124031014)), true, true, SegmentAnnotations(false, true)),
                SegmentEvent(((0.18749999999999997, 0.0), (0.3062015503875969, -0.37984496124031014)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((0.0, -0.6153846153846154), (0.5, -1.0)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((0.3062015503875969, -0.37984496124031014), (0.5, -1.0)), true, true, SegmentAnnotations(false, true)),
                SegmentEvent(((0.3062015503875969, -0.37984496124031014), (0.8, 0.0)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((0.18749999999999997, 0.0), (0.8, 0.0)), true, true, SegmentAnnotations(false, true)),
                ]
            @test annotated_segments == expected
        end

        @testset "rectangles intersect annotations" begin
            rectangle = [
                SegmentEvent(((0.0, 0.0), (3.0, 3.0)), true, true, SegmentAnnotations(false, true)),
                SegmentEvent(((0.0, 0.0), (4.0, -4.0)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((4.0, -4.0), (7.0, -1.0)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((3.0, 3.0), (7.0, -1.0)), true, true, SegmentAnnotations(false, true)),
            ]
            rectangle_horiz = [
                SegmentEvent(((-1.0, 0.0), (-1.0, 3.0)), true, false, SegmentAnnotations(false, true)),
                SegmentEvent(((-1.0, 0.0), (12.0, 0.0)), true, false, SegmentAnnotations(true, false)),
                SegmentEvent(((12.0, 0.0), (12.0, 3.0)), true, false, SegmentAnnotations(true, false)),
                SegmentEvent(((-1.0, 3.0), (12.0, 3.0)), true, false, SegmentAnnotations(false, true)),
            ]
            queue = SegmentEvent{Float64}[]
            for ev in vcat(rectangle, rectangle_horiz)
                add_annotated_segment!(queue, ev)
            end
            annotated_segments = event_loop!(queue; self_intersection=false)
            expected = [
                SegmentEvent(((-1.0, 0.0), (-1.0, 3.0)), true, false, SegmentAnnotations(false, true), SegmentAnnotations(false, false)),
                SegmentEvent(((-1.0, 0.0), (0.0, 0.0)), true, false, SegmentAnnotations(true, false), SegmentAnnotations(false, false)),
                SegmentEvent(((0.0, 0.0), (3.0, 3.0)), true, true, SegmentAnnotations(false, true), SegmentAnnotations(true, true)),
                SegmentEvent(((-1.0, 3.0), (3.0, 3.0)), true, false, SegmentAnnotations(false, true), SegmentAnnotations(false, false)),
                SegmentEvent(((0.0, 0.0), (4.0, -4.0)), true, true, SegmentAnnotations(true, false), SegmentAnnotations(false, false)),
                SegmentEvent(((0.0, 0.0), (6.0, 0.0)), true, false, SegmentAnnotations(true, false), SegmentAnnotations(true, true)),
                SegmentEvent(((3.0, 3.0), (6.0, 0.0)), true, true, SegmentAnnotations(false, true), SegmentAnnotations(true, true)),
                SegmentEvent(((4.0, -4.0), (7.0, -1.0)), true, true, SegmentAnnotations(true, false), SegmentAnnotations(false, false)),
                SegmentEvent(((6.0, 0.0), (7.0, -1.0)), true, true, SegmentAnnotations(false, true), SegmentAnnotations(false, false)),
                SegmentEvent(((6.0, 0.0), (12.0, 0.0)), true, false, SegmentAnnotations(true, false), SegmentAnnotations(false, false)),
                SegmentEvent(((12.0, 0.0), (12.0, 3.0)), true, false, SegmentAnnotations(true, false), SegmentAnnotations(false, false)),
                SegmentEvent(((3.0, 3.0), (12.0, 3.0)), true, false, SegmentAnnotations(false, true), SegmentAnnotations(false, false)),
            ]
            @test annotated_segments == expected
        end

        @testset "rectangle self-intersect annotations" begin
            self_intersect = [
                SegmentEvent(((0.0, 0.0), (2.0, 2.0)), true, true, SegmentAnnotations(false, true)),
                SegmentEvent(((0.0, 0.0), (4.0, 0.0)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((2.0, 2.0), (4.0, 0.0)), true, true, SegmentAnnotations(false, true)),
                SegmentEvent(((4.0, 0.0), (6.0, -2.0)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((6.0, -2.0), (8.5, -0.0)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((4.0, 0.0), (8.5, -0.0)), true, true, SegmentAnnotations(false, true)),
                SegmentEvent(((8.5, -0.0), (11.0, 0.0)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((11.0, 0.0), (11.0, 2.0)), true, true, SegmentAnnotations(true, false)),
                SegmentEvent(((8.5, -0.0), (11.0, 2.0)), true, true, SegmentAnnotations(false, true)),
            ]
            rectangle_horiz = [
                SegmentEvent(((-1.0, 0.0), (-1.0, 3.0)), true, false, SegmentAnnotations(false, true)),
                SegmentEvent(((-1.0, 0.0), (12.0, 0.0)), true, false, SegmentAnnotations(true, false)),
                SegmentEvent(((12.0, 0.0), (12.0, 3.0)), true, false, SegmentAnnotations(true, false)),
                SegmentEvent(((-1.0, 3.0), (12.0, 3.0)), true, false, SegmentAnnotations(false, true)),
            ]
            queue = SegmentEvent{Float64}[]
            for ev in vcat(self_intersect, rectangle_horiz)
                add_annotated_segment!(queue, ev)
            end
            annotated_segments = event_loop!(queue; self_intersection=false)
            expected = [
                SegmentEvent(((-1.0, 0.0), (-1.0, 3.0)), true, false, SegmentAnnotations(false, true), SegmentAnnotations(false, false)),
                SegmentEvent(((-1.0, 0.0), (0.0, 0.0)), true, false, SegmentAnnotations(true, false), SegmentAnnotations(false, false)),
                SegmentEvent(((0.0, 0.0), (2.0, 2.0)), true, true, SegmentAnnotations(false, true), SegmentAnnotations(true, true)),
                SegmentEvent(((0.0, 0.0), (4.0, 0.0)), true, false, SegmentAnnotations(true, false), SegmentAnnotations(true, false)),
                SegmentEvent(((2.0, 2.0), (4.0, 0.0)), true, true, SegmentAnnotations(false, true), SegmentAnnotations(true, true)),
                SegmentEvent(((4.0, 0.0), (6.0, -2.0)), true, true, SegmentAnnotations(true, false), SegmentAnnotations(false, false)),
                SegmentEvent(((6.0, -2.0), (8.5, -0.0)), true, true, SegmentAnnotations(true, false), SegmentAnnotations(false, false)),
                SegmentEvent(((4.0, 0.0), (8.5, -0.0)), true, false, SegmentAnnotations(true, false), SegmentAnnotations(false, true)),
                SegmentEvent(((8.5, -0.0), (11.0, 0.0)), true, false, SegmentAnnotations(true, false), SegmentAnnotations(true, false)),
                SegmentEvent(((11.0, 0.0), (11.0, 2.0)), true, true, SegmentAnnotations(true, false), SegmentAnnotations(true, true)),
                SegmentEvent(((8.5, -0.0), (11.0, 2.0)), true, true, SegmentAnnotations(false, true), SegmentAnnotations(true, true)),
                SegmentEvent(((11.0, 0.0), (12.0, 0.0)), true, false, SegmentAnnotations(true, false), SegmentAnnotations(false, false)),
                SegmentEvent(((12.0, 0.0), (12.0, 3.0)), true, false, SegmentAnnotations(true, false), SegmentAnnotations(false, false)),
                SegmentEvent(((-1.0, 3.0), (12.0, 3.0)), true, false, SegmentAnnotations(false, true), SegmentAnnotations(false, false)),
            ]
            @test annotated_segments == expected
        end
    end

    @testset "selection criteria" begin
        function check_selection_criteria(calc_criteria, criteria::Vector)
            for idx in 0:15
                above1 = idx > 7
                below1 = (idx % 8) > 3
                above2 = (idx % 4) > 1
                below2 = (idx % 2) > 0
                expect_filled = calc_criteria(above1, below1, above2, below2)
                #println(idx + 1, " $above1 $below1 $above2 $below2 - ", expect_filled)
                @assert expect_filled ? criteria[idx + 1] != BLANK : criteria[idx + 1] == BLANK
            end
            true
        end

        @testset "intersection" begin
            function is_intersect(above1, below1, above2, below2)
                is_intersect = ((above1 & above2) ⊻ (below1 & below2))
                is_intersect_segments = (above1 & below2) ⊻ (below1 & above2)
                (is_intersect || is_intersect_segments)
            end

            @test check_selection_criteria(
                is_intersect,
                PolygonAlgorithms.INTERSECTION_CRITERIA
            )
        end
        @testset "union" begin
            @test check_selection_criteria(
                (above1, below1, above2, below2) -> (above1 | above2) ⊻ (below1 | below2),
                PolygonAlgorithms.UNION_CRITERIA
            )
        end

        @testset "difference" begin
            @test check_selection_criteria(
                (above1, below1, above2, below2) -> (above1 & !above2) ⊻ (below1 & !below2),
                PolygonAlgorithms.DIFFERENCE_CRITERIA
            )
        end

        @testset "xor" begin
            @test check_selection_criteria(
                (above1, below1, above2, below2) -> (above1 ⊻ above2) ⊻ (below1 ⊻ below2),
                PolygonAlgorithms.XOR_CRITERIA
            )
        end
    end

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
    end

    @testset "holes" begin
        @testset "two squares with holes" begin
            regions = [
                [
                    SegmentEvent(((2.0, 1.0), (2.0, 2.0)), false, true, SegmentAnnotations(false, true)),
                    SegmentEvent(((1.0, 1.0), (2.0, 1.0)), false, true, SegmentAnnotations(false, true)),
                    SegmentEvent(((1.0, 1.0), (1.0, 2.0)), true,  true, SegmentAnnotations(true, false)),
                    SegmentEvent(((1.0, 1.0), (1.0, 2.0)), false, true, SegmentAnnotations(true, false)),
                ],
                [
                    SegmentEvent(((3.0, 0.0), (3.0, 3.0)), false, true, SegmentAnnotations(true, false)),
                    SegmentEvent(((0.0, 0.0), (3.0, 0.0)), false, true, SegmentAnnotations(true, false)),
                    SegmentEvent(((0.0, 0.0), (0.0, 3.0)), true,  true, SegmentAnnotations(false, true)),
                    SegmentEvent(((0.0, 0.0), (0.0, 3.0)), false, true, SegmentAnnotations(false, true)),
                ],
                [
                    SegmentEvent(((6.0, 1.0), (6.0, 2.0)), false, true, SegmentAnnotations(false, true)),
                    SegmentEvent(((5.0, 1.0), (6.0, 1.0)), false, true, SegmentAnnotations(false, true)),
                    SegmentEvent(((5.0, 1.0), (5.0, 2.0)), true,  true, SegmentAnnotations(true, false)),
                    SegmentEvent(((5.0, 1.0), (5.0, 2.0)), false, true, SegmentAnnotations(true, false)),
                ],
                [
                   SegmentEvent(((7.0, 0.0), (7.0, 3.0)), false, true, SegmentAnnotations(true, false))
                    SegmentEvent(((4.0, 0.0), (7.0, 0.0)), false, true, SegmentAnnotations(true, false))
                    SegmentEvent(((4.0, 0.0), (4.0, 3.0)), true,  true, SegmentAnnotations(false, true))
                    SegmentEvent(((4.0, 0.0), (4.0, 3.0)), false, true, SegmentAnnotations(false, true))
                ]
            ]
            classification = map(is_hole, regions)
            expected = [true, false, true, false]
            @test classification == expected
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
end
