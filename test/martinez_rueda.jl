using Test
using PolygonAlgorithms: SegmentEvent, SegmentAnnotations
using PolygonAlgorithms: add_annotated_segment!, convert_to_event_queue, check_and_divide_intersection!, event_loop!, find_transition
using PolygonAlgorithms: chain_segments

@testset "Martinez-Rueda algorithm" begin
    @testset "rectangle" begin
        rectangle = [
            (3.0, 3.0), (7.0, -1.0), (4.0, -4.0), (0.0, 0.0)
        ]
        event_queue = convert_to_event_queue(rectangle)
        expected = [
            SegmentEvent(((0.0, 0.0), (4.0, -4.0)), true),
            SegmentEvent(((0.0, 0.0), (3.0, 3.0)), true),
            SegmentEvent(((0.0, 0.0), (3.0, 3.0)), false),
            SegmentEvent(((3.0, 3.0), (7.0, -1.0)), true),
            SegmentEvent(((0.0, 0.0), (4.0, -4.0)), false),
            SegmentEvent(((4.0, -4.0), (7.0, -1.0)), true),
            SegmentEvent(((4.0, -4.0), (7.0, -1.0)), false),
            SegmentEvent(((3.0, 3.0), (7.0, -1.0)), false),
        ]
        @test event_queue == expected
    end

    @testset "vertical line" begin
        triangle = [
            (3.0, 3.0), (3.0, -4.0), (0.0, 0.0)
        ]
        event_queue = convert_to_event_queue(triangle)
        expected = [
            SegmentEvent(((0.0, 0.0), (3.0, -4.0)), true),
            SegmentEvent(((0.0, 0.0), (3.0, 3.0)), true),
            SegmentEvent(((0.0, 0.0), (3.0, -4.0)), false),
            SegmentEvent(((3.0, -4.0), (3.0, 3.0)), true),
            SegmentEvent(((3.0, -4.0), (3.0, 3.0)), false),
            SegmentEvent(((0.0, 0.0), (3.0, 3.0)), false),        
        ]
        @test event_queue == expected
    end

    @testset "find_transitions" begin
        # setup
        status = [
            SegmentEvent(((0.0, 5.0), (5.0, 5.0)), true),
            SegmentEvent(((0.0, 3.0), (5.0, 3.0)), true),
            SegmentEvent(((0.0, 1.0), (5.0, 1.0)), true),
        ]
        # test
        ev6 = SegmentEvent(((0.0, 6.0), (5.0, 6.0)), true)
        idx, prev, next = find_transition(status, ev6)
        @test isnothing(prev) && next == status[1]
        ev34 = SegmentEvent(((0.0, 3.0), (5.0, 4.0)), true)
        idx, prev, next = find_transition(status, ev34)
        @test prev == status[1] && next == status[2]
        ev32 = SegmentEvent(((0.0, 3.0), (5.0, 2.0)), true)
        idx, prev, next = find_transition(status, ev32)
        @test prev == status[2] && next == status[3]
        ev2 = SegmentEvent(((0.0, 2.0), (5.0, 2.0)), true)
        idx, prev, next = find_transition(status, ev2)
        @test prev == status[2] && next == status[3]
        ev0 = SegmentEvent(((0.0, 0.0), (5.0, 0.0)), true)
        idx, prev, next = find_transition(status, ev0)
        @test prev == status[3] && isnothing(next)
    end

    @testset "sweep status head" begin
        # setup
        status = [SegmentEvent(((0.0, 5.0), (5.0, 5.0)), true)]
        # test
        ev6 = SegmentEvent(((0.0, 6.0), (5.0, 6.0)), true)
        idx, prev, next = find_transition(status, ev6)
        @test isnothing(prev) && next == status[1]
        ev2 = SegmentEvent(((0.0, 2.0), (5.0, 2.0)), true)
        idx, prev, next = find_transition(status, ev2)
        @test prev == status[1] && isnothing(next)
    end

    @testset "sweep status rand lines" begin
        # lines originally from the intersection of two randomly generated spiky polygon
        segments = [
            SegmentEvent(((0.5824, 0.5157), (0.9235, 0.39180)), true)
            SegmentEvent(((0.5824, 0.5157), (0.9756, 0.3288)), true)
            SegmentEvent(((0.5732, 0.7060), (0.7169, 0.6187)), true)
            SegmentEvent(((0.5658, 0.6475), (0.6321, 0.6197)), true)
            SegmentEvent(((0.5658, 0.6475), (0.7169, 0.6187)), true)
            SegmentEvent(((0.5480, 0.9443), (0.7438, 0.6878)), true)
            SegmentEvent(((0.5480, 0.9443), (0.7455, 0.8477)), true)
            SegmentEvent(((0.5334, 0.9450), (0.7438, 0.6878)), true)
            SegmentEvent(((0.5218, 0.2095), (0.5839, 0.1114)), true)
            SegmentEvent(((0.5122, 0.9013), (0.6366, 0.9417)), true)
            SegmentEvent(((0.4867, 0.5998), (0.6095, 0.6079)), true)
            SegmentEvent(((0.4771, 0.4543), (0.5913, 0.4110)), true)
            SegmentEvent(((0.4767, 0.5368), (0.6615, 0.5250)), true)
            SegmentEvent(((0.4544, 0.0837), (0.6762, 0.4073)), true)
            SegmentEvent(((0.4544, 0.0837), (0.7364, 0.3326)), true)
            SegmentEvent(((0.4127, 0.4774), (0.6615, 0.5250)), true)
        ]
        sweep_status = SegmentEvent{Float64}[]
        for event in segments
            idx, above, below = find_transition(sweep_status, event)
            insert!(sweep_status, idx, event)
        end
        expected = [
            SegmentEvent(((0.5480, 0.9443), (0.7455, 0.8477)), true),
            SegmentEvent(((0.5480, 0.9443), (0.7438, 0.6878)), true),
            SegmentEvent(((0.5334, 0.9450), (0.7438, 0.6878)), true),
            SegmentEvent(((0.5122, 0.9013), (0.6366, 0.9417)), true),
            SegmentEvent(((0.5732, 0.7060), (0.7169, 0.6187)), true),
            SegmentEvent(((0.5658, 0.6475), (0.7169, 0.6187)), true),
            SegmentEvent(((0.5658, 0.6475), (0.6321, 0.6197)), true),
            SegmentEvent(((0.4867, 0.5998), (0.6095, 0.6079)), true),
            SegmentEvent(((0.5824, 0.5157), (0.9235, 0.3918)), true),
            SegmentEvent(((0.5824, 0.5157), (0.9756, 0.3288)), true),
            SegmentEvent(((0.4767, 0.5368), (0.6615, 0.525)), true),
            SegmentEvent(((0.4771, 0.4543), (0.5913, 0.411)), true),
            SegmentEvent(((0.4127, 0.4774), (0.6615, 0.525)), true),
            SegmentEvent(((0.5218, 0.2095), (0.5839, 0.1114)), true),
            SegmentEvent(((0.4544, 0.0837), (0.6762, 0.4073)), true),
            SegmentEvent(((0.4544, 0.0837), (0.7364, 0.3326)), true),
        ]
        @test sweep_status == expected
    end

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
            SegmentEvent(((-0.3062015503875969, -0.3798449612403101), (-0.18750000000000003, 0.0)), true, true, SegmentAnnotations(true, false)),
            SegmentEvent(((-0.8, 0.0), (-0.18750000000000003, 0.0)), true, true, SegmentAnnotations(false, true)),
            SegmentEvent(((-0.5, -1.0), (-0.0, -0.6153846153846154)), true, true, SegmentAnnotations(true, false)),
            SegmentEvent(((-0.3062015503875969, -0.3798449612403101), (-0.0, -0.6153846153846154)), true, true, SegmentAnnotations(false, true)),
            SegmentEvent(((-0.18750000000000003, 0.0), (0.0, 0.6)), true, true, SegmentAnnotations(false, true)),
            SegmentEvent(((-0.18750000000000003, 0.0), (0.1875, 0.0)), true, true, SegmentAnnotations(true, false)),
            SegmentEvent(((0.0, 0.6), (0.1875, 0.0)), true, true, SegmentAnnotations(false, true)),
            SegmentEvent(((-0.0, -0.6153846153846154), (0.30620155038759694, -0.3798449612403101)), true, true, SegmentAnnotations(false, true)),
            SegmentEvent(((0.1875, 0.0), (0.30620155038759694, -0.3798449612403101)), true, true, SegmentAnnotations(true, false)),
            SegmentEvent(((-0.0, -0.6153846153846154), (0.5, -1.0)), true, true, SegmentAnnotations(true, false)),
            SegmentEvent(((0.30620155038759694, -0.3798449612403101), (0.5, -1.0)), true, true, SegmentAnnotations(false, true)),
            SegmentEvent(((0.30620155038759694, -0.3798449612403101), (0.8, 0.0)), true, true, SegmentAnnotations(true, false)),
            SegmentEvent(((0.1875, 0.0), (0.8, 0.0)), true, true, SegmentAnnotations(false, true)),
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

    @testset "chain segments - rectangle" begin
        segments = [
            SegmentEvent(((0.0, 0.0), (4.0, -4.0)), true),
            SegmentEvent(((0.0, 0.0), (3.0, 3.0)), true),
            SegmentEvent(((3.0, 3.0), (7.0, -1.0)), true),
            SegmentEvent(((4.0, -4.0), (7.0, -1.0)), true),
        ]
        regions = chain_segments(segments)
        expected = [[(7.0, -1.0), (3.0, 3.0), (0.0, 0.0), (4.0, -4.0)]]
        @test regions == expected
    end

    @testset "chain segments - improper" begin
        segments = [
            SegmentEvent(((0.0, 0.0), (2.0, 2.0)), true),
            SegmentEvent(((0.0, 0.0), (2.0, 2.0)), true),
            SegmentEvent(((2.0, 2.0), (5.0, 1.0)), true),
            SegmentEvent(((2.0, 2.0), (3.0, 5.0)), true),
            SegmentEvent(((3.0, 5.0), (5.0, 1.0)), true),
        ]
        expected = [(0.0, 0.0), (2.0, 2.0), (3.0, 5.0), (5.0, 1.0), (2.0, 2.0)]
        @test_throws AssertionError chain_segments(segments)
    end
end
