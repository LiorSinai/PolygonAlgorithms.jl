using Test
using PolygonAlgorithms: SegmentEvent, SegmentAnnotations
using PolygonAlgorithms: compare_events, convert_to_event_queue
using PolygonAlgorithms: find_transition, any_intersect

@testset "line Sweep" verbose=false begin
    @testset "queueing" begin
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

            lowered_triangle = [
                (3.0, 3.0), (3.0, -4.0), (0.0, -5.0)
            ]
            event_queue = convert_to_event_queue(lowered_triangle)
            expected = [
                SegmentEvent(((0.0, -5.0), (3.0, -4.0)), true),
                SegmentEvent(((0.0, -5.0), (3.0, 3.0)), true),
                SegmentEvent(((0.0, -5.0), (3.0, -4.0)), false),
                SegmentEvent(((3.0, -4.0), (3.0, 3.0)), true),
                SegmentEvent(((3.0, -4.0), (3.0, 3.0)), false),
                SegmentEvent(((0.0, -5.0), (3.0, 3.0)), false),        
            ]
            @test event_queue == expected

            mirror_triangle = [
                (0.0, -4.0), (0.0, 3.0), (3.0, -5.0)
            ]
            event_queue = convert_to_event_queue(mirror_triangle)
            expected = [
                SegmentEvent(((0.0, -4.0), (3.0, -5.0)), true),   
                SegmentEvent(((0.0, -4.0), (0.0, 3.0)), true),    
                SegmentEvent(((0.0, -4.0), (0.0, 3.0)), false),   
                SegmentEvent(((0.0, 3.0), (3.0, -5.0)), true),
                SegmentEvent(((0.0, -4.0), (3.0, -5.0)), false),  
                SegmentEvent(((0.0, 3.0), (3.0, -5.0)), false),        
            ]
            @test event_queue == expected
        end
    end

    @testset "find_transitions" begin
        @testset "horizontal lines" begin
            # setup
            status = [
                SegmentEvent(((0.0, 5.0), (5.0, 5.0)), true),
                SegmentEvent(((0.0, 3.0), (5.0, 3.0)), true),
                SegmentEvent(((0.0, 1.0), (5.0, 1.0)), true),
            ]
            # test
            ev6 = SegmentEvent(((0.0, 6.0), (5.0, 6.0)), true)
            idx = find_transition(status, ev6)
            @test idx == 1
            ev34 = SegmentEvent(((0.0, 3.0), (5.0, 4.0)), true)
            idx = find_transition(status, ev34)
            @test idx == 2
            ev32 = SegmentEvent(((0.0, 3.0), (5.0, 2.0)), true)
            idx = find_transition(status, ev32)
            @test idx == 3
            ev2 = SegmentEvent(((0.0, 2.0), (5.0, 2.0)), true)
            idx = find_transition(status, ev2)
            @test idx == 3
            ev0 = SegmentEvent(((0.0, 0.0), (5.0, 0.0)), true)
            idx = find_transition(status, ev0)
            @test idx == 4
        end

        @testset "sweep status head" begin
            # setup
            status = [SegmentEvent(((0.0, 5.0), (5.0, 5.0)), true)]
            # test
            ev6 = SegmentEvent(((0.0, 6.0), (5.0, 6.0)), true)
            idx = find_transition(status, ev6)
            @test idx == 1
            ev2 = SegmentEvent(((0.0, 2.0), (5.0, 2.0)), true)
            idx = find_transition(status, ev2)
            @test idx == 2
        end

        @testset "size discrepancy" begin
            # from spiral-star example. 
            status = [
                SegmentEvent(((-10.0, -12.0), (0.0, -2.0)), true) # star leg
            ]
            # test
            ev = SegmentEvent(((-8.338, -10.337), (-8.239, -10.426)), true) # tiny segment of spiral
            idx = find_transition(status, ev)
            @test idx == 1
        end
    end

    @testset "order sweep status" begin
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
                idx = find_transition(sweep_status, event)
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
                SegmentEvent(((0.4767, 0.5368), (0.6615, 0.525)), true),
                SegmentEvent(((0.5824, 0.5157), (0.9235, 0.3918)), true),
                SegmentEvent(((0.5824, 0.5157), (0.9756, 0.3288)), true),
                SegmentEvent(((0.4127, 0.4774), (0.6615, 0.525)), true),
                SegmentEvent(((0.4771, 0.4543), (0.5913, 0.411)), true),
                SegmentEvent(((0.5218, 0.2095), (0.5839, 0.1114)), true),
                SegmentEvent(((0.4544, 0.0837), (0.6762, 0.4073)), true),
                SegmentEvent(((0.4544, 0.0837), (0.7364, 0.3326)), true),
            ]
            @test sweep_status == expected
        end

        @testset "sweep status rand lines 2" begin
            segments = [
                SegmentEvent(((0.01190, 0.9843), (0.2696, 0.7221)), true),
                SegmentEvent(((0.02439, 0.6210), (0.1274, 0.6140)), true),
                SegmentEvent(((0.01190, 0.9843), (0.1166, 0.9564)), true),
                SegmentEvent(((0.02706, 0.3861), (0.0599, 0.4676)), true),
                SegmentEvent(((0.00284, 0.3700), (0.0727, 0.1129)), true),
                SegmentEvent(((0.02201, 0.5214), (0.3929, 0.5001)), true),
            ]
            sweep_status = SegmentEvent{Float64}[]
            for event in segments
                idx = find_transition(sweep_status, event)
                insert!(sweep_status, idx, event)
            end
            expected =[
                SegmentEvent(((0.0119, 0.9843), (0.1166, 0.9564)), true),
                SegmentEvent(((0.0119, 0.9843), (0.2696, 0.7221)), true),
                SegmentEvent(((0.02439, 0.621), (0.1274, 0.614)), true),
                SegmentEvent(((0.02201, 0.5214), (0.3929, 0.5001)), true),
                SegmentEvent(((0.02706, 0.3861), (0.0599, 0.4676)), true),
                SegmentEvent(((0.00284, 0.37), (0.0727, 0.1129)), true),
            ]
            @test sweep_status == expected
        end

        @testset "sweep status almost colinear" begin
            sweep_status =[
                SegmentEvent(((0.7481773880710287, 0.2525213105681051), (0.891378804016348, 0.0914800094089202)), true),
                SegmentEvent(((0.7674956225378292, 0.19603318345420395), (0.93007516434479, 0.04850368886951295)), true),
                SegmentEvent(((0.7307442328548945, 0.22237367849003833), (0.879152990839225, 0.03969376626137666)), true),
                SegmentEvent(((0.7456555133038324, 0.17560361123176102), (0.879152990839225, 0.03969376626137666)), true)
            ]
            ev = SegmentEvent(((0.7692174020749689, 0.1749860342238233), (0.8076322984621129, 0.25945610393589347)), true)
            # with rtol=1e-3, this is linear
            idx = find_transition(sweep_status, ev)
            @test idx == 4
        end

        @testset "sweep status steps" begin
            # lines originally from the intersection of two randomly generated spiky polygon
            segments = [
                SegmentEvent(((0.0, 4.0), (3.0, 4.0)), true),
                SegmentEvent(((0.0, 2.0), (3.0, 2.0)), true),
                SegmentEvent(((3.0, 2.0), (3.0, -2.0)), true),
                SegmentEvent(((3.0, -2.0), (6.0, -2.0)), true),
                SegmentEvent(((3.0, -4.0), (6.0, -4.0)), true),
            ]
            sweep_status = SegmentEvent{Float64}[]
            for event in segments
                idx = find_transition(sweep_status, event)
                insert!(sweep_status, idx, event)
            end
            expected = [
                SegmentEvent(((0.0, 4.0), (3.0, 4.0)), true),
                SegmentEvent(((0.0, 2.0), (3.0, 2.0)), true),
                SegmentEvent(((3.0, 2.0), (3.0, -2.0)), true),
                SegmentEvent(((3.0, -2.0), (6.0, -2.0)), true),
                SegmentEvent(((3.0, -4.0), (6.0, -4.0)), true),
            ]
            @test sweep_status == expected
        end
    end

    @testset "any intersect" begin
        @testset "2 segments" begin
            segments = [
                ((1.0, 1.0), (2.0, 3.0)), ((2.0, 2.0), (5.0, 0.0))
            ]
            @test !any_intersect(segments...)
            # touch on line
            segments = [
                ((1.0, 1.0), (3.0, 3.0)), ((2.0, 2.0), (5.0, 0.0))
            ]
            @test any_intersect(segments...)
            # intersect
            segments = [
                ((1.0, 1.0), (4.0, 3.0)), ((2.0, 2.0), (5.0, 0.0))
            ]
            @test any_intersect(segments...)
            # same origin
            segments = [
                ((1.0, 1.0), (2.0, 3.0)), ((1.0, 1.0), (5.0, 0.0))
            ]
            @test any_intersect(segments...)
            @test !any_intersect(segments...; exclude_connected=true)
        end

        @testset "random lines" begin 
            segments = [
                ((0.162, 0.019), (0.214, 0.350)),
                ((0.643, 0.899), (0.813, 0.840)),
                ((0.060, 0.509), (0.253, 0.895)),
                ((0.296, 0.552), (0.300, 0.691)),
                ((0.353, 0.664), (0.364, 0.605)),
            ]
            @test !any_intersect(segments...)
            segments = [
                ((0.549, 0.831), (0.743, 0.401)),
                ((0.271, 0.306), (0.590, 0.558)),
                ((0.610, 0.515), (0.895, 0.327)),
                ((0.017, 0.434), (0.476, 0.946)),
                ((0.385, 0.793), (0.461, 0.063)),
            ]
            @test any_intersect(segments...)
            # end event intersects
            segments = [
                ((0.05, 0.106), (0.375, 0.307)),
                ((0.13, 0.780), (0.610, 0.989)),
                ((0.01, 0.252), (0.782, 0.591)),
                ((0.60, 0.876), (0.809, 0.374)),
                ((0.21, 0.617), (0.428, 0.544)),
            ]
            @test any_intersect(segments...)
        end

        @testset "parallel lines" begin
            # horizontal
            segments = [
                ((0.0, 0.2), (1.0, 0.2)),
                ((0.0, 0.4), (1.0, 0.4)),
                ((0.0, 0.6), (1.0, 0.6)),
                ((0.0, 0.8), (1.0, 0.8)),
                ((0.0, 1.0), (1.0, 1.0)),
            ]
            @test !any_intersect(segments...)
            push!(segments, ((0.8, 0.3), (0.9, 0.5)))
            @test any_intersect(segments...)
            # vertical
            segments = [
                ((0.2, 0.0), (0.2, 1.0)),
                ((0.4, 0.0), (0.4, 1.0)),
                ((0.6, 0.0), (0.6, 1.0)),
                ((0.8, 0.0), (0.8, 1.0)),
                ((1.0, 0.0), (1.0, 1.0)),
            ]
            @test !any_intersect(segments...)
            push!(segments, ((0.75, 0.1), (0.85, 0.2)));
            @test any_intersect(segments...)
        end

        @testset "connected" begin
            segments = [
                ((3.0, 3.0), (7.0, -1.0)),
                ((4.0, -4.0), (7.0, -1.0)),
                ((0.0, 0.0), (4.0, -4.0)),
                ((0.0, 0.0), (3.0, 3.0)),
            ]
            @test any_intersect(segments...)
            @test !any_intersect(segments...; exclude_connected=true)
        end
    end
end
