using PolygonAlgorithms: directed_graph_from_segments, map_connections, compute_graph_faces
using PolygonAlgorithms: SegmentEvent

@testset "graph face computation segments" begin
    diamond = [
        SegmentEvent(((1.0, 1.0), (3.0, 4.0)), true),
        SegmentEvent(((3.0, 4.0), (5.0, 1.0)), true),
        SegmentEvent(((5.0, 1.0), (3.0, -2.0)), true),
        SegmentEvent(((3.0, -2.0), (1.0, 1.0)), true),
    ];
    cross = [
        SegmentEvent(((1.0, 1.0), (3.0, 1.0)), true),
        SegmentEvent(((3.0, 1.0), (6.0, 1.0)), false),
        SegmentEvent(((3.0, 1.0), (3.0, 4.0)), true),
        SegmentEvent(((3.0, 1.0), (3.0, -2.0)), true),
    ];
    improper = [
        SegmentEvent(((1.0, 1.0), (3.0, 4.0)), true),
        SegmentEvent(((3.0, 4.0), (5.0, 1.0)), true),
        SegmentEvent(((5.0, 1.0), (7.0, 2.0)), true), # connection
        SegmentEvent(((7.0, 2.0), (8.0, 3.0)), true), # tail
        SegmentEvent(((8.0, 3.0), (8.0, 1.0)), true),
        SegmentEvent(((8.0, 1.0), (7.0, 2.0)), true), # close triangle
        SegmentEvent(((5.0, 1.0), (3.0, -2.0)), true),
        SegmentEvent(((3.0, -2.0), (1.0, 1.0)), true),
    ];

    @testset "cross from segments" begin
        graph = directed_graph_from_segments(cross)
        expected = Dict(
            (3.0, 1.0)  => [
                SegmentEvent(((3.0, 1.0), (1.0, 1.0))),
                SegmentEvent(((3.0, 1.0), (3.0, 4.0))),
                SegmentEvent(((3.0, 1.0), (6.0, 1.0))),
                SegmentEvent(((3.0, 1.0), (3.0, -2.0))),
            ],
            (1.0, 1.0)  => [SegmentEvent(((1.0, 1.0), (3.0, 1.0)))],
            (6.0, 1.0)  => [SegmentEvent(((6.0, 1.0), (3.0, 1.0)))],
            (3.0, 4.0)  => [SegmentEvent(((3.0, 4.0), (3.0, 1.0)))],
            (3.0, -2.0) => [SegmentEvent(((3.0, -2.0), (3.0, 1.0)))],
        )
        @test graph == expected
    end

    @testset "diamond from segments" begin
        graph = directed_graph_from_segments(diamond)
        expected = Dict(
            (1.0, 1.0)  => [SegmentEvent(((1.0, 1.0), (3.0, 4.0))), SegmentEvent(((1.0, 1.0), (3.0, -2.0)))],
            (5.0, 1.0)  => [SegmentEvent(((5.0, 1.0), (3.0, 4.0))), SegmentEvent(((5.0, 1.0), (3.0, -2.0)))],
            (3.0, 4.0)  => [SegmentEvent(((3.0, 4.0), (5.0, 1.0))), SegmentEvent(((3.0, 4.0), (1.0, 1.0)))],
            (3.0, -2.0) => [SegmentEvent(((3.0, -2.0), (1.0, 1.0))), SegmentEvent(((3.0, -2.0), (5.0, 1.0)))],
        )
        @test graph == expected
    end

    @testset "improper from segments" begin
        graph = directed_graph_from_segments(improper)
        expected = Dict(
            (8.0, 1.0) => [SegmentEvent(((8.0, 1.0), (7.0, 2.0))), SegmentEvent(((8.0, 1.0), (8.0, 3.0)))],
            (1.0, 1.0) => [SegmentEvent(((1.0, 1.0), (3.0, 4.0))), SegmentEvent(((1.0, 1.0), (3.0, -2.0)))],
            (8.0, 3.0) => [SegmentEvent(((8.0, 3.0), (8.0, 1.0))), SegmentEvent(((8.0, 3.0), (7.0, 2.0)))],
            (5.0, 1.0) => [SegmentEvent(((5.0, 1.0), (3.0, 4.0))), SegmentEvent(((5.0, 1.0), (7.0, 2.0))), SegmentEvent(((5.0, 1.0), (3.0, -2.0)))],
            (7.0, 2.0) => [SegmentEvent(((7.0, 2.0), (8.0, 3.0))), SegmentEvent(((7.0, 2.0), (8.0, 1.0))), SegmentEvent(((7.0, 2.0), (5.0, 1.0)))],
            (3.0, 4.0) => [SegmentEvent(((3.0, 4.0), (5.0, 1.0))), SegmentEvent(((3.0, 4.0), (1.0, 1.0)))],
            (3.0, -2.0) => [SegmentEvent(((3.0, -2.0), (1.0, 1.0))), SegmentEvent(((3.0, -2.0), (5.0, 1.0))),],
        )
        @test graph == expected
    end

    @testset "connections - cross" begin
        graph = Dict(
            (3.0, 1.0)  => [
                SegmentEvent(((3.0, 1.0), (1.0, 1.0))),
                SegmentEvent(((3.0, 1.0), (3.0, 4.0))),
                SegmentEvent(((3.0, 1.0), (6.0, 1.0))),
                SegmentEvent(((3.0, 1.0), (3.0, -2.0))),
            ],
            (1.0, 1.0)  => [SegmentEvent(((1.0, 1.0), (3.0, 1.0)))],
            (6.0, 1.0)  => [SegmentEvent(((6.0, 1.0), (3.0, 1.0)))],
            (3.0, 4.0)  => [SegmentEvent(((3.0, 4.0), (3.0, 1.0)))],
            (3.0, -2.0) => [SegmentEvent(((3.0, -2.0), (3.0, 1.0)))],
        )
        connections = map_connections(graph)
        expected = Dict(
            ((6.0, 1.0), (3.0, 1.0))  => SegmentEvent(((3.0, 1.0), (3.0, 4.0)), true),
            ((3.0, 4.0), (3.0, 1.0))  => SegmentEvent(((3.0, 1.0), (1.0, 1.0)), true),
            ((3.0, 1.0), (1.0, 1.0))  => SegmentEvent(((1.0, 1.0), (3.0, 1.0)), true),
            ((3.0, -2.0), (3.0, 1.0)) => SegmentEvent(((3.0, 1.0), (6.0, 1.0)), true),
            ((3.0, 1.0), (6.0, 1.0))  => SegmentEvent(((6.0, 1.0), (3.0, 1.0)), true),
            ((3.0, 1.0), (3.0, -2.0)) => SegmentEvent(((3.0, -2.0), (3.0, 1.0)), true),
            ((3.0, 1.0), (3.0, 4.0))  => SegmentEvent(((3.0, 4.0), (3.0, 1.0)), true),
            ((1.0, 1.0), (3.0, 1.0))  => SegmentEvent(((3.0, 1.0), (3.0, -2.0)), true),
        )
        @test connections == expected
    end

    @testset "connections - diamond" begin
        graph = Dict(
            (1.0, 1.0)  => [SegmentEvent(((1.0, 1.0), (3.0, 4.0))), SegmentEvent(((1.0, 1.0), (3.0, -2.0)))],
            (5.0, 1.0)  => [SegmentEvent(((5.0, 1.0), (3.0, 4.0))), SegmentEvent(((5.0, 1.0), (3.0, -2.0)))],
            (3.0, 4.0)  => [SegmentEvent(((3.0, 4.0), (5.0, 1.0))), SegmentEvent(((3.0, 4.0), (1.0, 1.0)))],
            (3.0, -2.0) => [SegmentEvent(((3.0, -2.0), (1.0, 1.0))), SegmentEvent(((3.0, -2.0), (5.0, 1.0)))],
        )
        connections = map_connections(graph)
        expected = Dict(
            ((3.0, 4.0), (5.0, 1.0))  => SegmentEvent(((5.0, 1.0), (3.0, -2.0)), true),
            ((5.0, 1.0), (3.0, 4.0))  => SegmentEvent(((3.0, 4.0), (1.0, 1.0)), true,),
            ((3.0, -2.0), (1.0, 1.0)) => SegmentEvent(((1.0, 1.0), (3.0, 4.0)), true,),
            ((3.0, -2.0), (5.0, 1.0)) => SegmentEvent(((5.0, 1.0), (3.0, 4.0)), true,),
            ((1.0, 1.0), (3.0, -2.0)) => SegmentEvent(((3.0, -2.0), (5.0, 1.0)), true),
            ((5.0, 1.0), (3.0, -2.0)) => SegmentEvent(((3.0, -2.0), (1.0, 1.0)), true),
            ((3.0, 4.0), (1.0, 1.0))  => SegmentEvent(((1.0, 1.0), (3.0, -2.0)), true),
            ((1.0, 1.0), (3.0, 4.0))  => SegmentEvent(((3.0, 4.0), (5.0, 1.0)), true,),
        )
        @test connections == expected
    end

    @testset "faces - cross" begin
        graph = Dict(
            (3.0, 1.0)  => [
                SegmentEvent(((3.0, 1.0), (1.0, 1.0))),
                SegmentEvent(((3.0, 1.0), (3.0, 4.0))),
                SegmentEvent(((3.0, 1.0), (6.0, 1.0))),
                SegmentEvent(((3.0, 1.0), (3.0, -2.0))),
            ],
            (1.0, 1.0)  => [SegmentEvent(((1.0, 1.0), (3.0, 1.0)))],
            (6.0, 1.0)  => [SegmentEvent(((6.0, 1.0), (3.0, 1.0)))],
            (3.0, 4.0)  => [SegmentEvent(((3.0, 4.0), (3.0, 1.0)))],
            (3.0, -2.0) => [SegmentEvent(((3.0, -2.0), (3.0, 1.0)))],
        )
        faces = compute_graph_faces(graph)
        expected = [[
            SegmentEvent(((3.0, 1.0), (3.0, 4.0)), true, true,),
            SegmentEvent(((3.0, 4.0), (3.0, 1.0)), true, true,),
            SegmentEvent(((3.0, 1.0), (1.0, 1.0)), true, true,),
            SegmentEvent(((1.0, 1.0), (3.0, 1.0)), true, true,),
            SegmentEvent(((3.0, 1.0), (3.0, -2.0)), true, true),
            SegmentEvent(((3.0, -2.0), (3.0, 1.0)), true, true),
            SegmentEvent(((3.0, 1.0), (6.0, 1.0)), true, true,),
            SegmentEvent(((6.0, 1.0), (3.0, 1.0)), true, true,),
        ]]
        @test faces == expected
    end

    @testset "faces - diamond" begin
        graph = Dict(
            (1.0, 1.0)  => [SegmentEvent(((1.0, 1.0), (3.0, 4.0))), SegmentEvent(((1.0, 1.0), (3.0, -2.0)))],
            (5.0, 1.0)  => [SegmentEvent(((5.0, 1.0), (3.0, 4.0))), SegmentEvent(((5.0, 1.0), (3.0, -2.0)))],
            (3.0, 4.0)  => [SegmentEvent(((3.0, 4.0), (5.0, 1.0))), SegmentEvent(((3.0, 4.0), (1.0, 1.0)))],
            (3.0, -2.0) => [SegmentEvent(((3.0, -2.0), (1.0, 1.0))), SegmentEvent(((3.0, -2.0), (5.0, 1.0)))],
        )
        faces = compute_graph_faces(graph)
        expected = [
            # exterior
            [
                SegmentEvent(((5.0, 1.0),(3.0, -2.0))),
                SegmentEvent(((3.0, -2.0), (1.0, 1.0))),
                SegmentEvent(((1.0, 1.0), (3.0, 4.0))),
                SegmentEvent(((3.0, 4.0), (5.0, 1.0))),
            ],
            # interior
            [
                SegmentEvent(((3.0, 4.0), (1.0, 1.0))),
                SegmentEvent(((1.0, 1.0), (3.0, -2.0))),
                SegmentEvent(((3.0, -2.0), (5.0, 1.0))),
                SegmentEvent(((5.0, 1.0),(3.0, 4.0))),
            ],
        ]
        @test faces == expected
    end

    @testset "faces - diamond with diagonal" begin
        graph = Dict(
            (1.0, 1.0)  => [SegmentEvent(((1.0, 1.0), (3.0, -2.0))), SegmentEvent(((1.0, 1.0), (5.0, 1.0))), SegmentEvent(((1.0, 1.0), (3.0, 4.0)))],
            (5.0, 1.0)  => [SegmentEvent(((5.0, 1.0), (3.0, -2.0))), SegmentEvent(((5.0, 1.0), (3.0, 4.0))), SegmentEvent(((5.0, 1.0), (1.0, 1.0)))],
            (3.0, 4.0)  => [SegmentEvent(((3.0, 4.0), (1.0, 1.0))), SegmentEvent(((3.0, 4.0), (5.0, 1.0)))],
            (3.0, -2.0) => [SegmentEvent(((3.0, -2.0), (5.0, 1.0))), SegmentEvent(((3.0, -2.0), (1.0, 1.0)))],
        )
        faces = compute_graph_faces(graph)
        expected = [
            # exterior
            [
                SegmentEvent(((5.0, 1.0),(3.0, -2.0))),
                SegmentEvent(((3.0, -2.0), (1.0, 1.0))),
                SegmentEvent(((1.0, 1.0), (3.0, 4.0))),
                SegmentEvent(((3.0, 4.0), (5.0, 1.0))),
            ],
            # interior
            [
                SegmentEvent(((3.0, 4.0), (1.0, 1.0))),
                SegmentEvent(((1.0, 1.0), (5.0, 1.0))),
                SegmentEvent(((5.0, 1.0), (3.0, 4.0))),
            ],
            [
                SegmentEvent(((5.0, 1.0), (1.0, 1.0))),
                SegmentEvent(((1.0, 1.0), (3.0, -2.0))),
                SegmentEvent(((3.0, -2.0), (5.0, 1.0))),
            ],
        ]
        @test faces == expected
    end

    @testset "faces - improper" begin
        # two separate segments
        graph = directed_graph_from_segments(improper[[1, 2, 4, 5, 6, 7, 8]])
        faces = compute_graph_faces(graph)
        expected = [
            # interior diamond
            [
                SegmentEvent(((5.0, 1.0), (3.0, -2.0))),
                SegmentEvent(((3.0, -2.0), (1.0, 1.0))),
                SegmentEvent(((1.0, 1.0), (3.0, 4.0))),
                SegmentEvent(((3.0, 4.0), (5.0, 1.0)))
            ],
            # exterior diamond
            [
                SegmentEvent(((5.0, 1.0), (3.0, 4.0)), ),
                SegmentEvent(((3.0, 4.0), (1.0, 1.0)), ),
                SegmentEvent(((1.0, 1.0), (3.0, -2.0)),),
                SegmentEvent(((3.0, -2.0), (5.0, 1.0)),),
            ],
            # interior triangle
            [
                SegmentEvent(((8.0, 3.0), (8.0, 1.0))),
                SegmentEvent(((8.0, 1.0), (7.0, 2.0))),
                SegmentEvent(((7.0, 2.0), (8.0, 3.0))),
            ],
            # exterior triangle
            [
                SegmentEvent(((8.0, 1.0), (8.0, 3.0))),
                SegmentEvent(((8.0, 3.0), (7.0, 2.0))),
                SegmentEvent(((7.0, 2.0), (8.0, 1.0))),
            ],
        ]
        @test faces == expected
        # with tail segment
        graph = directed_graph_from_segments(improper[[1, 2, 3, 4, 5, 7, 8]])
        faces = compute_graph_faces(graph)
        expected = [
            # interior
            [
                SegmentEvent(((5.0, 1.0), (3.0, -2.0)), true),
                SegmentEvent(((3.0, -2.0), (1.0, 1.0)), true),
                SegmentEvent(((1.0, 1.0), (3.0, 4.0)),  true),
                SegmentEvent(((3.0, 4.0), (5.0, 1.0)),  true),
            ],
            # exterior
            [
                SegmentEvent(((5.0, 1.0), (7.0, 2.0)), true),    
                SegmentEvent(((7.0, 2.0), (8.0, 3.0)), true),
                SegmentEvent(((8.0, 3.0), (8.0, 1.0)), true),
                SegmentEvent(((8.0, 1.0), (8.0, 3.0)), true),
                SegmentEvent(((8.0, 3.0), (7.0, 2.0)), true),
                SegmentEvent(((7.0, 2.0), (5.0, 1.0)), true),
                SegmentEvent(((5.0, 1.0), (3.0, 4.0)), true),
                SegmentEvent(((3.0, 4.0), (1.0, 1.0)), true),
                SegmentEvent(((1.0, 1.0), (3.0, -2.0)), true),
                SegmentEvent(((3.0, -2.0), (5.0, 1.0)), true),
            ]
        ]
        @test faces == expected
        # two connected segments
        graph = directed_graph_from_segments(improper)
        faces = compute_graph_faces(graph)
        expected = [
            # interior diamond
            [
                SegmentEvent(((5.0, 1.0), (3.0, -2.0)), true),
                SegmentEvent(((3.0, -2.0), (1.0, 1.0)), true),
                SegmentEvent(((1.0, 1.0), (3.0, 4.0)),  true),
                SegmentEvent(((3.0, 4.0), (5.0, 1.0)),  true),
            ],
            # exterior
            [
                SegmentEvent(((5.0, 1.0), (7.0, 2.0)), true),
                SegmentEvent(((7.0, 2.0), (8.0, 1.0)), true,),
                SegmentEvent(((8.0, 1.0), (8.0, 3.0)), true,),
                SegmentEvent(((8.0, 3.0), (7.0, 2.0)), true,),
                SegmentEvent(((7.0, 2.0), (5.0, 1.0)), true,),
                SegmentEvent(((5.0, 1.0), (3.0, 4.0)), true,),
                SegmentEvent(((3.0, 4.0), (1.0, 1.0)), true,),
                SegmentEvent(((1.0, 1.0), (3.0, -2.0)), true),
                SegmentEvent(((3.0, -2.0), (5.0, 1.0)), true),
            ],
            # interior triangle
            [
                SegmentEvent(((8.0, 3.0), (8.0, 1.0)), true),
                SegmentEvent(((8.0, 1.0), (7.0, 2.0)), true),
                SegmentEvent(((7.0, 2.0), (8.0, 3.0)), true),
            ],
        ]
        @test faces == expected
    end
end