"""
    directed_graph_from_segments(segments; digits=6)

Create a directed graph of half-edges from a list of segments.
Each segment creates two half-edges.
Each node is a point with the list of connected half edges,
where each half edge points to another node and are sorted clockwise around the point.
"""
function directed_graph_from_segments(events::AbstractVector{<:SegmentEvent{T}}; digits::Int=6) where {T}
    graph = Dict{Point2D{T}, Vector{SegmentEvent{T}}}() # point => half_edges
    for ev in events
        # tolerance
        segment = (round.(ev.segment[1], digits=digits), round.(ev.segment[2], digits=digits))
        # create segments
        current = SegmentEvent(segment, true, true, ev.self_annotations, SegmentAnnotations())
        other = SegmentEvent(reverse(segment), true, true, ev.self_annotations, SegmentAnnotations())
        push!(get!(() -> SegmentEvent{T}[], graph, current[1]), current)
        push!(get!(() -> SegmentEvent{T}[], graph, other[1]), other)
    end
    # sort clockwise around the node
    for node in keys(graph)
        sort!(graph[node], lt=(seg1, seg2) -> !isless_polar_angle(seg1[2], seg2[2], seg1[1]))
    end
    graph
end

"""
    map_connections(graph::Dict)

If Z is the next neighbor clockwise around the node, then the half-edge
that comes before (X, Y) is the half-edge (Z, X).
"""
function map_connections(graph::Dict{Tuple{T, T}, <:AbstractVector{SegmentEvent{T}}}) where {T}
    connections = Dict{Segment2D{T}, SegmentEvent{T}}()
    for (node, half_edges) in graph
        for (idx, edge) in enumerate(half_edges)
            next_idx = idx == length(half_edges) ? 1 : idx + 1
            next_edge = half_edges[next_idx]
            connections[(next_edge[2], node)] = edge
        end
    end
    connections
end

function walk_connections(connections::Dict{Segment2D{T}, SegmentEvent{T}}) where {T}
    chains = Vector{SegmentEvent{T}}[]
    visited = Set{Segment2D{T}}()
    for (start, next_edge) in connections
        @debug("At $start")
        if !(start in visited)
            @debug("   at $(next_edge.segment)")
            push!(visited, start)
            chain = [next_edge]
            while next_edge.segment != start
                if next_edge.segment in visited
                    # Theoretically this should not happen
                    # This avoid the infinite loop in the unlikely scenario it does
                    @warn "Cycle detected at edge $(next_edge.segment)"
                    break
                end
                @debug("   at $(next_edge.segment)")
                push!(visited, next_edge.segment)
                next_edge = connections[next_edge.segment]
                push!(chain, next_edge)
            end
            push!(chains, chain)
        end
    end
    chains
end

"""
    compute_graph_faces(graph::Dict)

Compute all faces in a directed graph.

The graph is dictionary of `point => half-edges`, with half-edges sorted clockwise
around the point.

This will return all exterior faces as counter-clockwise, and all
interior faces as clockwise.

Source
- https://stackoverflow.com/a/67169488. Posted by templatetypedef
"""
function compute_graph_faces(graph::Dict{Tuple{T, T}, <:AbstractVector{SegmentEvent{T}}}) where {T}
    connections = map_connections(graph)
    walk_connections(connections)
end