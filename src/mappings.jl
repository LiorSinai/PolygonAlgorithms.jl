## Mappings between different representations of polygons

"""
    is_hole(segments::Vector{<:AnnotatedSegment}, [counter_clockwise])

A polygon is a hole if its outer face is filled but not its inner face.

For robustness, the annotations of all segments are checked
and a decision is made based on the majority outcome.

If the direction `counter_clockwise` is not given, the segments must be in a chain.
That is, `segments_{i}[2] == segments_{i+1}[2]` for all segments.

For a counter-clockwise polygon, the inner face is to the left of each segment.
For a clockwise polygon, the inner face is to the right of each segment.

Filled above/below annotations are converted to left/right annotations using the normal
vector.
This is `(-Δy, Δx)` for a counter-clockwise polygon and `(Δy, -Δx)` for a clockwise polygon.
"""
function is_hole(polygon::AbstractVector{<:AnnotatedSegment{T}}) where {T}
    n = length(polygon)
    @assert all(polygon[i][2] == polygon[(i % n) + 1][1] for i in 1:n)
    points = [segment[1] for segment in polygon]
    counter_clockwise = is_counter_clockwise(points)
    is_hole(polygon, counter_clockwise)
end

function is_hole(polygon::AbstractVector{<:AnnotatedSegment{T}}, counter_clockwise::Bool) where {T}
    # for robustness, use majority voting for all points
    votes_face = 0
    votes_hole = 0
    for segment in polygon
        Δx = segment[2][1] - segment[1][1]
        ann = segment.self_annotations
        if (Δx == 0) || (ann.fill_above == ann.fill_below)
            # skip vertical segments or filled both sides or filled on neither
            continue
        end
        if counter_clockwise
            # inner face on left with normal (-Δy, Δx)
            # Therefore Δx > 0 ? above : below
            filled_inner = Δx > 0 ? ann.fill_above : ann.fill_below
        else
            # inner face on right with normal (Δy, -Δx)
            # Therefore -Δx < 0 ? below : above
            filled_inner = Δx > 0 ? ann.fill_below : ann.fill_above
        end
        if filled_inner
            votes_face += 1
        else
            votes_hole += 1
        end
    end
    votes_hole > votes_face
end

function events_to_paths(
    events::Union{Vector{<:SegmentEvent}, Vector{<:AnnotatedSegment}}
    ; digits::Integer=6
    )
    # event → segments
    graph = directed_graph_from_segments(events; digits=digits)
    # faces → paths
    faces = compute_graph_faces(graph)
    polygons = map(segments -> map(event -> event[1], segments), faces)
    # paths → exteriors, holes
    exteriors = is_counter_clockwise.(polygons)
    not_holes = .!is_hole.(faces[exteriors], true) # ignore exteriors of holes
    interiors = .!exteriors
    holes = is_hole.(faces[interiors], false) # ignore repeated interiors
    polygons[exteriors][not_holes], polygons[interiors][holes]
end

"""
    match_holes_polygons(polygons::Vector{<:Polygon}, holes::Vector{<:Tuple})

An algorithm for matching holes to polygons.
Returns the index of each parent for each hole.

Assumes that the polygons do not intersect.

For every hole, match to a polygon that contains the hole.
If there are multiple polygons possible, the polygon with the least area is chosen.
If no polygons are found, return the hole as a `Polygon`.

In the best case there is one polygon or one hole. Then this runs in `O(1)` time.
In the worst case, none of the holes match to a polygon.
Then this runs in `O(phn)` time where `p` is the number polygons,
`h` is the number of holes and `n` is the average number of vertices defining each polygon.
"""
function match_holes_polygons(
    polygons::Vector{<:Polygon},
    holes::Vector{<:Path2D}
    ; atol::AbstractFloat=default_atol
    )
    if length(polygons) == 1
        return fill(1, length(holes))
    elseif isempty(holes)
        return Int[]
    end
    areas = map(area_polygon, polygons)
    # sort by ascending areas. Therefore smallest parent is matched first.
    # Trade off is the hole may be tried in many smaller polygons first.
    idxs = sortperm(areas)
    parents = zeros(Int, length(holes))
    for (idx_h, candidate) in enumerate(holes)
        for (idx_p, parent) in zip(idxs, polygons[idxs])
            # Assume that no segments intersect after the Martínez-Rueda algorithm.
            # Then only need to check a point.
            found = contains(parent.exterior, candidate[1]; atol=atol, on_border_is_inside=true)
            if found
                parents[idx_h] = idx_p
                break
            end
        end
    end
    parents
end

function paths_to_polygons(
    exteriors::Vector{<:Vector{<:Point2D}},
    holes::Vector{<:Vector{<:Point2D}},
    ; atol::AbstractFloat=default_rtol
    )
    polygons = Polygon.(exteriors)
    parents = match_holes_polygons(polygons, holes; atol=atol)
    for (idx, hole) in zip(parents, holes)
        if idx != 0
            push!(polygons[idx].holes, hole)
        end
    end
    polygons
end

function events_to_polygons(
    events::Union{Vector{<:SegmentEvent}, Vector{<:AnnotatedSegment}}
    ; digits::Integer=6, atol::AbstractFloat=default_atol
    )
    exteriors, holes = events_to_paths(events; digits=digits)
    paths_to_polygons(exteriors, holes; atol=atol)
end
