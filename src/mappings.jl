## Mappings between different representations of polygons

"""
    is_hole(segments::Vector{<:AnnotatedSegment}, [counter_clockwise])

A polygon is a hole if its outer face is filled but not its inner face.

For robustness, the annotations of all segments are checked
and a decision is made based on the majority outcome.
The time complexity is therefore `O(n)` where `n` is the number of segments.

If the direction `counter_clockwise` is not given, the segments must be in a chain.
That is, `segments_{i}[2] == segments_{i+1}[1]` for all segments.

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
    # only the direction is needed, so don't treat thin polygons as zero area
    counter_clockwise = is_counter_clockwise(points; atol=eps(T))
    is_hole(polygon, counter_clockwise)
end

function is_hole(polygon::AbstractVector{<:AnnotatedSegment{T}}, counter_clockwise::Bool) where {T}
    # for robustness, use majority voting for all points
    votes_face = 0
    votes_hole = 0
    for segment in polygon
        Δx = segment[2][1] - segment[1][1]
        ann = segment.self_annotations
        if is_vertical(segment.segment) || (ann.fill_above == ann.fill_below)
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

"""
Segments are converted back to paths and polygons by computing the faces of a graph.

For every non-zero area polygon, every area in the graph is counted twice.
The graph will always have one counter-clockwise exterior and one or more clockwise
interiors. 

The following selection strategies are used to return half the faces with no duplicates:
- `MERGE_FACES`: selects all exterior faces that are not holes with interior faces that are holes, resulting in connected graphs being merged together. This returns fewer, larger polygons with holes.
- `SPLIT_FACES`: As exteriors, select all interior faces that are not holes that are on the exterior and reverse to make them counter-clockwise. As holes, select interior faces that are holes that are not on the exterior. This splits connected graphs across multiple interiors and returns smaller polygons with fewer explicit holes.
- `CLOCKWISE_FACES`: all clockwise faces.
- `COUNTER_CLOCKWISE_FACES`: all counter-clockwise faces.
"""
@enum FaceSelectionStrategy MERGE_FACES SPLIT_FACES CLOCKWISE_FACES COUNTER_CLOCKWISE_FACES

"""
    segments_to_paths(segments; atol=default_atol, face_selection=SPLIT_FACES)

Assumes there are no intersections between the segments.

The segments are joined in a graph. 
Each face with non-zero area will return two faces, a counter-clockwise face with positive area (exterior)
and a clockwise face with negative area (interior).
Faces thinner than `atol` are considered to have zero area.
This is estimated as `|area| ≤ atol * perimeter / 2`. See `_first_moment_with_error`.

Returns `(exteriors, holes, lines)`, where `lines` are the zero-area faces.

The following selection strategies are used to return half the faces with no duplicates:
- `MERGE_FACES`: selects all exterior faces that are not holes with interior faces that are holes, resulting in connected graphs being merged together. This returns fewer, larger polygons with holes.
- `SPLIT_FACES`: As exteriors, select all interior faces that are not holes that are on the exterior and reverse to make them counter-clockwise. As holes, select interior faces that are holes that are not on the exterior. This splits connected graphs across multiple interiors and returns smaller polygons with fewer explicit holes.
- `CLOCKWISE_FACES`: all clockwise faces.
- `COUNTER_CLOCKWISE_FACES`: all counter-clockwise faces.

Time complexity is `O(n*log(n))` in the worst case where `n` is the number of segments, and `O(n)` if the
number of segments at each vertex is bounded. See `directed_graph_from_segments` and `compute_graph_faces`.
Classifying the faces is `O(n)` because the faces contain `2n` half-edges in total.

`SPLIT_FACES` iteratively expands the exterior from the faces that touch it.
This usually takes a few passes, but in the worst case it is `O(f*n)` where `f` is the number of interior faces,
if each pass only adds one face.
"""
function segments_to_paths(
    segments::Union{AbstractVector{<:SegmentEvent{T}}, AbstractVector{<:AnnotatedSegment{T}}}
    ;
    atol::AbstractFloat=default_atol,
    face_selection::FaceSelectionStrategy=SPLIT_FACES,
    ) where T
    # event → segments
    graph = directed_graph_from_segments(segments; digits=decimal_tolerance(atol))
    # faces → paths
    faces = compute_graph_faces(graph)
    polygons = map(segments -> map(event -> event[1], segments), faces)
    # A face is zero area if it is thinner than atol. See _first_moment_with_error.
    moments_errors = map(p -> _first_moment_with_error(p; atol=atol), polygons)
    are_exteriors = map(((moment, error),) -> moment > error, moments_errors) # counter-clockwise
    are_interiors = map(((moment, error),) -> moment < -error, moments_errors) # clockwise
    are_zeros = map(((moment, error),) -> abs(moment) <= error, moments_errors)
    # paths → exteriors, holes
    if face_selection == CLOCKWISE_FACES
        exteriors = polygons[are_interiors]
        holes = empty(exteriors)
    elseif face_selection == COUNTER_CLOCKWISE_FACES
        exteriors = polygons[are_exteriors]
        holes = empty(exteriors)
    else
        moments = first_moment.(polygons)
        are_exteriors = moments .> eps(T) # counter-clockwise
        are_interiors = moments .< -eps(T) # clockwise
        not_holes = .!is_hole.(faces[are_exteriors], true) # ignore exteriors of holes
        if face_selection == MERGE_FACES
            are_holes = is_hole.(faces[are_interiors], false) # ignore repeated interiors
            exteriors = polygons[are_exteriors][not_holes]
            holes = polygons[are_interiors][are_holes]
        else # SPLIT_FACES
            are_holes = is_hole.(faces[are_interiors], false)
            interiors = polygons[are_interiors]
            inner_faces = interiors[.!are_holes]
            holes = interiors[are_holes]
            # for exterior faces, choose only inner faces on the exterior
            # as well as any connected exterior faces that are possibly inside the exteriors in holes
            exterior_points = Set{Point2D{T}}()
            for points in polygons[are_exteriors][not_holes]
                union!(exterior_points, points)
            end
            exteriors = empty(inner_faces)
            expand_frontier = true
            while expand_frontier
                expand_frontier = false
                for idx in length(inner_faces):-1:1
                    candidate = inner_faces[idx]
                    if any(pt -> (pt in exterior_points), candidate)
                        exterior = popat!(inner_faces, idx)
                        push!(exteriors, reverse!(exterior))
                        union!(exterior_points, exterior)
                        expand_frontier = true
                    end
                end
            end
            # for holes, choose only inner faces not on the exterior
            # because holes on the exterior are self-evident 
            filter!(pts -> !any(pt -> pt in exterior_points, pts), holes)
        end
    end
    exteriors, holes, polygons[are_zeros]
end

"""
    match_interiors(polygons; atol=default_atol)

For each polygon, find the index of the largest polygon that contains it, or `0` if there is none.
Polygons that are already contained in another polygon cannot be parents.
This results in alternating levels of exteriors and holes, where the top level polygons are exteriors.

Assumes that the polygons do not intersect, so only one point needs to be checked.
Returns a vector of parent indices.

Time complexity is `O(p*log(p) + p²n)` where `p` is the number of polygons and `n` is the average number
of vertices per polygon, because each polygon may be checked against every larger polygon.
In the worst case, where many vertices of a polygon lie on the border of a candidate parent,
finding a vertex to check costs `O(n²)` per pair instead of `O(n)`.
"""
function match_interiors(
    polygons::Vector{<:Path2D}
    ; atol::AbstractFloat=default_atol,
    )
    areas = map(area_polygon, polygons)
    # sort by descending areas. Therefore largest parent is matched first.
    idxs = sortperm(areas, rev=true)
    parents = zeros(Int, length(polygons))
    for (rank, idx1) in enumerate(idxs)
        polygon = polygons[idx1]
        for idx2 in idxs[1:(rank - 1)] # can only be in larger polygons
            if parents[idx2] != 0
                continue
            end
            # Assume that no segments intersect. Then only need to check a single point
            parent = polygons[idx2]
            j = 1
            while (j < length(polygon)) && on_border(parent, polygon[j]; atol=atol)
                j += 1
            end
            found = contains(parent, polygon[j]; atol=atol, on_border_is_inside=false)
            if found
                parents[idx1] = idx2
                break
            end
        end
    end
    parents
end

"""
    match_holes_polygons(polygons::Vector, holes::Vector; atol=default_atol)

An algorithm for matching holes to polygons.
Returns the index of the parent for each hole, or `0` if no parent was found.

Assumes that the polygons do not intersect.

For every hole, match to a polygon that contains the hole.
If there are multiple polygons possible, the polygon with the least area is chosen.

Let `p` be the number of polygons, `h` the number of holes and `n` the average number of vertices
per polygon.

In the best case there is one polygon or no holes. Then this runs in `O(h)` time,
because every hole is assigned to the one polygon without checking.
Otherwise, calculating and sorting the areas costs `O(pn + p*log(p))`
and each hole is checked against the polygons in order of increasing area, costing `O(n)` per check.
In the worst case none of the holes match to a polygon, and this runs in `O(phn)` time.
"""
function match_holes_polygons(
    polygons::Vector{<:Path2D},
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
        for idx_p in idxs
            # Assume that no segments intersect.
            # Then only need to check a point.
            found = contains(polygons[idx_p], candidate[1]; atol=atol, on_border_is_inside=true)
            if found
                parents[idx_h] = idx_p
                break
            end
        end
    end
    parents
end

"""
    paths_to_polygons(exteriors, holes; atol=default_atol)

Create a `Polygon` for each exterior and assign each hole to its parent using `match_holes_polygons`.
Holes without a parent are dropped with a warning.

Time complexity is dominated by `match_holes_polygons`, which is `O(phn)` in the worst case where `p` is the
number of exteriors, `h` is the number of holes and `n` is the average number of vertices per exterior.
"""
function paths_to_polygons(
    exteriors::Vector{<:Vector{<:Point2D}},
    holes::Vector{<:Vector{<:Point2D}},
    lines::Vector{<:Vector{<:Point2D}},
    ; atol::AbstractFloat=default_atol
    )
    polygons = Polygon.(exteriors)
    parents = match_holes_polygons(exteriors, holes; atol=atol)
    for (idx, hole) in zip(parents, holes)
        if idx != 0
            push!(polygons[idx].holes, hole)
        else
            @warn "Hole $hole has no parent and will be dropped."
        end
    end
    vcat(polygons, Polygon.(lines))
end

"""
    segments_to_polygons(segments; atol=default_atol, face_selection=SPLIT_FACES)

Convert segments to a list of `Polygon`s with holes.

The segments are first converted to exteriors and holes with `segments_to_paths`
and then each hole is matched to its parent with `paths_to_polygons`.
See `segments_to_paths` for details on `face_selection`.

Assumes there are no intersections between the segments.

Time complexity is the sum of `segments_to_paths` and `paths_to_polygons`:
`O(n*log(n))` to compute the faces from `n` segments, plus up to `O(phm)` to match `h` holes to `p` exteriors
with an average of `m` vertices each.
"""
function segments_to_polygons(
    segments::Union{AbstractVector{<:SegmentEvent}, AbstractVector{<:AnnotatedSegment}}
    ; 
    atol::AbstractFloat=default_atol,
    face_selection::FaceSelectionStrategy=SPLIT_FACES,
    )
    exteriors, holes, lines = segments_to_paths(segments; atol=atol, face_selection=face_selection)
    paths_to_polygons(exteriors, holes, lines; atol=atol)
end
