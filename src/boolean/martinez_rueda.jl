@enum AnnotationFill BLANK=0 ABOVE=1 BELOW=2 EMPTY=3

"""
    is_hole(polygon::Vector{<:SegmentEvent})

A necessary and sufficient condition for a polygon to be classified as a hole is that
its lowest segment must be filled below and not above.
"""
function is_hole(polygon::AbstractVector{<:SegmentEvent{T}}) where {T}
    points = [segment[1] for segment in polygon]
    counter_clockwise = is_counter_clockwise(points)
    # for robustness, use majority voting for all points
    votes_face = 0
    votes_hole = 0
    for segment in polygon
        Δx = segment[2][1] - segment[1][1]
        ann = segment.self_annotations
        if (Δx == 0) || (ann.fill_above == ann.fill_below)
            # skip vertical segments || filled both sides or neither
            continue
        end
        if counter_clockwise
            # inner normal is left normal (-Δy, Δx)
            # Therefore Δx > 0 ? above : below
            filled_inner = Δx > 0 ? ann.fill_above : ann.fill_below
        else
            # inner normal is right normal (Δy, -Δx)
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

function segments_to_paths(segments::Vector{<:SegmentEvent}; digits::Integer=6)
    # segments → faces → paths
    graph = directed_graph_from_segments(segments; digits=digits)
    faces = compute_graph_faces(graph)
    polygons = map(segments -> map(event -> event.point, segments), faces)
    # paths → exteriors, holes
    exteriors = is_counter_clockwise.(polygons)
    not_holes = .!is_hole.(faces[exteriors]) # ignore exteriors of holes
    interiors = .!exteriors
    holes = is_hole.(faces[interiors]) # ignore repeated interiors
    polygons[exteriors][not_holes], polygons[interiors][holes]
end

"""
    martinez_rueda_algorithm(
    selection_criteria, subject, others...
    ; atol=default_atol, rtol=default_rtol
    )

The Martínez-Rueda-Feito polygon clipping algorithm.
Returns regions and edges of intersection.
It runs in `O((n+m+k)log(n+m))` time where `n` and `m` are the number of vertices of `polygon1` 
and `polygon2` respectively and `k` is the total number of intersections between all segments.
Use `intersect_convex` for convex polygons for an `O(n+m)` algorithm.

The input polygons can be:
- A list of points: `Vector{Tuple{Float64, Float64}}`.
- Polygons: `Vector{Polygon{T}}`.
- A Segment event queue: `Vector{<:SegmentEvent{Float64}}`. The core algorithm uses this representation.
- The `subject` can also be a list of polygons. This enables the output to be an input. 
    Note that in this case the `subject` list is treated as a single polygon.
    If any of the polygons overlap, this is equivalent to passing a self-intersecting polygon and
    some areas might be classified as holes according to the even-odd rule

Description:
- Operates at a segment level and is an extension of the Bentley-Ottman line intersection algorithm.
    Segments are scanned from left to right, bottom to top. 
- The key assumption is that only the segments immediately above and below the current segment need to be inspected for intersections.
    This makes the algorithm fast but also sensitive to determining these segments correctly.
- The segment that is immediately below (or empty space) is used to determine the fill annotations for the current segment.
- Once all annotations are done, the desired segments can be selected that match a given criteria.
- These segments are then chained together to form the polygons.

Limitations
1. It is sensitive to numeric inaccuracies e.g. a line that is almost vertical or tiny regions 
    of intersection.

References 
- paper: https://www.researchgate.net/publication/220163820_A_new_algorithm_for_computing_Boolean_operations_on_polygons
- article: https://sean.fun/a/polygon-clipping-pt2/
- article source code: https://github.com/velipso/polybooljs
"""
function martinez_rueda_algorithm(
    selection_criteria::Vector{AnnotationFill},
    subject::Path2D{T},
    others::Vararg{Path2D{T}},
    ; atol::AbstractFloat=default_atol, options...
    ) where T
    event_queue_base = convert_to_event_queue(subject; primary=true, atol=atol)
    event_queue_others = map(p -> convert_to_event_queue(p; primary=false, atol=atol), others)
    segments = martinez_rueda_algorithm(
        selection_criteria, event_queue_base, event_queue_others...; atol=atol, options...
    )
    exteriors, holes = segments_to_paths(segments; digits=significant_digits(atol))
    vcat(exteriors, holes)
end

# Multiple subjects
function martinez_rueda_algorithm(
    selection_criteria::Vector{AnnotationFill},
    subjects::AbstractVector{<:Path2D{T}},
    others::Vararg{Path2D{T}},
    ; atol::AbstractFloat=default_atol, options...
    ) where T
    subject_queue = SegmentEvent{T}[]
    map(p -> convert_to_event_queue!(subject_queue, p; primary=true, atol=atol), subjects)
    event_queue_others = map(p -> convert_to_event_queue(p; primary=false, atol=atol), others)
    segments = martinez_rueda_algorithm(
        selection_criteria, subject_queue, event_queue_others...; atol=atol, options...
    )
    exteriors, holes = segments_to_paths(segments; digits=significant_digits(atol))
    vcat(exteriors, holes)
end

# Polygons with holes
function martinez_rueda_algorithm(
    selection_criteria::Vector{AnnotationFill},
    subject::Polygon{T},
    others::Vararg{Polygon{T}},
    ; atol::AbstractFloat=default_atol, options...
    ) where T
    event_queue_base = convert_to_event_queue(subject.exterior; primary=true, atol=atol)
    for hole in subject.holes
        convert_to_event_queue!(event_queue_base, hole; primary=true, atol=atol)
    end
    event_queue_others = map(p -> convert_to_event_queue(p.exterior; primary=false, atol=atol), others)
    for (queue, other) in zip(event_queue_others, others)
        for hole in other.holes
            convert_to_event_queue!(queue, hole; primary=false, atol=atol)
        end
    end
    segments = martinez_rueda_algorithm(
        selection_criteria, event_queue_base, event_queue_others...; atol=atol, options...
    )
    exteriors, holes = segments_to_paths(segments; digits=significant_digits(atol))
    paths_to_polygons(exteriors, holes; atol=atol)
end

# Multiple subjects with holes
function martinez_rueda_algorithm(
    selection_criteria::Vector{AnnotationFill},
    subjects::AbstractVector{<:Polygon{T}},
    clips::Vararg{Polygon{T}},
    ; atol::AbstractFloat=default_atol, options...
    ) where T
    subject_queue = SegmentEvent{T}[]
    map(p -> convert_to_event_queue!(subject_queue, p.exterior; primary=true, atol=atol), subjects)
    for subject in subjects
        for hole in subject.holes
            convert_to_event_queue!(subject_queue, hole; primary=true, atol=atol)
        end
    end
    event_queue_clips = map(p -> convert_to_event_queue(p.exterior; primary=false, atol=atol), clips)
    for (queue, other) in zip(event_queue_clips, clips)
        for hole in other.holes
            convert_to_event_queue!(queue, hole; primary=false, atol=atol)
        end
    end
    segments = martinez_rueda_algorithm(
        selection_criteria, subject_queue, event_queue_clips...; atol=atol, options...
    )
    exteriors, holes = segments_to_paths(segments; digits=significant_digits(atol))
    paths_to_polygons(exteriors, holes; atol=atol)
end

# Core algorithm: SegmentEvent input
function martinez_rueda_algorithm(
    selection_criteria::Vector{AnnotationFill},
    subject::Vector{<:SegmentEvent{T}},
    polygons::Vararg{Vector{<:SegmentEvent{T}}},
    ; atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_rtol
    ) where T
    base_annotated_segments = event_loop!(subject; self_intersection=true, atol=atol, rtol=rtol)
    for polygon in polygons
        annotated_segments = event_loop!(polygon; self_intersection=true, atol=atol, rtol=rtol)
        queue = SegmentEvent{T}[]
        for ev in vcat(base_annotated_segments, annotated_segments)
            add_annotated_segment!(queue, ev)
        end
        combined_annotated_segments = event_loop!(queue; self_intersection=false, atol=atol, rtol=rtol)
        # for consistent reporting, swap annotations so that self annotations are always the primary
        for ev in combined_annotated_segments
            if !ev.primary
                temp = ev.self_annotations
                ev.self_annotations = ev.other_annotations
                ev.other_annotations = temp
            end
        end
        base_annotated_segments = apply_selection_criteria(combined_annotated_segments, selection_criteria)
    end
    base_annotated_segments
end

function add_annotated_segment!(queue::Vector{<:SegmentEvent}, ev::SegmentEvent)
    pt1 = ev.segment[1]
    pt2 = ev.segment[2]
    forward = _compare_points(pt1, pt2)
    if forward == 0
        return queue # zero length segment
    end
    start = forward < 0 ? pt1 : pt2
    end_ = forward < 0 ? pt2 : pt1
    segment = (start, end_)
    add_segment_event!(queue, segment, ev.primary, ev.self_annotations, ev.other_annotations)
end

#############################################################
##                    Event Loop                           ##
#############################################################

function event_loop!(
    queue::Vector{SegmentEvent{T}}
    ; self_intersection::Bool, atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_atol
    ) where T # eventLoop
    annotated_segments = SegmentEvent{T}[]
    sweep_status = SegmentEvent{T}[] # current events in a vertical line, top to bottom.
    while !(isempty(queue))
        head = queue[1]
        queue_length = length(queue)
        status_length = length(sweep_status)
        @debug("[event_loop!] ($(queue_length), $(status_length)): $(head)")
        if head.is_start # then check for intersections and add to sweep status
            idx = find_transition(sweep_status, head; atol=atol)
            above = idx == 1 ? nothing : sweep_status[idx - 1]
            below = (idx > length(sweep_status)) ? nothing : sweep_status[idx]
            @debug("[event_loop!] transition idx=$idx")
            @debug("[event_loop!] above=$above")
            @debug("[event_loop!] below=$below")
            check_and_divide_intersection!(queue, head, above, self_intersection; atol=atol, rtol=rtol)
            if queue[1] != head
                continue # either head was removed or something was inserted ahead of it
            end
            check_and_divide_intersection!(queue, head, below, self_intersection; atol=atol, rtol=rtol)
            if queue[1] != head
                continue # either head was removed or something was inserted ahead of it
            end
            if self_intersection
                calculate_self_annotations!(head, below)
            else
                calculate_other_annotations!(head, below)
            end
            insert!(sweep_status, idx, head)
        else # event is ending, so remove it from the status
            idx = find_transition(sweep_status, head.other; atol=atol)
            if !(0 < idx <= length(sweep_status) && sweep_status[idx] === head.other)
                @warn "$(head.other) was not in the expected location in the sweep status. " * 
                    "Falling back to linear search. This might result in incorrect annotations and hence open chains."
                idx = findfirst(x -> x === head.other, sweep_status)
                @assert(
                    !isnothing(idx),
                    "$(head.other) is missing from the sweep_status. The start event should always be processed before the end event."
                )
            end
            if (idx != 1) && (idx != length(sweep_status))
                # there will be 2 new adjacent edges, so check the intersection between them
                check_and_divide_intersection!(
                    queue, sweep_status[idx - 1], sweep_status[idx + 1], self_intersection
                    ; atol=atol, rtol=rtol)
            end
            push!(annotated_segments, copy_segment(head.other, head.other.primary))
            popat!(sweep_status, idx)
        end
        popfirst!(queue)
    end
    annotated_segments
end

function check_and_divide_intersection!(
    queue::Vector{<:SegmentEvent}, ev1::SegmentEvent, ev2::Nothing, self_intersection::Bool
    ; atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_rtol
    )
    queue
end

function check_and_divide_intersection!(
    queue::Vector{<:SegmentEvent},
    ev1::SegmentEvent,
    ev2::SegmentEvent,
    self_intersection::Bool
    ; atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_rtol
    )
    pt = intersect_geometry(ev1.segment, ev2.segment; atol=atol, rtol=rtol)
    if isnothing(pt)
        @debug("no intersection or parallel lines at $ev1 -- $ev2")
        # Lines might be on top of each other 
        ori2_start = get_orientation(ev1.segment[1], ev1.segment[2], ev2.segment[1]; atol=atol)
        ori2_end = get_orientation(ev1.segment[1], ev1.segment[2], ev2.segment[2]; atol=atol)
        if (ori2_start == COLINEAR) && (ori2_end == COLINEAR)
            divide_coincident_intersection!(queue, ev1, ev2, self_intersection; atol=atol)
        end
        return queue
    else
        divide_intersection!(queue, ev1, ev2, pt; atol=atol)
    end
end

function divide_intersection!(queue::Vector{<:SegmentEvent}, ev1::SegmentEvent, ev2::SegmentEvent, pt::Nothing; atol=1e-6)
    queue
end

function divide_intersection!(
    queue::Vector{<:SegmentEvent},
    ev1::SegmentEvent,
    ev2::SegmentEvent,
    pt::Point2D
    ; atol::AbstractFloat=default_atol
    ) # checkIntersection
    @debug("[divide_intersection!] $(ev1.segment) -- $(ev2.segment) at $(pt)")
    at_start1, at_end1, along1 = classify_intersection(ev1.segment, pt; atol=atol)
    at_start2, at_end2, along2 = classify_intersection(ev2.segment, pt; atol=atol)
    @debug("[divide_intersection!] $at_start1 $at_end1 $along1")
    @debug("[divide_intersection!] $at_start2 $at_end2 $along2")
    if along1 && along2
        divide_event!(queue, ev1, pt; atol=atol)
        divide_event!(queue, ev2, pt; atol=atol)
    elseif along1
        if at_start2
            divide_event!(queue, ev1, ev2.segment[1]; atol=atol)
        elseif at_end2
            divide_event!(queue, ev1, ev2.segment[2]; atol=atol)
        end
    elseif along2
        if at_start1
            divide_event!(queue, ev2, ev1.segment[1]; atol=atol)
        elseif at_end1
            divide_event!(queue, ev2, ev1.segment[2]; atol=atol)
        end
    end
    queue
end

function divide_coincident_intersection!(
    queue::Vector{<:SegmentEvent}, ev1::SegmentEvent, ev2::SegmentEvent, self_intersection::Bool
    ; atol::AbstractFloat=default_atol
    )
    # This assumes:
    # - ev1 is on top of or to the right of ev2, because events are processed left to right.
    # - both points of ev2 are colinear with ev1 .
    @debug("[divide_coincident_intersection!] $(ev1) -- $ev2")
    start1_on_end2 = is_same_point(ev1.segment[1], ev2.segment[2]; atol=atol)
    end1_on_start2 = is_same_point(ev1.segment[2], ev2.segment[1]; atol=atol)
    if start1_on_end2 || end1_on_start2
        return queue # segments touch at endpoints, so no further divisions
    end
    starts_equal = is_same_point(ev1.segment[1], ev2.segment[1]; atol=atol)
    ends_equal = is_same_point(ev1.segment[2], ev2.segment[2]; atol=atol)
    @debug("[divide_coincident_intersection!] starts_equal=$starts_equal")
    @debug("[divide_coincident_intersection!] ends_equal=$ends_equal")
    if starts_equal && ends_equal
        # segments are equal. Keep the second one
        return merge_same_segments!(queue, ev1, ev2, self_intersection)
    end
    start1_between = !starts_equal && on_segment(ev1.segment[1], ev2.segment; atol=atol)
    end1_between = !ends_equal && on_segment(ev1.segment[2], ev2.segment; atol=atol)
    end2_between = !ends_equal && on_segment(ev2.segment[2], ev1.segment; atol=atol)
    if starts_equal
        if end1_between
            # (a1)---(a2)
            # (b1)----x------(b2)
            divide_event!(queue, ev2, ev1.segment[2]; atol=atol)
        elseif end2_between
            # (a1)----x-----(a2)
            # (b1)---(b2)
            divide_event!(queue, ev1, ev2.segment[2]; atol=atol)
        else # are these segment colinear?
            return queue
        end
        # duplicate a1->x, so remove ev1
        return merge_same_segments!(queue, ev1, ev2, self_intersection)
    elseif start1_between
        if !ends_equal # then make a2 equal to b2
            if end1_between
                #         (a1)---(a2)
                #  (b1)-----------x-----(b2)
                divide_event!(queue, ev2, ev1.segment[2]; atol=atol)
            elseif end2_between
                #         (a1)----x-----(a2)
                #  (b1)----------(b2)
                divide_event!(queue, ev1, ev2.segment[2]; atol=atol);
            else # are these segments colinear?
                return queue
            end
        end
        #         (a1)---(a2)
        #  (b1)----x-----(b2)
        # equal segment a1-b2 isn't in the status stack yet, so don't return it
        divide_event!(queue, ev2, ev1.segment[1]; atol=atol);
    end
    queue
end

"""
    divide_event!(queue, ev, pt; atol=1e-6)

Divide an event `ev` and `ev.other` in `queue` into 4:
```
--x-->  to  --> x-->
```
"""
function divide_event!(queue::Vector{<:SegmentEvent}, ev::SegmentEvent, pt::Point2D; atol::AbstractFloat=default_atol) # eventDivide
    # assumes pt lies on ev.segment
    new_segment = (pt, ev.segment[2])
    @debug("[divide_event!] new_segment=$(new_segment)")
    e1, e2 = update_end!(ev, pt; atol=atol)
    # fix position of end in queue
    pop_key!(queue, e2)
    insert_in_order!(queue, e2; lt=compare_events)
    # add new segment at the end. Reset other_annotations
    add_segment_event!(queue, new_segment, ev.primary, ev.self_annotations, SegmentAnnotations())
end

"""
    update_end!(queue, ev, pt; atol=1e-6)

Slides an end backwards.
```
    (start)------------(end)    to:
    (start)---(end)
```
"""
function update_end!(ev::SegmentEvent, end_point::Point2D; atol::AbstractFloat=default_atol)
    # Assumes ev is a start event.        
    @assert ev.is_start
    ev.segment = (ev.segment[1], end_point)
    ev.other_point = end_point
    other = ev.other
    other.segment = (ev.segment[1], end_point)
    other.point = end_point
    if abs(ev.segment[1][1] - end_point[1]) <= atol &&
        (ev.segment[1][2] > end_point[2]) && ev.is_start
        @warn "Reversing direction for new vertical segment: $(ev.segment)."
        ev.is_start = false
        other.is_start = true
    end
    ev, other
end

function merge_same_segments!(queue::Vector{<:SegmentEvent}, discard::SegmentEvent, survive::SegmentEvent, self_intersection::Bool)
    @debug("[merge_same_segments!] discard=$discard")
    @debug("[merge_same_segments!] survive=$survive")
    pop_key!(queue, discard)
    pop_key!(queue, discard.other)
    if self_intersection
        # fill status is calculated bottom to top, so surviving's fill_below cannot change
        # however, surviving fill_above will be whatever the discarded's one would have been
        toggle = isnothing(discard.self_annotations.fill_below) ? true : discard.self_annotations.fill_above != discard.self_annotations.fill_below
        if toggle
            @assert !isnothing(survive.self_annotations.fill_above) "missing self_annotations in surviving segment: $(survive)" # preempt !nothing error
            survive.self_annotations.fill_above = !survive.self_annotations.fill_above
        end
    elseif discard.primary != survive.primary # merge two segments that belong to different polygons
        # each segment has distinct knowledge, so no special logic is needed
        # note that this can only happen once per segment in this phase, because we are guaranteed that all self-intersections are gone
        survive.other_annotations = discard.self_annotations
    else # merge two segments that belong to the same polygon
        if isnothing(survive.other_annotations.fill_above)
            @assert !isnothing(dicard.self_annotations.fill_above) "missing self_annotations in discarded segment: $(dicard)"
            survive.other_annotations = discard.other_annotations
        end
    end
    @debug("[merge_same_segments!] survive=$survive")
    queue
end

function calculate_self_annotations!(ev::SegmentEvent, below::Union{Nothing, SegmentEvent})
    # if a new segment, than toggle, else use existing knowledge
    @debug("[calculate_self_annotations!] event: $(ev)")
    @debug("[calculate_self_annotations!] below: $(below)")
    toggle = isnothing(ev.self_annotations.fill_below) ? true : ev.self_annotations.fill_above != ev.self_annotations.fill_below
    if isnothing(below)
        ev.self_annotations.fill_below = false
    else
        @assert !isnothing(below.self_annotations.fill_above) "missing annotations below: $(below)" # preempt !nothing error
        ev.self_annotations.fill_below = below.self_annotations.fill_above # below should already be filled
    end
    if toggle
        ev.self_annotations.fill_above = !ev.self_annotations.fill_below
    else
        ev.self_annotations.fill_above = ev.self_annotations.fill_below
    end
    @debug("[calculate_self_annotations!] self_annotations: $(ev.self_annotations)")
    ev.self_annotations
end

function calculate_other_annotations!(ev::SegmentEvent, below::Nothing)
    @debug("[calculate_other_annotations!] ev=$ev")
    @debug("[calculate_other_annotations!] below=$below")
    if isnothing(ev.other_annotations.fill_above)
        # if nothing is below the event, it cannot be in the other polygon
        inside = false
        ev.other_annotations.fill_above = inside
        ev.other_annotations.fill_below = inside
    end
    @debug("[calculate_other_annotations!] ev.other_annotations=$(ev.other_annotations)")
    ev.other_annotations
end
    
function calculate_other_annotations!(ev::SegmentEvent, below::SegmentEvent)
    @debug("[calculate_other_annotations!] ev=$ev")
    @debug("[calculate_other_annotations!] below=$below")
    if isnothing(ev.other_annotations.fill_above)
        # something is below ev, so copy the below segment's other polygon's above
        inside = (ev.primary == below.primary) ? 
            below.other_annotations.fill_above : below.self_annotations.fill_above
        ev.other_annotations.fill_above = inside
        ev.other_annotations.fill_below = inside
    end
    @debug("[calculate_other_annotations!] ev.other_annotations=$(ev.other_annotations)")
    ev.other_annotations
end

#############################################################
##                 Segment Selection                       ##
#############################################################

#=
Table
    Above1  Below1 Above2 Below2
01. No      No     No     No
02. No      No     No     Yes
03. No      No     Yes    No
04. No      No     Yes    Yes
05. No      Yes    No     No
06. No      Yes    No     Yes
07. No      Yes    Yes    No
08. No      Yes    Yes    Yes
09. Yes     No     No     No
10. Yes     No     No     Yes
11. Yes     No     Yes    No
12. Yes     No     Yes    Yes
13. Yes     Yes    No     No
14. Yes     Yes    No     Yes
15. Yes     Yes    Yes    No
16. Yes     Yes    Yes    Yes
=#

# For the 16 rows, indicate if (1) in the output shape (non-blank) and (2) which side is filled in that output shape

INTERSECTION_CRITERIA = [
    BLANK, BLANK, BLANK, BLANK,
    BLANK, BELOW, EMPTY, BELOW,
    BLANK, EMPTY, ABOVE, ABOVE,
    BLANK, BELOW, ABOVE, BLANK,
] # both below, both above, but not all 4

INTERSECTION_SEGMENT_CRITERIA = [
    BLANK, BLANK, BLANK, BLANK,
    BLANK, BLANK, EMPTY, BLANK,
    BLANK, EMPTY, BLANK, BLANK,
    BLANK, BLANK, BLANK, BLANK,
] # above and below different segments

UNION_CRITERIA = [
    BLANK, BELOW, ABOVE, BLANK,
    BELOW, BELOW, BLANK, BLANK,
    ABOVE, BLANK, ABOVE, BLANK,
    BLANK, BLANK, BLANK, BLANK,
] # filled only above/only below. (above1 | above2) ⊻ (below1 | below2)

DIFFERENCE_CRITERIA = [
    BLANK, BLANK, BLANK, BLANK,
    BELOW, BLANK, BELOW, BLANK,
    ABOVE, ABOVE, BLANK, BLANK,
    BLANK, ABOVE, BELOW, BLANK,
] # primary - secondary. (above1 && !above2) ⊻  (below1 && !below2)

XOR_CRITERIA = [
    BLANK, BELOW, ABOVE, BLANK,
    BELOW, BLANK, BLANK, ABOVE,
    ABOVE, BLANK, BLANK, BELOW,
    BLANK, ABOVE, BELOW, BLANK,
] # (above1 ⊻ above2) ⊻ (below1 ⊻ below2)

function apply_selection_criteria(annotated_segments::Vector{<:SegmentEvent{T}}, criteria::Vector{AnnotationFill}) where T
    result = SegmentEvent{T}[]
    for ev in annotated_segments
        index = (ev.self_annotations.fill_above ? 9 : 1) +
                (ev.self_annotations.fill_below ? 4 : 0) + 
                (ev.other_annotations.fill_above ? 2 : 0) + 
                (ev.other_annotations.fill_below ? 1 : 0)
        if criteria[index] != BLANK
            new_segment = copy_segment(ev, true)
            new_segment.self_annotations.fill_above = criteria[index] == ABOVE
            new_segment.self_annotations.fill_below = criteria[index] == BELOW
            new_segment.other_annotations.fill_above = nothing
            new_segment.other_annotations.fill_below = nothing
            push!(result, new_segment)
        end
    end
    result
end

#############################################################
##                  Polygons and Holes                     ##
#############################################################

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

"""
    match_holes_polygons(polygons::Vector{<:Polygon}, holes::Vector{<:Tuple})

An algorithm for matching holes to polygons.
Returns the index of each parent for each hole.

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
        found = false
        for (idx_p, parent) in zip(idxs, polygons[idxs])
            # Assume that no segments intersect after the Martínez-Rueda algorithm.
            # Then only need to check a point not on the exterior.
            j = 1
            while (j < length(candidate)) && on_border(parent.exterior, candidate[j])
                j += 1
            end
            found = contains(parent.exterior, candidate[j]; atol=atol, on_border_is_inside=false)
            if found
                parents[idx_h] = idx_p
                break
            end
        end
    end
    parents
end
