"""
    AnnotationFill

Outcome for a segment in a selection criteria table. See `apply_selection_criteria`.
- `BLANK`: the segment is not part of the output.
- `ABOVE`: the segment is part of the output and the output is filled above it.
- `BELOW`: the segment is part of the output and the output is filled below it.
- `EMPTY`: the segment is part of the output but neither side is filled e.g. a shared edge
    between two polygons that only touch or an extended edge.
"""
@enum AnnotationFill BLANK=0 ABOVE=1 BELOW=2 EMPTY=3

"""
    martinez_rueda_algorithm(
    selection_criteria, subject, others...
    ; atol=default_atol, rtol=default_rtol, face_selection=SPLIT_FACES, fill_rule=EVEN_ODD
    )

The Martínez-Rueda-Feito polygon clipping algorithm.
Returns regions and edges of intersection.
It runs in `O((n+m+k)log(n+m))` time where `n` and `m` are the number of vertices of `polygon1` 
and `polygon2` respectively and `k` is the total number of intersections between all segments.
Use `intersect_convex` for convex polygons for an `O(n+m)` algorithm.

See `PolygonAlgorithms.segments_to_paths` for details on the `face_selection` parameter.

The input polygons can be:
- A list of points: `Vector{Tuple{Float64, Float64}}`.
- Polygons: `Vector{Polygon{T}}`.
- A Segment event queue: `Vector{<:SegmentEvent{Float64}}`. The core algorithm uses this representation.
- The `subject` can also be a list of polygons. In this case, all polygons are treated as a single polygon.
    If any of the polygons overlap, this is like passing a self-intersecting polygon and
    some areas might be classified as holes according to the `fill_rule`..

Keyword arguments:
- `atol`, `rtol`: tolerances for point comparisons, orientations and intersections.
- `fill_rule`: rule for which regions of self-intersecting polygons are filled. See `FillRule`.
- `face_selection`: how to convert the resultant segments back into polygons. See `FaceSelectionStrategy`.
    Not used for `SegmentEvent` inputs.

`selection_criteria` determines the boolean operation. It is one of `INTERSECTION_CRITERIA`, `UNION_CRITERIA`,
`DIFFERENCE_CRITERIA` or `XOR_CRITERIA`. See `apply_selection_criteria`.

Description:
- Operates at a segment level and is an extension of the Bentley-Ottmann line intersection algorithm.
    Segments are scanned from left to right, bottom to top. 
- The key assumption is that only the segments immediately above and below the current segment need to be inspected for intersections.
    This makes the algorithm fast but also sensitive to determining these segments correctly.
- The segment that is immediately below (or empty space) is used to determine the fill annotations for the current segment.
- Once all annotations are done, the desired segments can be selected based on a given criteria.
- These segments are then chained back to paths and polygons.
    This requires casting to a grid to match starting and end points of segments.
    This is achieved by rounding any decimal places, by default to the 6th decimal place.
    See `PolygonAlgorithms.segments_to_paths` and `PolygonAlgorithms.segments_to_polygons` for more detail. 

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
    ;
    atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_rtol,
    face_selection::FaceSelectionStrategy=SPLIT_FACES,
    options...
    ) where T
    subject_events = convert_to_events(subject; primary=true, atol=atol)
    other_events = map(
        p -> convert_to_events(p; primary=false, atol=atol)
        , others
    )
    segments = martinez_rueda_algorithm(
        selection_criteria, subject_events, other_events...
        ; atol=atol, rtol=rtol, options...
    )
    exteriors, holes = segments_to_paths(
        segments;
        atol=atol,
        face_selection=face_selection
    )
    vcat(exteriors, holes)
end

# Multiple subjects
function martinez_rueda_algorithm(
    selection_criteria::Vector{AnnotationFill},
    subjects::AbstractVector{<:Path2D{T}},
    others::Vararg{Path2D{T}},
    ;
    atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_rtol,
    face_selection::FaceSelectionStrategy=SPLIT_FACES,
    options...
    ) where T
    subject_events = SegmentEvent{T}[]
    for p in subjects
        push!(
            subject_events,
            convert_to_events(p; primary=true, atol=atol)...
        )
    end
    other_events = map(
        p -> convert_to_events(p; primary=false, atol=atol),
        others
    )
    segments = martinez_rueda_algorithm(
        selection_criteria, subject_events, other_events...
        ; atol=atol, rtol=rtol, options...
    )
    exteriors, holes = segments_to_paths(
        segments;
        atol=atol,
        face_selection=face_selection
    )
    vcat(exteriors, holes)
end

# Polygons with holes
function martinez_rueda_algorithm(
    selection_criteria::Vector{AnnotationFill},
    subject::Polygon{T},
    others::Vararg{Polygon{T}},
    ;
    atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_rtol,
    face_selection::FaceSelectionStrategy=SPLIT_FACES,
    options...
    ) where T
    subject_events = convert_to_events(subject.exterior; atol=atol)
    for p in subject.holes
        push!(
            subject_events,
            convert_to_events(p; primary=true, atol=atol)...
        )
    end
    other_events = map(
        p -> convert_to_events(p.exterior; primary=false, atol=atol),
        others
    )
    for (other_event, other) in zip(other_events, others)
        for p in other.holes
            push!(other_event, convert_to_events(p; primary=false, atol=atol)...)
        end
    end
    segments = martinez_rueda_algorithm(
        selection_criteria, subject_events, other_events...
        ; atol=atol, rtol=rtol, options...
    )
    segments_to_polygons(
        segments
        ; atol=atol, rtol=rtol,
        face_selection=face_selection,
    )
end

# Multiple subjects with holes
function martinez_rueda_algorithm(
    selection_criteria::Vector{AnnotationFill},
    subjects::AbstractVector{<:Polygon{T}},
    others::Vararg{Polygon{T}},
    ;
    atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_rtol,
    face_selection::FaceSelectionStrategy=SPLIT_FACES,
    options...
    ) where T
    subject_events = SegmentEvent{T}[]
    for subject in subjects
        for p in vcat([subject.exterior], subject.holes)
            push!(subject_events, convert_to_events(p; primary=true, atol=atol)...)
        end
    end
    other_events = map(
        p -> convert_to_events(p.exterior; primary=false, atol=atol),
        others
    )
    for (other_event, other) in zip(other_events, others)
        for p in other.holes
            push!(other_event, convert_to_events(p; primary=false, atol=atol)...)
        end
    end
    segments = martinez_rueda_algorithm(
        selection_criteria, subject_events, other_events...
        ; atol=atol, rtol=rtol, options...
    )
    segments_to_polygons(
        segments
        ; atol=atol, rtol=rtol,
        face_selection=face_selection,
    )
end

# Core algorithm: SegmentEvent input
function martinez_rueda_algorithm(
    selection_criteria::Vector{AnnotationFill},
    subject::Vector{<:SegmentEvent{T}},
    polygons::Vararg{Vector{<:SegmentEvent{T}}},
    ; atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_rtol,
    fill_rule::FillRule=EVEN_ODD
    ) where T
    base_annotated_segments = event_loop!(subject; self_intersection=true, atol=atol, rtol=rtol, fill_rule=fill_rule)
    for polygon in polygons
        annotated_segments = event_loop!(polygon; self_intersection=true, atol=atol, rtol=rtol, fill_rule=fill_rule)
        queue = SegmentEvent{T}[]
        for ev in vcat(base_annotated_segments, annotated_segments)
            push!(queue, create_segment_pair(ev)...)
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

# Core algorithm for 1 polygon
function martinez_rueda_algorithm(
    selection_criteria::Vector{AnnotationFill},
    subject::Vector{<:SegmentEvent{T}},
    ; atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_rtol,
    fill_rule::FillRule=EVEN_ODD
    ) where T
    base_annotated_segments = event_loop!(subject; self_intersection=true, atol=atol, rtol=rtol, fill_rule=fill_rule)
    for seg in base_annotated_segments
        # apply_selection_criteria will fail (and should fail) if annotations are nothing
        seg.other_annotations.fill_above = false
        seg.other_annotations.fill_below = false
    end
    selected_segments = apply_selection_criteria(base_annotated_segments, selection_criteria)
    selected_segments
end

#############################################################
##                    Event Loop                           ##
#############################################################

"""
    event_loop!(queue; self_intersection, atol=default_atol, rtol=default_rtol, fill_rule=EVEN_ODD)

The core line sweep of the Martinez-Rueda algorithm.
Consumes `queue` and returns a list of annotated segments (copies of the start events).

Events are processed in order from left to right.
For each start event, the segments immediately above and below it in the sweep status are checked for intersections.
If one is found, the segments are divided and the event is processed again later.
Otherwise the event's annotations are calculated from the segment below it and it is inserted into the sweep status.
For each end event, the segment is removed from the sweep status and the two newly adjacent segments
are checked for intersections.

The loop is used in two phases:
1. `self_intersection=true`: process segments from a single polygon.
    Calculates `self_annotations` based on the `fill_rule`.
2. `self_intersection=false`: process segments from two polygons, which have already been through phase 1.
"""
function event_loop!(
    queue::Vector{SegmentEvent{T}}
    ; self_intersection::Bool,
    atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_rtol,
    fill_rule::FillRule=EVEN_ODD
    ) where T # eventLoop
    annotated_segments = SegmentEvent{T}[]
    sweep_status = SegmentEvent{T}[] # current events in a vertical line, top to bottom.
    sort_event_queue!(queue; atol=atol, rtol=rtol)
    while !(isempty(queue))
        head = queue[1]
        queue_length = length(queue)
        status_length = length(sweep_status)
        @debug("[event_loop!] ($(queue_length), $(status_length)): $(head)")
        if head.is_start # then check for intersections and add to sweep status
            idx = find_transition(sweep_status, head; atol=atol, rtol=rtol)
            above = idx == 1 ? nothing : sweep_status[idx - 1]
            below = (idx > length(sweep_status)) ? nothing : sweep_status[idx]
            @debug("[event_loop!] transition idx=$idx")
            @debug("[event_loop!] above=$above")
            @debug("[event_loop!] below=$below")
            check_and_divide_intersection!(
                queue, head, above, self_intersection
                ; atol=atol, rtol=rtol, fill_rule=fill_rule
            )
            if queue[1] !== head
                continue # either head was removed or something was inserted ahead of it
            end
            check_and_divide_intersection!(
                queue, head, below, self_intersection
                ; atol=atol, rtol=rtol, fill_rule=fill_rule
            )
            if queue[1] !== head
                continue # either head was removed or something was inserted ahead of it
            end
            if self_intersection
                calculate_self_annotations!(head, below; fill_rule=fill_rule)
            else
                calculate_other_annotations!(head, below)
            end
            insert!(sweep_status, idx, head)
        else # event is ending, so remove it from the status
            idx = find_transition(sweep_status, head.other; atol=atol, rtol=rtol)
            if !(0 < idx <= length(sweep_status) && sweep_status[idx] === head.other)
                @warn "$(head.other) was not in the expected location in the sweep status. " * 
                    "Falling back to linear search."
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
                    ; atol=atol, rtol=rtol, fill_rule=fill_rule
                )
            end
            push!(annotated_segments, copy(head.other))
            popat!(sweep_status, idx)
        end
        popfirst!(queue)
    end
    annotated_segments
end

function check_and_divide_intersection!(
    queue::Vector{<:SegmentEvent}, ev1::SegmentEvent, ev2::Nothing, self_intersection::Bool
    ; atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_rtol, fill_rule::FillRule=EVEN_ODD
    )
    queue
end

"""
    check_and_divide_intersection!(
        queue, ev1, ev2, self_intersection
        ; atol=default_atol, rtol=default_rtol, fill_rule=EVEN_ODD
    )

Check if the segments of start events `ev1` and `ev2` intersect.
If they intersect at a single point, divide them at that point with `divide_intersection!`.
If they are colinear and overlap, divide and merge them with `divide_coincident_intersection!`.
If `ev2` is `nothing`, do nothing.

Returns `queue`, which is modified in place.
"""
function check_and_divide_intersection!(
    queue::Vector{<:SegmentEvent},
    ev1::SegmentEvent,
    ev2::SegmentEvent,
    self_intersection::Bool
    ; atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_rtol,
    fill_rule::FillRule=EVEN_ODD
    )
    pt = intersect_geometry(ev1.segment, ev2.segment; atol=atol, rtol=rtol)
    if isnothing(pt)
        @debug("no intersection or parallel lines at $ev1 -- $ev2")
        # Lines might be on top of each other 
        ori2_start = get_orientation(ev1.segment[1], ev1.segment[2], ev2.segment[1]; atol=atol, rtol=rtol)
        ori2_end = get_orientation(ev1.segment[1], ev1.segment[2], ev2.segment[2]; atol=atol, rtol=rtol)
        if (ori2_start == COLINEAR) && (ori2_end == COLINEAR)
            divide_coincident_intersection!(
                queue, ev1, ev2, self_intersection
                ; fill_rule=fill_rule, atol=atol, rtol=rtol
            )
        end
        return queue
    else
        divide_intersection!(queue, ev1, ev2, pt, self_intersection; atol=atol, rtol=rtol)
    end
end

function divide_intersection!(
    queue::Vector{<:SegmentEvent}, ev1::SegmentEvent, ev2::SegmentEvent, pt, self_intersection::Bool=false
    ; atol=default_atol, rtol=default_rtol
    )
    queue
end

"""
    divide_intersection!(queue, ev1, ev2, pt, [self_intersection=false]; atol=default_atol, rtol=default_rtol)

Divide the segments of `ev1` and `ev2` at their intersection point `pt` with `divide_event!`.

A segment is only divided if `pt` lies strictly between its end points.
If `pt` is near an end point of the other segment, that end point is used instead of `pt`
so that the divided segments share exactly the same vertex.
Nothing is done if `pt` is `nothing`.
"""
function divide_intersection!(
    queue::Vector{<:SegmentEvent},
    ev1::SegmentEvent,
    ev2::SegmentEvent,
    pt::Point2D,
    self_intersection::Bool=false,
    ; atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_rtol
    ) # checkIntersection
    _divide_event!(args...) = divide_event!(args...; atol=atol, rtol=rtol)
    @debug("[divide_intersection!] $(ev1.segment) -- $(ev2.segment) at $(pt)")
    at_start1, at_end1, along1 = classify_intersection(ev1.segment, pt; atol=atol)
    at_start2, at_end2, along2 = classify_intersection(ev2.segment, pt; atol=atol)
    @debug("[divide_intersection!] $at_start1 $at_end1 $along1")
    @debug("[divide_intersection!] $at_start2 $at_end2 $along2")
    if along1 && along2
        _divide_event!(queue, ev1, pt, self_intersection)
        _divide_event!(queue, ev2, pt, self_intersection)
    elseif along1
        if at_start2
            _divide_event!(queue, ev1, ev2.segment[1], self_intersection)
        elseif at_end2
            _divide_event!(queue, ev1, ev2.segment[2], self_intersection)
        end
    elseif along2
        if at_start1
            _divide_event!(queue, ev2, ev1.segment[1], self_intersection)
        elseif at_end1
            _divide_event!(queue, ev2, ev1.segment[2], self_intersection)
        end
    end
    queue
end

"""
    divide_coincident_intersection!(
        queue, ev1, ev2, self_intersection
        ; fill_rule=EVEN_ODD, atol=default_atol, rtol=default_rtol
    )

Handle two colinear segments that overlap.
The segments are divided so that the overlapping part is represented by two identical segments,
which are then merged into one with `merge_same_segments!`.

Assumes:
- `ev1` is on top of or to the right of `ev2`, because events are processed left to right.
- Both points of `ev2` are colinear with `ev1`.

If the segments only touch at an end point, nothing is done.
"""
function divide_coincident_intersection!(
    queue::Vector{<:SegmentEvent}, ev1::SegmentEvent, ev2::SegmentEvent, self_intersection::Bool
    ; fill_rule::FillRule=EVEN_ODD, atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_rtol
    )
    # This assumes:
    # - ev1 is on top of or to the right of ev2, because events are processed left to right.
    # - both points of ev2 are colinear with ev1
    _divide_event!(args...) = divide_event!(args...; atol=atol, rtol=rtol)
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
        return merge_same_segments!(queue, ev1, ev2, self_intersection; fill_rule=fill_rule)
    end
    start1_between = !starts_equal && on_segment(ev1.segment[1], ev2.segment; atol=atol)
    end1_between = !ends_equal && on_segment(ev1.segment[2], ev2.segment; atol=atol, rtol=rtol)
    end2_between = !ends_equal && on_segment(ev2.segment[2], ev1.segment; atol=atol, rtol=rtol)
    if starts_equal
        if end1_between
            # (a1)---(a2)
            # (b1)----x------(b2)
            _divide_event!(queue, ev2, ev1.segment[2], self_intersection)
        elseif end2_between
            # (a1)----x-----(a2)
            # (b1)---(b2)
            _divide_event!(queue, ev1, ev2.segment[2], self_intersection)
        else # are these segment colinear?
            return queue
        end
        # duplicate a1->x, so remove ev1
        return merge_same_segments!(queue, ev1, ev2, self_intersection; fill_rule=fill_rule)
    elseif start1_between
        if !ends_equal # then make a2 equal to b2
            if end1_between
                #         (a1)---(a2)
                #  (b1)-----------x-----(b2)
                _divide_event!(queue, ev2, ev1.segment[2], self_intersection)
            elseif end2_between
                #         (a1)----x-----(a2)
                #  (b1)----------(b2)
                _divide_event!(queue, ev1, ev2.segment[2], self_intersection);
            else # are these segments colinear?
                return queue
            end
        end
        #         (a1)---(a2)
        #  (b1)----x-----(b2)
        # equal segment a1-b2 isn't in the status stack yet, so don't return it
        _divide_event!(queue, ev2, ev1.segment[1], self_intersection);
    end
    queue
end

"""
    divide_event!(queue, ev, pt, [self_intersection=false]; atol=default_atol, rtol=default_rtol)

Divide a start event `ev` and its end event `ev.other` in `queue` into 4 events at the point `pt`:
```
--x-->  to  --> x-->
```

The first segment keeps `ev` with its end moved to `pt` (see `update_end!`) and its end event is re-sorted in the queue.
A new pair of events is added for the second segment from `pt` to the original end.

The new segment's `other_annotations` are reset.
If `self_intersection` is true, its fill annotations are also reset because they will be recalculated,
but its windings are kept because `ev` may represent multiple merged coincident edges.
Otherwise it shares `ev`'s `self_annotations`.

Assumes `pt` lies on `ev.segment`.
"""
function divide_event!(
    queue::Vector{<:SegmentEvent}, ev::SegmentEvent, pt::Point2D,
    self_intersection::Bool=false
    ; atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_rtol
    ) # eventDivide
    # assumes pt lies on ev.segment
    lt(a, b) = compare_events(a,b; atol=atol, rtol=rtol)
    new_segment = (pt, ev.segment[2])
    @debug("[divide_event!] new_segment=$(new_segment)")
    e1, e2 = update_end!(ev, pt)
    # fix position of end in queue
    pop_key!(queue, e2)
    insert_in_order!(queue, e2; lt=lt)
    # add new segment at the end. Reset other_annotations. Reset self-annotations if self-intersection.
    if self_intersection
        # the fill can be recalculated, but the windings must be kept because ev might
        # represent multiple merged coincident edges. See merge_same_segments!
        self_annotations = SegmentAnnotations(
            forward=ev.self_annotations.forward,
            winding_top_to_bottom=ev.self_annotations.winding_top_to_bottom,
            winding_left_to_right=ev.self_annotations.winding_left_to_right,
        )
    else
        self_annotations = ev.self_annotations
    end
    start_event, end_event = create_segment_pair(
        new_segment, ev.primary,
        self_annotations,
        SegmentAnnotations()
    )
    insert_in_order!(queue, start_event; lt=lt)
    insert_in_order!(queue, end_event; lt=lt)
end

"""
    update_end!(ev, end_point)

Slides the end of the segment of start event `ev` backwards to `end_point`,
updating both `ev` and `ev.other`. Returns `(ev, ev.other)`.

The position of `ev.other` in the event queue is not updated.
```
    (start)------------(end)    to:
    (start)---(end)
```
"""
function update_end!(ev::SegmentEvent, end_point::Point2D)
    # Assumes ev is a start event.        
    @assert ev.is_start
    ev.segment = (ev.segment[1], end_point)
    ev.other_point = end_point
    other = ev.other
    other.segment = (ev.segment[1], end_point)
    other.point = end_point
    ev, other
end

"""
    merge_same_segments!(queue, discard, survive, self_intersection; fill_rule=EVEN_ODD)

Merge two identical segments by removing `discard` (and `discard.other`) from the queue
and combining its annotations into `survive`.
Assumes `discard` is on top of `survive`.

There are three cases:
1. `self_intersection`: the segments belong to the same polygon. The windings of `discard` are added
    to `survive` and `survive.fill_above` is recalculated, so that `survive` represents all the coincident edges.
2. The segments belong to different polygons: `discard`'s self annotations become `survive`'s other annotations.
3. The segments belong to the same polygon in the combining phase: `survive` takes `discard`'s other annotations
    if it does not have any.
"""
function merge_same_segments!(
    queue::Vector{<:SegmentEvent}, discard::SegmentEvent, survive::SegmentEvent, self_intersection::Bool
    ; fill_rule::FillRule=EVEN_ODD
    )
    # discard is assumed to be on top of survive
    @debug("[merge_same_segments!] discard=$discard")
    @debug("[merge_same_segments!] survive=$survive")
    pop_key!(queue, discard)
    pop_key!(queue, discard.other)
    if self_intersection
        # fill status is calculated bottom to top, so survive's winding_below and fill_below are already correct.
        # Add discard's windings to survive's, so that survive's winding change accounts for all its coincident edges.
        @assert !isnothing(survive.self_annotations.fill_above) "missing self_annotations in surviving segment: $(survive)"
        annotations = survive.self_annotations
        annotations.winding_top_to_bottom = get_winding_top_to_bottom!(survive) + get_winding_top_to_bottom!(discard)
        annotations.winding_left_to_right = get_winding_left_to_right!(survive) + get_winding_left_to_right!(discard)
        if fill_rule == EVEN_ODD
            is_odd_count = isodd(get_winding_change!(survive))
            annotations.fill_above = is_odd_count ? !annotations.fill_below : annotations.fill_below
        else
            winding_above = annotations.winding_below + get_winding_change!(survive)
            annotations.fill_above = is_filled(winding_above, fill_rule)
        end
    elseif discard.primary != survive.primary # merge two segments that belong to different polygons
        # each segment has distinct knowledge, so no special logic is needed
        # note that this can only happen once per segment in this phase, because we are guaranteed that all self-intersections are gone
        survive.other_annotations = discard.self_annotations
    else # merge two segments that belong to the same polygon
        if isnothing(survive.other_annotations.fill_above)
            @assert !isnothing(discard.self_annotations.fill_above) "missing self_annotations in discarded segment: $(discard)"
            survive.other_annotations = discard.other_annotations
        end
    end
    @debug("[merge_same_segments!] survive=$survive\n")
    queue
end

"""
    calculate_self_annotations!(ev, below; fill_rule)

Calculate `fill_above` and `fill_below` of `ev.self_annotations` from the segment immediately `below` it
in the sweep status, or from empty space if `below` is `nothing`.

- `EVEN_ODD`: `fill_below` is the `fill_above` of the segment below.
    `fill_above` is toggled if `ev` represents an odd number of edges.
- Winding rules (`NON_ZERO`, `POSITIVE`, `NEGATIVE`): the winding number below is the winding number above the segment below.
    The winding number above is that plus the change in winding across `ev`.
    Both are converted to fills with `is_filled`.
"""
function calculate_self_annotations!(
    ev::SegmentEvent,
    below::Union{Nothing, SegmentEvent}
    ;  fill_rule::FillRule
    )
    @debug("[calculate_self_annotations!] event: $(ev)")
    @debug("[calculate_self_annotations!] below: $(below)")
    if fill_rule == EVEN_ODD
        if isnothing(below)
            ev.self_annotations.fill_below = false
        else
            @assert !isnothing(below.self_annotations.fill_above) "missing annotations below: $(below)" # preempt !nothing error
            ev.self_annotations.fill_below = below.self_annotations.fill_above # below should already be filled
        end
        # Each edge toggles the fill. Usually ev is one edge, but it may represent multiple coincident edges
        # (see merge_same_segments!). Its winding change is the sum of their ±1 changes, so it is odd if
        # and only if the number of edges is odd.
        is_odd_count = isodd(get_winding_change!(ev))
        ev.self_annotations.fill_above = is_odd_count ? !ev.self_annotations.fill_below : ev.self_annotations.fill_below
    else
        # the winding from the region below to the bottom
        #   = winding of the region above the below segment
        #   = (winding_below + winding_top_to_bottom) of the below segment
        winding_below = isnothing(below) ? 0 :
            below.self_annotations.winding_below + get_winding_top_to_bottom!(below)
        ev.self_annotations.winding_below = winding_below
        # the winding from the region above to the bottom = winding_below + Δwinding
        winding_above = winding_below + get_winding_change!(ev)
        ev.self_annotations.fill_above = is_filled(winding_above, fill_rule)
        ev.self_annotations.fill_below = is_filled(winding_below, fill_rule)
    end
    @debug("[calculate_self_annotations!] self_annotations: $(ev.self_annotations)")
    ev.self_annotations
end

"""
    calculate_other_annotations!(ev, below)

Calculate whether `ev` is inside the other polygon, based on the segment immediately `below` it.
Because segments have already been divided, the other polygon is on both sides of `ev` or neither,
so `fill_above == fill_below`.

If `below` is `nothing`, `ev` is outside the other polygon.
If `below` is from the same polygon, copy its `other_annotations.fill_above`;
otherwise copy its `self_annotations.fill_above`.

Existing annotations (e.g. from merging coincident segments) are not overwritten.
"""
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

"""
    INTERSECTION_CRITERIA

Selection criteria for the intersection: regions filled by both polygons.
"""
INTERSECTION_CRITERIA = [
    BLANK, BLANK, BLANK, BLANK,
    BLANK, BELOW, EMPTY, BELOW,
    BLANK, EMPTY, ABOVE, ABOVE,
    BLANK, BELOW, ABOVE, BLANK,
] # both below, both above, but not all 4

"""
    INTERSECTION_SEGMENT_CRITERIA

Selection criteria for only the zero-area parts of the intersection..
The polygons may still overlap elsewhere.
"""
INTERSECTION_SEGMENT_CRITERIA = [
    BLANK, BLANK, BLANK, BLANK,
    BLANK, BLANK, EMPTY, BLANK,
    BLANK, EMPTY, BLANK, BLANK,
    BLANK, BLANK, BLANK, BLANK,
]

"""
    UNION_CRITERIA

Selection criteria for the union: regions filled by either polygon.
"""
UNION_CRITERIA = [
    BLANK, BELOW, ABOVE, BLANK,
    BELOW, BELOW, BLANK, BLANK,
    ABOVE, BLANK, ABOVE, BLANK,
    BLANK, BLANK, BLANK, BLANK,
] # filled only above/only below. (above1 | above2) ⊻ (below1 | below2)

"""
    DIFFERENCE_CRITERIA

Selection criteria for the difference: regions filled by the primary polygon but not the secondary polygon.
"""
DIFFERENCE_CRITERIA = [
    BLANK, BLANK, BLANK, BLANK,
    BELOW, BLANK, BELOW, BLANK,
    ABOVE, ABOVE, BLANK, BLANK,
    BLANK, ABOVE, BELOW, BLANK,
] # primary - secondary. (above1 && !above2) ⊻  (below1 && !below2)

"""
    XOR_CRITERIA

Selection criteria for the exclusive or: regions filled by exactly one of the polygons.
"""
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
            new_segment = SegmentEvent(
                ev.segment,
                ev.is_start,
                true,
                SegmentAnnotations(
                    criteria[index] == ABOVE,
                    criteria[index] == BELOW,
                ),
                SegmentAnnotations(),
            )
            push!(result, new_segment)
        end
    end
    result
end

