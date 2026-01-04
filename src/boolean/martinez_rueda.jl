@enum AnnotationFill BLANK=0 ABOVE=1 BELOW=2 EMPTY=3

"""
    martinez_rueda_algorithm(selection_criteria, base, others...; atol=default_atol)

The Martínez-Rueda-Feito polygon clipping algorithm.
Returns regions and edges of intersection.
It runs in `O((n+m+k)log(n+m))` time where `n` and `m` are the number of vertices of `polygon1` 
and `polygon2` respectively and `k` is the total number of intersections between all segments.
Use `intersect_convex` for convex polygons for an `O(n+m)` algorithm.

The `base` and `others` must both be either:
- Polygons: `Vector{Tuple{Float64, Float64}}`.
- Multi-polygons: `Vector{Vector{Tuple{Float64, Float64}}}`. Note: these are not interpreted as holes.
- Segment event queue: `Vector{<:SegmentEvent{Float64}}`. The core algorithm uses this representation.

Description:
- Operates at a segment level and is an extension of the Bentley-Ottman line intersection algorithm.
Segments are scanned from left to right, bottom to top. 
- The key assumption is that only the segments immediately above and below the current segment need to be inspected for intersections.
This makes the algorithm fast but also sensitive to determining these segments correctly.
- The segment that is immediately below (or empty space) is used to determine the fill annotations for the current segment.
- Once all annotations are done, the desired segments can be selected that match a given criteria.

Limitations
1. It can fail for improper polygons: polygons with lines sticking out.
2. It is sensitive to numeric inaccuracies e.g. a line that is almost vertical or tiny regions 
of intersection.

References 
- paper: https://www.researchgate.net/publication/220163820_A_new_algorithm_for_computing_Boolean_operations_on_polygons
- article: https://sean.fun/a/polygon-clipping-pt2/
- article source code: https://github.com/velipso/polybooljs
"""
function martinez_rueda_algorithm(
    selection_criteria::Vector{AnnotationFill},
    base::Path2D{T},
    others::Vararg{Path2D{T}},
    ; atol::AbstractFloat=default_atol
    ) where T
    event_queue_base = convert_to_event_queue(base; primary=true, atol=atol)
    event_queue_others = map(p -> convert_to_event_queue(p; primary=false, atol=atol), others)
    martinez_rueda_algorithm(selection_criteria, event_queue_base, event_queue_others...; atol=atol)
end

function martinez_rueda_algorithm(
    selection_criteria::Vector{AnnotationFill},
    base::Polygon{T},
    others::Vararg{Polygon{T}},
    ; atol::AbstractFloat=default_atol
    ) where T
    event_queue_base = convert_to_event_queue(base.exterior; primary=true, atol=atol)
    for hole in base.holes
        convert_to_event_queue!(event_queue_base, hole; primary=true, atol=atol)
    end
    event_queue_others = map(p -> convert_to_event_queue(p.exterior; primary=false, atol=atol), others)
    for (queue, other) in zip(event_queue_others, others)
        for hole in other.holes
            convert_to_event_queue!(queue, hole; primary=false, atol=atol)
        end
    end
    #TODO return Vector{Polygon}. This will require:
    # 1. Return Vector{Vector{SegmentEvent}}. Then convert to Path2D/Poylgon as required.
    # 2. An indicator if a polygon is a hole or not.
    # 3. Sorting of the resulting regions into holes/polygons.
    #    Caveat: Handle polygons inside the hole of another polygon.
    martinez_rueda_algorithm(selection_criteria, event_queue_base, event_queue_others...; atol=atol)
end

function martinez_rueda_algorithm(
    selection_criteria::Vector{AnnotationFill},
    subjects::AbstractVector{<:Path2D{T}},
    clips::AbstractVector{<:Path2D{T}},
    ; atol::AbstractFloat=default_atol
    ) where T
    event_queue_subjects = map(p -> convert_to_event_queue(p; primary=true, atol=atol), subjects)
    event_queue_clips = map(p -> convert_to_event_queue(p; primary=false, atol=atol), clips)
    martinez_rueda_algorithm(selection_criteria, event_queue_subjects..., event_queue_clips...; atol=atol)
end

function martinez_rueda_algorithm(
    selection_criteria::Vector{AnnotationFill},
    base::Vector{<:SegmentEvent{T}},
    polygons::Vararg{Vector{<:SegmentEvent{T}}},
    ; atol::AbstractFloat=default_atol
    ) where T
    base_annotated_segments = event_loop!(base; self_intersection=true, atol=atol)
    for polygon in polygons
        annotated_segments = event_loop!(polygon; self_intersection=true, atol=atol)
        queue = SegmentEvent{T}[]
        for ev in vcat(base_annotated_segments, annotated_segments)
            add_annotated_segment!(queue, ev)
        end
        combined_annotated_segments = event_loop!(queue; self_intersection=false, atol=atol)
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
    empty_segments, regions_segments = separate(is_empty_segment, base_annotated_segments)
    regions = chain_segments(regions_segments; atol=atol, check_closes=true)
    segment_chains = chain_segments(empty_segments; atol=atol, check_closes=false)
    # It is also possible to attach some segment_chains to regions.
    # This will give consistent results with the Weiler-Atherton implementation.
    # For now, skipping this step.
    vcat(regions, segment_chains)
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
    ; self_intersection::Bool, atol::AbstractFloat=default_atol
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
            check_and_divide_intersection!(queue, head, above, self_intersection; atol=atol)
            if queue[1] != head
                continue # either head was removed or something was inserted ahead of it
            end
            check_and_divide_intersection!(queue, head, below, self_intersection; atol=atol)
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
                    ; atol=atol)
            end
            push!(annotated_segments, copy_segment(head.other, head.other.primary))
            popat!(sweep_status, idx)
        end
        popfirst!(queue)
    end
    annotated_segments
end

function check_and_divide_intersection!(
    queue::Vector{<:SegmentEvent}, ev1::SegmentEvent, ev2::Nothing, self_intersection::Bool; atol=1e-6
    )
    queue
end

function check_and_divide_intersection!(
    queue::Vector{<:SegmentEvent},
    ev1::SegmentEvent,
    ev2::SegmentEvent,
    self_intersection::Bool
    ; atol::AbstractFloat=default_atol
    )
    pt = intersect_geometry(ev1.segment, ev2.segment)
    if isnothing(pt)
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
##                 Segment Chaining                        ##
#############################################################

is_empty_segment(ev::SegmentEvent) = (ev.self_annotations.fill_above == false) && (ev.self_annotations.fill_below == false)
struct SegmentChainCandidate{T}
    chain_idx::Int
    match_chain_start::Bool
    match_segment_start::Bool
    match_point::Point2D{T}
    other_point::Point2D{T}
end

function chain_segments(
    segments::AbstractVector{SegmentEvent{T}}
    ; atol::AbstractFloat=default_atol, check_closes::Bool=true
    ) where T
    # Note: if any of the regions intersect at a vertex, than this is not guaranteed to give consistent results
    # They might be joined into one region or presented as separate regions.
    # This algorithm can fail if the polygon is improper (it has lines jutting out)
    chains = Vector{Point2D{T}}[] # this is the same type as regions
    regions = Vector{Point2D{T}}[]
    processed = Set{Segment2D{T}}()
    for event in segments
        if event.segment in processed
            continue
        end
        push!(processed, event.segment)
        candidates = SegmentChainCandidate{T}[]
        @debug("[chain_segment]: event=$event")
        @debug("[chain_segment]: candidates=$candidates")
        for (chain_idx, chain) in enumerate(chains)
            insert_matching_candidate!(candidates, chain, chain_idx, event.segment; atol=atol)
        end
        if length(candidates) == 0 # start a new open chain
            chain = [event.segment[1], event.segment[2]]
            @debug("[chain_segment]: new chain")
            push!(chains, chain)
        elseif length(candidates) == 1 # check if it closes else append to chain
            candidate = candidates[1]
            chain = chains[candidate.chain_idx]
            if check_closes && closes_chain(chain, candidate; atol=atol)
                popat!(chains, candidate.chain_idx)
                push!(regions, chain)
                @debug("[chain_segment]: closed chain")
            else
                append_candidate!(chain, candidate; atol=atol)
                @debug("[chain_segment]: appended chain")
            end
        elseif length(candidates) == 2 # join two chains together
            cand1 = candidates[1]
            cand2 = candidates[2]
            @assert cand1.match_segment_start != cand2.match_segment_start "Same point of segment $(cand1.segment) linked to two open chains"
            chain1 = chains[cand1.chain_idx]
            chain2 = chains[cand2.chain_idx]
            append_candidate!(chain1, cand1; atol=atol)
            new_chain = join_chains!(chain1, chain2, cand1.match_chain_start, cand2.match_chain_start)
            chains[cand1.chain_idx] = new_chain
            @debug("[chain_segment]: combined chains")
            deleteat!(chains, cand2.chain_idx)
        else # confused
            throw("Matched segment $(candidate.segment) to more than 2 chains.")
        end
    end
    # TODO: it might be possible to close some open chains
    # - it is improper: the beginning and end is a segment(s) jutting out, so it can be closed with a segment in processing
    if check_closes
        if !isempty(chains)
            fuzzy_close!(chains, regions; atol=atol)
        end
        @assert isempty(chains) "There are still open chains at the end of processing all segments."
        return regions
    else
        return chains
    end
end

function insert_matching_candidate!(
    candidates::Vector{<:SegmentChainCandidate},
    chain::Vector{<:Point2D},
    chain_idx::Int,
    segment::Segment2D
    ; atol::AbstractFloat=default_atol
    )
    is_match = false
    if is_same_point(chain[1], segment[1]; atol=atol)
        is_match = true
        match_chain_start = true
        match_idx = 1
        other_idx = 2
    elseif is_same_point(chain[1], segment[2]; atol=atol)
        is_match = true
        match_chain_start = true
        match_idx = 2
        other_idx = 1
    elseif is_same_point(chain[end], segment[1]; atol=atol)
        is_match = true
        match_chain_start = false
        match_idx = 1
        other_idx = 2
    elseif is_same_point(chain[end], segment[2]; atol=atol)
        is_match = true
        match_chain_start = false
        match_idx = 2
        other_idx = 1
    end
    if is_match
        candidate = SegmentChainCandidate(
            chain_idx, match_chain_start, match_idx == 1, segment[match_idx], segment[other_idx]
        )
        push!(candidates, candidate)
    end
end

function append_candidate!(chain::Vector{<:Tuple}, candidate::SegmentChainCandidate
    ; atol::AbstractFloat=default_atol)
    if candidate.match_chain_start
        if length(chain) > 1 && 
            get_orientation(candidate.other_point, chain[1], chain[2]; atol=atol) == COLINEAR
            popfirst!(chain)
        end
        insert!(chain, 1, candidate.other_point)
    else
        if length(chain) > 1 && 
            get_orientation(chain[end-1], chain[end], candidate.other_point; atol=atol) == COLINEAR
            pop!(chain)
        end
        push!(chain, candidate.other_point)
    end
end

function closes_chain(chain::Vector{<:Point2D}, candidate::SegmentChainCandidate; atol::AbstractFloat=default_atol)
    if candidate.match_chain_start
        return is_same_point(chain[end], candidate.other_point; atol=atol)
    else
        return is_same_point(chain[1], candidate.other_point; atol=atol)
    end
end

function fuzzy_close!(chains::Vector{<:Vector{<:Point2D}}, regions::Vector{<:Vector{<:Point2D}}; atol)
    for idx in reverse(eachindex(chains))
        if is_fuzzy_closed(chains[idx], length(regions) + 1; atol=atol)
            push!(regions, popat!(chains, idx))
        end
    end
    regions
end

function is_fuzzy_closed(chain::Vector{<:Point2D}, idx::Int; atol::AbstractFloat, rtol::AbstractFloat=1.0)
    if is_same_point(chain[1], chain[end]; atol=atol)
        return true
    end
    gap = norm(chain[1], chain[end])
    gaps = norm.(chain[1:(end-1)], chain[2:end])
    mean_gap = sum(gaps) / length(gaps)
    if gap / mean_gap <= rtol
        @warn("Region $idx was not closed, but it has a relatively small gap and will be considered closed.")
        return true
    end
    false
end

function join_chains!(chain1::Vector{<:Point2D}, chain2::Vector{<:Point2D}, match_chain1_start, match_chain2_start)
    # Note: with clever use of reverse! can change this to always modify chain1 in place for the same output
    if match_chain1_start && match_chain2_start
        # <--- --->
        return push!(reverse!(chain2), chain1...)
    elseif match_chain1_start && !match_chain2_start
        # <--- <----
        return push!(chain2, chain1...)
    elseif !match_chain1_start && match_chain2_start
        # ---> --->
        return push!(chain1, chain2...) 
    else # !match_chain1_start && !match_chain2_start
        # ----> <-----
        return push!(chain1, reverse!(chain2)...) 
    end
end