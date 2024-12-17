@enum AnnotationFill BLANK=0 ABOVE=1 BELOW=2 EMPTY=3

"""
    martinez_rueda_algorithm(polygon1::Vector{<:Point2D}, polygon2::Vector{<:Point2D})

This uses the Martinez-Rueda-Feito polygon clipping algorithm.
It runs in `O((n+m+k)log(n+m))` time where `n` and `m` are the number of vertices of `polygon1`` and `polygon2` respectively.
Use `intersect_convex` for convex polygons for an `O(n+m)` algorithm.

Limitations
1. It can fail for improper polygons: polygons with lines sticking out.

References 
- paper: https://www.researchgate.net/publication/220163820_A_new_algorithm_for_computing_Boolean_operations_on_polygons
- article: https://sean.fun/a/polygon-clipping-pt2/
- article source code: https://github.com/velipso/polybooljs
"""
function martinez_rueda_algorithm(
    polygon1::Polygon2D{T},
    polygon2::Polygon2D{T},
    selection_criteria::Vector{AnnotationFill}
    ; atol::Float64=1e-6
    ) where T
    event_queue1 = convert_to_event_queue(polygon1; primary=true, atol=atol)
    annotated_segments1 = event_loop!(event_queue1, self_intersection=true)
    event_queue2 = convert_to_event_queue(polygon2; primary=false, atol=atol)
    annotated_segments2 = event_loop!(event_queue2, self_intersection=true)
    queue = SegmentEvent{T}[]
    for ev in vcat(annotated_segments1, annotated_segments2)
        add_annotated_segment!(queue, ev)
    end
    annotated_segments3 = event_loop!(queue, self_intersection=false)
    # for consistent reporting, swap annotations so that self annotations are always the primary
    for ev in annotated_segments3
        if !ev.primary
            temp = ev.self_annotations
            ev.self_annotations = ev.other_annotations
            ev.other_annotations = temp
        end
    end
    selected = apply_selection_criteria(annotated_segments3, selection_criteria)
    empty_segments, regions_segments = separate(is_empty_segment, selected)
    regions = chain_segments(regions_segments; atol=atol, check_closes=true)
    segment_chains = chain_segments(empty_segments; atol=atol, check_closes=false)
    # It is also possible to attach some segment_chains to regions.
    # This will give consistent results with the Weiler-Atherton implementation.
    # For now, skipping this step.
    vcat(regions, segment_chains)
end

function insert_in_order!(vec::Vector{T}, data::T; lt=isless, rev=false) where T
    idx = searchsortedfirst(vec, data; lt=lt, rev=rev)
    insert!(vec, idx, data)
end

function pop_key!(vec::Vector{T}, key::T) where T
    idx = findfirst(x->x===key, vec)
    if isnothing(idx)
        throw(KeyError(key))
    end
    popat!(vec, idx)
end

#############################################################
##                     SegmentEvent                        ##
#############################################################

mutable struct SegmentAnnotations
    fill_above::Union{Nothing, Bool}
    fill_below::Union{Nothing, Bool}
end

SegmentAnnotations() = SegmentAnnotations(nothing, nothing)

==(ann1::SegmentAnnotations, ann2::SegmentAnnotations) = 
    (ann1.fill_above == ann2.fill_above) && (ann1.fill_below == ann2.fill_below)

mutable struct SegmentEvent{T}
    segment::Segment2D{T}
    is_start::Bool
    primary::Bool # primary or secondary polygon
    self_annotations::SegmentAnnotations # this polygon
    other_annotations::SegmentAnnotations # other polygon
    # convenience properties
    other::Union{Nothing,SegmentEvent{T}} # links to opposite event
    point::Point2D{T} # is_start ? segment[1] : segment[2]
    other_point::Point2D{T} # is_start ? segment[2] : segment[1]
end

function SegmentEvent(
    segment::Segment2D,
    is_start::Bool,
    primary::Bool=true,
    self_annotations::SegmentAnnotations=SegmentAnnotations(),
    other_annotations::SegmentAnnotations=SegmentAnnotations(),
    )
    point = is_start ? segment[1] : segment[2]
    other_point = is_start ? segment[2] : segment[1]
    SegmentEvent(segment, is_start, primary, self_annotations, other_annotations, nothing, point, other_point)
end

function copy_segment(event::SegmentEvent) 
    # make a copy independent of the original event.
    SegmentEvent(
        deepcopy(event.segment), event.is_start, event.primary, deepcopy(event.self_annotations), deepcopy(event.other_annotations)
    )
end

"""
    ==(ev1::SegmentEvent, ev2::SegmentEvent)

Test equality by comparing attributes but excluding the `other` link.
Use `===` instead if the `other` link is required.
"""
==(ev1::SegmentEvent, ev2::SegmentEvent) = 
    (ev1.segment == ev2.segment) && (ev1.is_start == ev2.is_start) &&
    (ev1.primary == ev2.primary) &&
    (ev1.self_annotations == ev2.self_annotations) && (ev1.other_annotations == ev2.other_annotations)

function Base.show(io::IO, event::SegmentEvent)
    #print(io, typeof(event), "(")
    print(io, "SegmentEvent(")
    print(io, event.segment)
    print(io, ", ", event.is_start)
    print(io, ", ", event.primary)
    print(io, ", ", event.self_annotations)
    print(io, ", ", event.other_annotations)
    #print(io, ", ", event.other)
    #print(io, ", ", event.point)
    #print(io, ", ", event.other_point)
    print(io, ")")
end

#############################################################
##                  Initialise Events                      ##
#############################################################

function convert_to_event_queue(polygon::Polygon2D{T}; primary::Bool=true, atol::Float64=1e-6) where T
    # The event list reads all segments from left to right, end->start, top to bottom
    queue = SegmentEvent{T}[]
    pt2 = polygon[end]
    for i in eachindex(polygon)
        pt1 = pt2
        pt2 = polygon[i]
        forward = _compare_points(pt1, pt2; atol=atol)
        if forward == 0
            continue # zero length segment
        end
        start = forward < 0 ? pt1 : pt2
        end_ = forward < 0 ? pt2 : pt1
        segment = (start, end_)
        add_segment_event!(queue, segment, primary)
    end
    queue
end

function add_segment_event!(
    queue::Vector{<:SegmentEvent},
    segment::Segment2D,
    primary::Bool,
    shared_self_annotations::SegmentAnnotations=SegmentAnnotations(),
    shared_other_annotations::SegmentAnnotations=SegmentAnnotations(),
    )
    forward = _compare_points(segment[1], segment[2])
    @assert forward < 0 "invalid segment $(segment). Require start to be to the left or directly below the end."
    start_event = SegmentEvent(segment, true, primary, shared_self_annotations, shared_other_annotations)
    end_event = SegmentEvent(segment, false, primary, shared_self_annotations, shared_other_annotations)
    start_event.other = end_event
    end_event.other = start_event   
    insert_in_order!(queue, start_event; lt=compare_events)
    insert_in_order!(queue, end_event; lt=compare_events)
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

"""
    _compare_points(pt1, pt2)
Return:
- -1 if pt1 is smaller
- 0 if the same
- 1 if pt2 is smaller
"""
function _compare_points(pt1::Point2D{T}, pt2::Point2D{T}; atol::Float64=1e-6) where T # pointsCompare
    if abs(pt1[1] - pt2[1]) < atol # on a vertical line
        if abs(pt1[2] - pt2[2]) < atol # same point
            return 0
        end
        return pt1[2] < pt2[2] ? -1 : 1; # compare Y values
    end
    return pt1[1] < pt2[1] ? -1 : 1; # compare X values
end

"""
    compare_events(event, here)

Smaller events are to the left or bottom. Otherwise, end events come before the start  
"""
function compare_events(event::SegmentEvent, here::SegmentEvent) # eventCompare
    # Assumes events are left to right
    comp = _compare_points(event.point, here.point)
    if comp != 0
        return comp < 0
    end
    # Selected points are the same -> events on top of each other.
    comp = _compare_points(event.other_point, here.other_point)
    if comp === 0
        return false # equal segments
    end
    # Two events on top of each other.
    if event.is_start != here.is_start
        # favor the one that isn't the start
        return event.is_start ? false : true
    end
    # share a common start point ⋅< or a common end point >⋅
    # Manually calculate if the other point is above
    is_above_or_on(event.other_point, here.segment) ? false : true
end

#############################################################
##                    Event Loop                           ##
#############################################################

function event_loop!(queue::Vector{SegmentEvent{T}}; self_intersection::Bool) where T # eventLoop
    annotated_segments = SegmentEvent{T}[]
    sweep_status = SegmentEvent{T}[] # current events in a vertical line, top to bottom.
    while !(isempty(queue))
        head = queue[1]
        #queue_length = length(queue)
        #status_length = length(sweep_status)
        #println("\nevent_loop ($(queue_length), $(status_length)): $(head)")
        #println("   sweep_status=$sweep_status")
        if head.is_start # then check for intersections and add to sweep status
            idx, above, below = find_transition(sweep_status, head)
            #println("   transition idx=$idx")
            check_and_divide_intersection!(queue, head, above, self_intersection)
            if queue[1] != head
                continue # either head was removed or something was inserted ahead of it
            end
            check_and_divide_intersection!(queue, head, below, self_intersection)
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
            idx = findfirst(x->x===head.other, sweep_status)
            @assert !isnothing(idx) "$(head.other) is missing from the sweep_status. The start event should always be processed before the end event."
            if (idx != 1) && (idx != length(sweep_status))
                # there will be 2 new adjacent edges, so check the intersection between them
                check_and_divide_intersection!(queue, sweep_status[idx - 1], sweep_status[idx + 1], self_intersection)
            end
            push!(annotated_segments, copy_segment(head.other))
            popat!(sweep_status, idx)
        end
        popfirst!(queue)
    end
    annotated_segments
end

function find_transition(list::Vector{<:SegmentEvent}, event::SegmentEvent)
    idx = searchsortedfirst(list, event; lt=is_above)
    above = idx == 1 ? nothing : list[idx - 1]
    below = (idx > length(list)) ? nothing : list[idx]
    idx, above, below
end

function is_above(ev::SegmentEvent, other::SegmentEvent; atol::AbstractFloat=1e-6) # statusCompare
    # Critical function. May be source of errors that only emerge later.
    # Assumes segments always go left to right.
    # Project right most segment's start point on to the line through the segment and compare y values.
    seg1 = ev.segment
    seg2 = other.segment
    ori_start = get_orientation(seg1[1], seg2[1], seg1[2])
    if ori_start == COLINEAR
        return !is_above_or_on(seg2[2], seg1)
    end
    !is_above_or_on(seg2[1], seg1; atol=atol)
end

function check_and_divide_intersection!(queue::Vector{<:SegmentEvent}, ev1::SegmentEvent, ev2::Nothing, self_intersection::Bool; atol=1e-6)
    queue
end

function check_and_divide_intersection!(
    queue::Vector{<:SegmentEvent},
    ev1::SegmentEvent,
    ev2::SegmentEvent,
    self_intersection::Bool
    ; atol::Float64=1e-6
    )
    pt = intersect_geometry(ev1.segment, ev2.segment)
    if isnothing(pt)
        # Lines need to be exactly on top of each other 
        ori2_start = get_orientation(ev1.segment[1], ev1.segment[2], ev2.segment[1])
        ori2_end = get_orientation(ev1.segment[1], ev1.segment[2], ev2.segment[2])
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
    ; atol::Float64=1e-6
    ) # checkIntersection
    #println("   divide_intersection: $(ev1.segment) -- $(ev2.segment) at $(pt)")
    at_start1, at_end1, along1 = classify_intersection(ev1.segment, pt; atol=atol)
    at_start2, at_end2, along2 = classify_intersection(ev2.segment, pt; atol=atol)
    #println("   ", at_start1, " ", at_end1, " ", along1)
    #println("   ", at_start2, " ", at_end2, " ", along2)
    if along1 && along2
        divide_event!(queue, ev1, pt)
        divide_event!(queue, ev2, pt)
    elseif along1
        if at_start2
            divide_event!(queue, ev1, ev2.segment[1])
        elseif at_end2
            divide_event!(queue, ev1, ev2.segment[2])
        end
    elseif along2
        if at_start1
            divide_event!(queue, ev2, ev1.segment[1])
        elseif at_end1
            divide_event!(queue, ev2, ev1.segment[2])
        end
    end
    queue
end

function divide_coincident_intersection!(
    queue::Vector{<:SegmentEvent}, ev1::SegmentEvent, ev2::SegmentEvent, self_intersection::Bool
    ; atol::Float64=1e-6
    )
    # This assumes:
    # - ev1 is on top of or to the right of ev2, because events are processed left to right.
    # - both points of ev2 are colinear with ev1 .
    start1_on_end2 = is_same_point(ev1.segment[1], ev2.segment[2]; atol=atol)
    end1_on_start2 = is_same_point(ev1.segment[2], ev2.segment[1]; atol=atol)
    if start1_on_end2 || end1_on_start2
        return queue # segments touch at endpoints, so no further divisions
    end
    starts_equal = is_same_point(ev1.segment[1], ev2.segment[1]; atol=atol)
    ends_equal = is_same_point(ev1.segment[2], ev2.segment[2]; atol=atol)
    #println("   starts_equal=$starts_equal")
    #println("   ends_equal=$ends_equal")
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
            divide_event!(queue, ev2, ev1.segment[2])
        elseif end2_between
            # (a1)----x-----(a2)
            # (b1)---(b2)
            divide_event!(queue, ev1, ev2.segment[2])
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
                divide_event!(queue, ev2, ev1.segment[2])
            elseif end2_between
                #         (a1)----x-----(a2)
                #  (b1)----------(b2)
                divide_event!(queue, ev1, ev2.segment[2]);
            else # are these segments colinear?
                return queue
            end
        end
        #         (a1)---(a2)
        #  (b1)----x-----(b2)
        # equal segment a1-b2 isn't in the status stack yet, so don't return it
        divide_event!(queue, ev2, ev1.segment[1]);
    end
    queue
end

"""
    divide_event!(queue, ev, pt)

Divide an event `ev` and `ev.other` in `queue` into 4:
```
--x-->  to  --> x-->
```
"""
function divide_event!(queue::Vector{<:SegmentEvent}, ev::SegmentEvent, pt::Point2D) # eventDivide
    # assumes pt lies on ev.segment
    new_segment = (pt, ev.segment[2])
    # println("   divide_event: $(new_segment)")
    e1, e2 = update_end!(ev, pt)
    # fix position of end in queue
    pop_key!(queue, e2)
    insert_in_order!(queue, e2; lt=compare_events)
    # add new segment at the end. Reset other_annotations
    add_segment_event!(queue, new_segment, ev.primary, ev.self_annotations, SegmentAnnotations())
end

"""
    update_end!(queue, ev, pt)

Slides an end backwards.
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

function merge_same_segments!(queue::Vector{<:SegmentEvent}, discard::SegmentEvent, survive::SegmentEvent, self_intersection::Bool)
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
    queue
end

function calculate_self_annotations!(ev::SegmentEvent, below::Union{Nothing, SegmentEvent})
    # if a new segment, than toggle, else use existing knowledge
    #println("   event: $(ev)")
    #println("   below: $(below)")
    toggle = isnothing(ev.self_annotations.fill_below) ? true : ev.self_annotations.fill_above != ev.self_annotations.fill_below
    if isnothing(below)
        ev.self_annotations.fill_below = false # TODO primaryPolyInverted
    else
        @assert !isnothing(below.self_annotations.fill_above) "missing annotations below: $(below)" # preempt !nothing error
        ev.self_annotations.fill_below = below.self_annotations.fill_above # below should already be filled
    end
    if toggle
        ev.self_annotations.fill_above = !ev.self_annotations.fill_below
    else
        ev.self_annotations.fill_above = ev.self_annotations.fill_below
    end
    # println("   self_annotations: ", ev.self_annotations)
    ev.self_annotations
end

function calculate_other_annotations!(ev::SegmentEvent, below::Nothing)
    #println("calculate_other_annotations")
    #println("   ev=$ev")
    #println("   below=$below")
    if isnothing(ev.other_annotations.fill_above)
        # if nothing is below ev, only inside if the other polygon is inverted
        inside = false # TODO is_primary ? secondaryPolyInverted : primaryPolyInverted
        ev.other_annotations.fill_above = inside
        ev.other_annotations.fill_below = inside
    end
    #println("   ev.other_annotations=$(ev.other_annotations)")
    ev.other_annotations
end
    
function calculate_other_annotations!(ev::SegmentEvent, below::SegmentEvent)
    #println("calculate_other_annotations")
    #println("   ev=$ev")
    #println("   below=$below")
    if isnothing(ev.other_annotations.fill_above)
        # something is below ev, so copy the below segment's other polygon's above
        inside = (ev.primary == below.primary) ? 
            below.other_annotations.fill_above : below.self_annotations.fill_above
        ev.other_annotations.fill_above = inside
        ev.other_annotations.fill_below = inside
    end
    #println("   ev.other_annotations=$(ev.other_annotations)")
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

INTERSECTION_CRITERIA = [
    BLANK, BLANK, BLANK, BLANK,
    BLANK, BELOW, EMPTY, BELOW,
    BLANK, EMPTY, ABOVE, ABOVE,
    BLANK, BELOW, ABOVE, BLANK,
] # both below, both above, but not all 4

INTERSECTION_SEGMENT_CRITERIA = [
    BLANK, BLANK, ABOVE, BLANK,
    BLANK, BLANK, EMPTY, BLANK,
    BLANK, EMPTY, BLANK, BLANK,
    BLANK, BLANK, BLANK, BLANK,
] # above and below different segments

UNION_CRITERIA = [
    BLANK, BELOW, ABOVE, BLANK,
    BELOW, BELOW, BLANK, BLANK,
    ABOVE, BLANK, ABOVE, BLANK,
    BLANK, BLANK, BLANK, BLANK,
] # filled at once or twice on only 1 side

DIFFERENCE_CRITERIA = [
    BLANK, BLANK, BLANK, BLANK,
    BELOW, BLANK, BELOW, BLANK,
    ABOVE, ABOVE, BLANK, BLANK,
    BLANK, ABOVE, BELOW, BLANK,
] # primary - secondary. 

XOR_CRITERIA = [
    BLANK, BELOW, ABOVE, BLANK,
    BELOW, BLANK, BLANK, ABOVE,
    ABOVE, BLANK, BLANK, BELOW,
    BLANK, ABOVE, BELOW, BLANK,
]

function apply_selection_criteria(annotated_segments::Vector{<:SegmentEvent{T}}, criteria::Vector{AnnotationFill}) where T
    result = SegmentEvent{T}[]
    for ev in annotated_segments
        index = (ev.self_annotations.fill_above ? 9 : 1) +
                (ev.self_annotations.fill_below ? 4 : 0) + 
                (ev.other_annotations.fill_above ? 2 : 0) + 
                (ev.other_annotations.fill_below ? 1 : 0)
        if criteria[index] != BLANK
            new_segment = copy_segment(ev)
            new_segment.self_annotations.fill_above = criteria[index] == ABOVE
            new_segment.self_annotations.fill_below = criteria[index] == BELOW
            new_segment.other_annotations.fill_above = nothing
            new_segment.other_annotations.fill_below = nothing
            push!(result, copy_segment(new_segment))
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

function chain_segments(segments::AbstractVector{SegmentEvent{T}}; atol::Float64=1e-6, check_closes::Bool=true) where T
    # Note: if any of the regions intersect at a vertex, than this is not guaranteed to give consistent results
    # They might be joined into one region or presented as separate regions.
    # This algorithm can fail if the polygon is improper (it has lines jutting out)
    chains = Vector{Point2D{T}}[] # this is the same type as regions
    regions = Vector{Point2D{T}}[]
    processed = Set{Segment2D{T}}()
    for event in segments
        if event.segment in processed
            #println("processed")
            continue
        end
        push!(processed, event.segment)
        candidates = SegmentChainCandidate{T}[]
        #println("\nevent=$event")
        #println("candidates=$candidates")
        for (chain_idx, chain) in enumerate(chains)
            insert_matching_candidate!(candidates, chain, chain_idx, event.segment; atol=atol)
        end
        if length(candidates) == 0 # start a new open chain
            chain = [event.segment[1], event.segment[2]]
            #println("   new chain: $chain")
            push!(chains, chain)
        elseif length(candidates) == 1 # check if it closes else append to chain
            candidate = candidates[1]
            chain = chains[candidate.chain_idx]
            if check_closes && closes_chain(chain, candidate; atol=atol)
                popat!(chains, candidate.chain_idx)
                push!(regions, chain)
                #println("   closed chain: $chain")
            else
                append_candidate!(chain, candidate)
                #println("   appended chain: $chain")
            end
        elseif length(candidates) == 2 # join two chains together
            cand1 = candidates[1]
            cand2 = candidates[2]
            @assert cand1.match_segment_start != cand2.match_segment_start "Same point of segment $(cand1.segment) linked to two open chains"
            chain1 = chains[cand1.chain_idx]
            chain2 = chains[cand2.chain_idx]
            append_candidate!(chain1, cand1)
            new_chain = join_chains!(chain1, chain2, cand1.match_chain_start, cand2.match_chain_start)
            chains[cand1.chain_idx] = new_chain
            #println("   combined chains: $new_chain")
            deleteat!(chains, cand2.chain_idx)
        else # confused
            throw("Matched segment $(candidate.segment) to more than 2 chains.")
        end
    end
    # TODO: it might be possible to close some open chains
    # - it is improper: the beginning and end is a segment(s) jutting out, so it can be closed with a segment in processing
    if check_closes
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
    ; atol::Float64
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

function append_candidate!(chain::Vector{<:Point2D}, candidate::SegmentChainCandidate)
    if candidate.match_chain_start
        insert!(chain, 1, candidate.other_point)
    else
        push!(chain, candidate.other_point)
    end
end

function closes_chain(chain::Vector{<:Point2D}, candidate::SegmentChainCandidate; atol::Float64)
    if candidate.match_chain_start
        return is_same_point(chain[end], candidate.other_point; atol=atol)
    else
        return is_same_point(chain[1], candidate.other_point; atol=atol)
    end
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