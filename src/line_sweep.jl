#############################################################
##                  Initialise Events                      ##
#############################################################

"""
    convert_to_events(polygon; primary=true, atol=default_atol)

Convert each edge of `polygon` into a start and end `SegmentEvent` and insert them into a list.
See `convert_to_events!` to add to an existing list.

Segments are oriented so that they always go left to right, or bottom to top if vertical.
Zero length edges are skipped.

`primary` marks whether the segments belong to the subject (`true`) or the clip (`false`) polygon.
"""
function convert_to_events(
    polygon::Path2D{T}; options...
    ) where T
    # The event list reads all segments from left to right, end to start, bottom to top
    convert_to_events!(SegmentEvent{T}[], polygon; options...)
end

"""
    convert_to_events(queue, polygon; primary=true, atol=default_atol)

In-place version of `convert_to_events` which inserts the events into an existing `queue`.
"""
function convert_to_events!(
    events::Vector{<:SegmentEvent}, polygon::Path2D
    ; primary::Bool=true,
    atol::AbstractFloat=default_atol,
    )
    pt2 = polygon[end]
    for i in eachindex(polygon)
        pt1 = pt2
        pt2 = polygon[i]
        forward = _compare_points(pt1, pt2; atol=atol)
        if forward == Int8(0)
            continue # zero length segment
        end
        start = forward < Int8(0) ? pt1 : pt2
        tail = forward < Int8(0) ? pt2 : pt1
        segment = (start, tail)
        push!(events,
            create_segment_pair(
                segment, primary, SegmentAnnotations(forward=forward), SegmentAnnotations()
            )...
        )
    end
    events
end

function sort_event_queue!(events::Vector{<:SegmentEvent}; atol=default_atol)
    sort!(events; lt=(a, b)->compare_events(a,b; atol=atol))
end

"""
    create_segment_pair(
        segment, primary,
        [shared_self_annotations, shared_other_annotations]
        ; atol=default_atol
    )

Create a linked pair of start and end `SegmentEvent`s for `segment`.

The segment should already be ordered left to right (see `_compare_points`).
The two events share the same annotation objects, so updating the annotations on one updates the other.
"""
function create_segment_pair(
    segment::Segment2D,
    primary::Bool,
    shared_self_annotations::SegmentAnnotations=SegmentAnnotations(),
    shared_other_annotations::SegmentAnnotations=SegmentAnnotations(),
    )
    start_event = SegmentEvent(segment, true, primary, shared_self_annotations, shared_other_annotations)
    end_event = SegmentEvent(segment, false, primary, shared_self_annotations, shared_other_annotations)
    start_event.other = end_event
    end_event.other = start_event
    start_event, end_event
end

function create_segment_pair(ev::SegmentEvent)
    create_segment_pair(
        ev.segment, ev.primary,
        ev.self_annotations, ev.other_annotations
    )
end

"""
    _compare_points(pt1, pt2; atol=default_atol)

Compare points left to right, then bottom to top.
Points are compared by their `x` values unless they are vertically aligned (see `is_vertical`),
in which case they are compared by their `y` values.

Return:
- -1 if pt1 is smaller
- 0 if the same (vertically aligned and the y values are within `atol`)
- 1 if pt2 is smaller
"""
function _compare_points(pt1::Point2D{T}, pt2::Point2D{T}; atol::AbstractFloat=default_atol) where T # pointsCompare
    if is_vertical(pt1, pt2)
        if abs(pt1[2] - pt2[2]) <= atol # same point
            return Int8(0)
        end
        return pt1[2] < pt2[2] ? Int8(-1) : Int8(1); # compare Y values
    end
    return pt1[1] < pt2[1] ? Int8(-1) : Int8(1); # compare X values
end

"""
    compare_events(event, here; [atol])

Smaller events are to the left or bottom. Otherwise, end events come before the start.
Returns true if smaller.
"""
function compare_events(
    event::SegmentEvent, here::SegmentEvent
    ; atol::AbstractFloat=default_atol
    ) # eventCompare
    # Assumes events are left to right
    comp = _compare_points(event.point, here.point; atol=atol)
    if comp != Int8(0)
        return comp < Int8(0)
    end
    # Selected points are the same -> events on top of each other.
    comp = _compare_points(event.other_point, here.other_point; atol=atol)
    if comp === Int8(0)
        return false # equal segments
    end
    # Two events on top of each other.
    if event.is_start != here.is_start
        # favor the one that isn't the start
        return event.is_start ? false : true
    end
    # share a common start point ⋅< or a common end point >⋅
    # Manually calculate if the other point is above
    if is_vertical(here.segment)
        # projecting the point won't work.
        # instead, assume smaller segment leans towards the right
        return event.other_point[1] > here.segment[1][1]
    end
    is_above_or_on(event.other_point, here.segment; atol=atol) ? false : true
end

#############################################################
##                  Above                                  ##
#############################################################

"""
    is_above(event, other; atol=default_atol)

!!!!! Critical function. It may be source of errors that only emerge later.

Return `true` if `event` is strictly above `other`.

A segment is considered above another if:
    1. It is to the left of the other segment.
    2. And the start point of the other segment is orientated clockwise from it.
Or symmetrically:
    1. It is in line or to the right of the other segment.
    2. And its start point is orientated counter-clockwise from the other segment.

Assumes segments always go left to right.
"""
function is_above(
    ev::SegmentEvent{T}, other::SegmentEvent{T}
    ; atol::AbstractFloat=default_atol
    ) where T # statusCompare
    seg1 = ev.segment
    seg2 = other.segment
    if (seg2[1][1] - seg1[1][1]) > eps(T) # seg1 starts left of seg2 by more than eps
        orient = get_orientation(seg1[1], seg1[2], seg2[1]; atol=atol)
        if orient == COLINEAR
            orient = get_orientation(seg1[1], seg1[2], seg2[2]; atol=atol)
        end
        return orient == CLOCKWISE
    else
        orient = get_orientation(seg2[1], seg2[2], seg1[1]; atol=atol)
        if orient == COLINEAR
            orient = get_orientation(seg2[1], seg2[2], seg1[2]; atol=atol)
        end
        return orient == COUNTER_CLOCKWISE
    end
end

"""
    find_transition(sweep_status, event; atol=default_atol)

Binary search for the index in the sweep status (ordered top to bottom) where `event` should be inserted,
i.e. the index of the first segment that `event` is above according to `is_above`.
"""
function find_transition(
    list::Vector{<:SegmentEvent}, event::SegmentEvent
    ; atol::AbstractFloat=default_atol
    )
    searchsortedfirst(list, event; lt=(x, y) -> is_above(x, y; atol=atol))
end

"""
    is_vertex_intersection(segment1, segment2; atol=default_atol)

Return `true` if any end point of either segment lies on the other segment.
"""
function is_vertex_intersection(
    segment1::Segment2D, segment2::Segment2D
    ; atol::AbstractFloat=default_atol
    )
    fn(a, b) = on_segment(a, b; atol=atol)
    fn(segment1[1], segment2) || 
        fn(segment1[2], segment2) ||
        fn(segment2[1], segment1) ||
        fn(segment2[2], segment1)
end

#############################################################
##                  Line Sweep                             ##
#############################################################

"""
    any_intersect(segments::Segment2D...; atol=default_atol, include_vertices=true)
    any_intersect(queue::Vector{<:SegmentEvent}; atol=default_atol, include_vertices=true)

A line sweep algorithm for determining if any segment intersects with any other segment.

If `include_vertices=false`, intersections where an end point of one segment lies on another segment are ignored.
This is useful for checking that a polygon does not self-intersect, because adjacent edges always share a vertex.

Shamos-Hoey algorithm. It runs in `O(n*log(n))` time where `n` is the number of segments.

Reference:
- http://euro.ecom.cmu.edu/people/faculty/mshamos/1976GeometricIntersection.pdf
"""
function any_intersect(
    queue::Vector{SegmentEvent{T}}
    ;
    atol::AbstractFloat=default_atol,
    include_vertices::Bool=true
    ) where T
    sort_event_queue!(queue; atol=atol)
    sweep_status = SegmentEvent{T}[] # current events in a vertical line, top to bottom.
    for head in queue
        queue_length = length(queue)
        status_length = length(sweep_status)
        @debug("[any_intersect] ($(queue_length), $(status_length)): $(head)")
        if head.is_start # then check for intersections and add to sweep status
            idx = find_transition(sweep_status, head; atol=atol)
            above = idx == 1 ? nothing : sweep_status[idx - 1]
            below = (idx > length(sweep_status)) ? nothing : sweep_status[idx]
            @debug("[any_intersect] transition idx=$idx")
            @debug("[any_intersect] above=$above")
            @debug("[any_intersect] below=$below")
            if !isnothing(above) && do_intersect(head.segment, above.segment; atol=atol) &&
                (include_vertices || !is_vertex_intersection(head.segment, above.segment; atol=atol))
                return true
            elseif !isnothing(below) && do_intersect(head.segment, below.segment; atol=atol) &&
                (include_vertices || !is_vertex_intersection(head.segment, below.segment; atol=atol))
                return true
            end
            insert!(sweep_status, idx, head)
        else # event is ending, so remove it from the status
            idx = find_transition(sweep_status, head.other; atol=atol)
            if !(0 < idx <= length(sweep_status) && sweep_status[idx] === head.other)
                @warn "$(head.other) was not in the expected location in the sweep status. " * 
                    "Falling back to linear search. Results might not be reliable."
                idx = findfirst(x -> x === head.other, sweep_status)
                @assert(
                    !isnothing(idx),
                    "$(head.other) is missing from the sweep_status. The start event should always be processed before the end event."
                )
            end
            if (idx != 1) && (idx != length(sweep_status))
                # there will be 2 new adjacent edges, so check the intersection between them
                above = sweep_status[idx - 1]
                below = sweep_status[idx + 1]
                if do_intersect(above.segment, below.segment; atol=atol) && 
                    (include_vertices || !is_vertex_intersection(above.segment, below.segment; atol=atol))
                    return true
                end
            end
            popat!(sweep_status, idx)
        end
    end
    false
end

function any_intersect(
    segments::Vararg{Segment2D{T}}
    ; atol::AbstractFloat=default_atol,
    options...
    ) where T
    events = SegmentEvent{T}[]
    for segment in segments
        forward = _compare_points(segment[1], segment[2]; atol=atol)
        segment_ = forward < Int8(0) ? (segment[1], segment[2]) : (segment[2], segment[1])
        push!(events, create_segment_pair(segment_, true)...)
    end
    any_intersect(events; atol=atol, options...)
end
