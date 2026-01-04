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
    # Martinez-Rueda properties
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

function copy_segment(event::SegmentEvent, is_primary::Bool) 
    # make a copy independent of the original event.
    SegmentEvent(
        deepcopy(event.segment), event.is_start, is_primary, deepcopy(event.self_annotations), deepcopy(event.other_annotations)
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

function convert_to_event_queue(
    polygon::Path2D{T}; options...
    ) where T
    # The event list reads all segments from left to right, end to start, bottom to top
    convert_to_event_queue!(SegmentEvent{T}[], polygon; options...)
end

function convert_to_event_queue!(
    queue::Vector{<:SegmentEvent}, polygon::Path2D;
    primary::Bool=true, atol::AbstractFloat=default_atol
    )
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
    start_event = SegmentEvent(segment, true, primary, shared_self_annotations, shared_other_annotations)
    end_event = SegmentEvent(segment, false, primary, shared_self_annotations, shared_other_annotations)
    start_event.other = end_event
    end_event.other = start_event   
    insert_in_order!(queue, start_event; lt=compare_events)
    insert_in_order!(queue, end_event; lt=compare_events)
end

"""
    _compare_points(pt1, pt2; atol=default_atol)
Return:
- -1 if pt1 is smaller
- 0 if the same
- 1 if pt2 is smaller
"""
function _compare_points(pt1::Point2D{T}, pt2::Point2D{T}; atol::AbstractFloat=default_atol) where T # pointsCompare
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

Smaller events are to the left or bottom. Otherwise, end events come before the start.
Returns true if smaller.
"""
function compare_events(event::SegmentEvent, here::SegmentEvent; atol::AbstractFloat=default_atol) # eventCompare
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
    if abs(here.segment[1][1] - here.segment[2][1]) < atol # vertical
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
    is_above(event, other, [atol])


!!!!! Critical function. May be source of errors that only emerge later.

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
    ev::SegmentEvent, other::SegmentEvent
    ; atol::AbstractFloat=default_atol
    ) # statusCompare
    seg1 = ev.segment
    seg2 = other.segment
    if (seg1[1][1] < seg2[1][1])
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

function find_transition(
    list::Vector{<:SegmentEvent}, event::SegmentEvent
    ; atol::AbstractFloat=default_atol
    )
    searchsortedfirst(list, event; lt=(x, y) -> is_above(x, y; atol=atol))
end

function is_vertex_intersection(segment1::Segment2D, segment2::Segment2D; atol::AbstractFloat=default_atol)
    on_segment(segment1[1], segment2; atol=atol) || 
        on_segment(segment1[2], segment2; atol=atol) ||
        on_segment(segment2[1], segment1; atol=atol) ||
        on_segment(segment2[2], segment1; atol=atol)
end

#############################################################
##                  Line Sweep                             ##
#############################################################

"""
    any_intersect(segment::Segment2D, ...; atol=default_atol, include_vertices=true)
    any_intersect(queue::Vector{SegmentEvent}}; atol=default_atol, include_vertices=true)

A line sweep algorithm for determining if any segment intersects with any other segment.

Shamos-Hoey algorithm. It runs in `O(n*log(n))` time where `n` is the length of the segments.

Reference:
- http://euro.ecom.cmu.edu/people/faculty/mshamos/1976GeometricIntersection.pdf
"""
function any_intersect(
    queue::Vector{SegmentEvent{T}}
    ; atol::AbstractFloat=default_atol, include_vertices::Bool=true
    ) where T
    sweep_status = SegmentEvent{T}[] # current events in a vertical line, top to bottom.
    for head in queue
        queue_length = length(queue)
        status_length = length(sweep_status)
        @debug("[do_intersect] ($(queue_length), $(status_length)): $(head)")
        if head.is_start # then check for intersections and add to sweep status
            idx = find_transition(sweep_status, head; atol=atol)
            above = idx == 1 ? nothing : sweep_status[idx - 1]
            below = (idx > length(sweep_status)) ? nothing : sweep_status[idx]
            @debug("[event_loop!] transition idx=$idx")
            @debug("[event_loop!] above=$above")
            @debug("[event_loop!] below=$below")
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

function any_intersect(segments::Vararg{Segment2D{T}}; atol::AbstractFloat=default_atol, options...) where T
    queue = SegmentEvent{T}[]
    for segment in segments
        forward = _compare_points(segment[1], segment[2]; atol=atol)
        segment_ = forward < 0 ? (segment[1], segment[2]) : (segment[2], segment[1])
        add_segment_event!(queue, segment_, true)
    end
    any_intersect(queue; atol=atol, options...)
end
