"""
    SegmentAnnotations(
        [fill_above, fill_below]
        ; forward=-1, winding_top_to_bottom=nothing, winding_left_to_right=nothing, winding_below=0
    )

Annotations of a segment for one polygon, used by the Martinez-Rueda algorithm.
A value of `nothing` means it has not been calculated yet.

Fields:
- `fill_above`: whether the region directly above the segment is inside the polygon.
- `fill_below`: whether the region directly below the segment is inside the polygon.
- `forward`: the original direction of the segment, as returned by `_compare_points(pt1, pt2)`.
    `-1` if the segment was originally left to right, `1` if it was reversed.
- `winding_top_to_bottom`: change in winding number when crossing the segment from top to bottom.
    If the segment represents multiple coincident edges, this is the sum of their changes.
- `winding_left_to_right`: change in winding number when crossing the segment from left to right.
    Used for vertical segments.
- `winding_below`: the winding number of the region directly below the segment.
"""
mutable struct SegmentAnnotations
    fill_above::Union{Nothing, Bool}
    fill_below::Union{Nothing, Bool}
    # fill rule properties
    forward::Int8 # original direction of segment. -1 is current, 0 is zero-length, 1 is reversed
    winding_top_to_bottom::Union{Nothing,Int8}
    winding_left_to_right::Union{Nothing,Int8}
    winding_below::Int
end

function SegmentAnnotations(
    fill_above::Union{Nothing, Bool}=nothing,
    fill_below::Union{Nothing, Bool}=nothing
    ; forward::Int8=Int8(-1),
    winding_top_to_bottom::Union{Nothing,Int8}=nothing,
    winding_left_to_right::Union{Nothing,Int8}=nothing,
    winding_below::Int=0,
    )
    SegmentAnnotations(
        fill_above, fill_below,
        forward,
        winding_top_to_bottom, winding_left_to_right, winding_below
    )
end

==(ann1::SegmentAnnotations, ann2::SegmentAnnotations) = 
    (ann1.fill_above == ann2.fill_above) && (ann1.fill_below == ann2.fill_below)

"""
    SegmentEvent(
    segment, is_start, primary=true,
    [self_annotations, other_annotations]
    )

An event in a line sweep algorithm marking a change in state, either the
start of a segment or the end of it.

Each segment has two events, a start and an end, which are linked to each other with the `other` field.
Both events share the same `segment`, which goes left to right (see `_compare_points`).
`point` is the point where the event occurs and `other_point` is the opposite end of the segment.

For Martinez-Rueda, `primary` indicates if the segment belongs to the subject polygon (`true`)
or a clip polygon (`false`), and `self_annotations`/`other_annotations` hold the fill annotations
for its own polygon and the other polygon respectively. See `SegmentAnnotations`.

See `any_intersect` and `martinez_rueda_algorithm`.
"""
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

eltype(ev::SegmentEvent{T}) where T = T

function SegmentEvent(
    segment::Segment2D,
    is_start::Bool,
    primary::Bool=true,
    self_annotations::SegmentAnnotations=SegmentAnnotations(),
    other_annotations::SegmentAnnotations=SegmentAnnotations(),
    )
    point = is_start ? segment[1] : segment[2]
    other_point = is_start ? segment[2] : segment[1]
    SegmentEvent(
        segment, is_start, primary,
        self_annotations, other_annotations,
        nothing,
        point, other_point,
    )
end

copy(ann::SegmentAnnotations) = SegmentAnnotations(
    ann.fill_above, ann.fill_below,
    ann.forward,
    ann.winding_top_to_bottom, ann.winding_left_to_right, ann.winding_below,
)

"""
    copy(event::SegmentEvent)

Make a copy of `event` that is independent of the original.
The copy is not linked to an `other` event.
"""
function copy(event::SegmentEvent)
    # make a copy independent of the original event.
    SegmentEvent(
        event.segment, # immutable
        event.is_start,
        event.primary,
        copy(event.self_annotations),
        copy(event.other_annotations),
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

getindex(ev::SegmentEvent, idx::Integer) = ev.segment[idx]


"""
    AnnotatedSegment(segment, annotations)
    AnnotatedSegment(point1, point2, [annotations])

A directed segment with the fill annotations of its polygon. Unlike a `SegmentEvent`, the
segment can go in any direction.
"""
struct AnnotatedSegment{T}
    segment::Segment2D{T}
    self_annotations::SegmentAnnotations
end

function AnnotatedSegment(
    point1::Point2D,
    point2::Point2D,
    ann::SegmentAnnotations=SegmentAnnotations(),
    )
    AnnotatedSegment((point1, point2), ann)
end

getindex(seg::AnnotatedSegment, idx::Integer) = seg.segment[idx]

==(seg1::AnnotatedSegment, seg2::AnnotatedSegment) = 
    (seg1.segment == seg2.segment) &&
    (seg1.self_annotations == seg2.self_annotations)

reverse(segment::AnnotatedSegment) = AnnotatedSegment(reverse(segment.segment), segment.self_annotations)

isless(seg1::AnnotatedSegment, seg2::AnnotatedSegment) = isless(seg1.segment, seg2.segment) # used for sorting. See cyclic_set_equality
