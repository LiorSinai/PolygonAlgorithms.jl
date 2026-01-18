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
