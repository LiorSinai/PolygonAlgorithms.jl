#========= windings for fill rules =========#

"""
Fill rules determine which regions of a self-intersecting polygon are considered filled.

EVEN_ODD rule: alternating regions are filled.

For winding rules:
1. Construct a ray heading out from a given point P in any direction towards infinity.
1. Find all the intersections of the contour C with this ray. 
1. Score up the winding number as follows:
    - for every clockwise intersection (line heading left to right through the ray) subtract 1.
    - for every counter-clockwise intersection (line heading right to left through the ray) add 1.
1. Use the winding number to determine if filled or not.
    - NON_ZERO: non-zero windings are filled.
    - POSITIVE: only positive windings are filled (top to bottom ray).
    - NEGATIVE: only negative winding are filled (top to bottom ray).
References:
- https://www.angusj.com/clipper2/Docs/Units/Clipper/Types/FillRule.htm
- https://en.wikipedia.org/wiki/Nonzero-rule
"""
@enum FillRule EVEN_ODD NON_ZERO POSITIVE NEGATIVE

"""
    calc_winding_top_to_bottom(start, tail)

Change in winding number when a ray crosses the directed edge `start → tail` from top to bottom.

Returns `1` if the edge goes left to right, `-1` if it goes right to left
and `0` if it is vertical (within `eps(T)`).
"""
function calc_winding_top_to_bottom(start::Point2D{T}, tail::Point2D{T}) where T
    # going from top to bottom, we cross right to left from the ray's perspective
    # when the ray goes left to right
    dx = tail[1] - start[1]
    if abs(dx) <= eps(T)
        return Int8(0)
    end
    dx > Int8(0) ? Int8(1) : Int8(-1)
end

"""
    calc_winding_left_to_right(start, tail)

Change in winding number when a ray crosses the directed edge `start → tail` from left to right.

Returns `1` if the edge goes bottom to top, `-1` if it goes top to bottom
and `0` if it is horizontal (within `eps(T)`).
"""
function calc_winding_left_to_right(start::Point2D{T}, tail::Point2D{T}) where T
    dy = tail[2] - start[2]
    if abs(dy) <= eps(T)
        return Int8(0)
    end
    dy > Int8(0) ? Int8(1) : Int8(-1)
end

"""
    get_winding_top_to_bottom!(ev)

Return `ev.self_annotations.winding_top_to_bottom`, calculating it with `calc_winding_top_to_bottom`
if it is `nothing`.
The original direction of the edge is recovered from the `forward` annotation.
The calculated value is cached on both `ev` and `ev.other`.
"""
function get_winding_top_to_bottom!(ev::SegmentEvent)
    annotations = ev.self_annotations
    if isnothing(annotations.winding_top_to_bottom)
        pt1 = annotations.forward < Int(0) ? ev.segment[1] : ev.segment[2]
        pt2 = annotations.forward < Int(0) ? ev.segment[2] : ev.segment[1]
        annotations.winding_top_to_bottom = calc_winding_top_to_bottom(pt1, pt2)
        if !isnothing(ev.other)
            ev.other.self_annotations.winding_top_to_bottom = annotations.winding_top_to_bottom
        end
    end
    annotations.winding_top_to_bottom
end

"""
    get_winding_left_to_right!(ev)

Return `ev.self_annotations.winding_left_to_right`, calculating it with `calc_winding_left_to_right`
if it is `nothing`.
The original direction of the edge is recovered from the `forward` annotation.
The calculated value is cached on both `ev` and `ev.other`.
"""
function get_winding_left_to_right!(ev::SegmentEvent)
    annotations = ev.self_annotations
    if isnothing(annotations.winding_left_to_right)
        pt1 = annotations.forward < Int(0) ? ev.segment[1] : ev.segment[2]
        pt2 = annotations.forward < Int(0) ? ev.segment[2] : ev.segment[1]
        annotations.winding_left_to_right = calc_winding_left_to_right(pt1, pt2)
        if !isnothing(ev.other)
            ev.other.self_annotations.winding_left_to_right = annotations.winding_left_to_right
        end
    end
    annotations.winding_left_to_right
end

"""
    get_winding_change!(ev)

Change in winding number when crossing `ev` from the region below to the region above.
For a vertical segment the winding does not change along the y-axis, so instead use the change
when crossing it from left to right.

If `ev` represents multiple coincident edges, this is the sum of their changes.
"""
function get_winding_change!(ev::SegmentEvent)
    winding = get_winding_top_to_bottom!(ev)
    winding == 0 ? get_winding_left_to_right!(ev) : winding
end

"""
    is_filled(winding, fill_rule)

Whether a region with winding number `winding` is filled according to `fill_rule`.
The even-odd rule is equivalent to an odd winding number because each edge crossed changes 
the winding number by ±1.
"""
function is_filled(winding::Integer, fill_rule::FillRule)
    if fill_rule == EVEN_ODD
        return isodd(winding)
    elseif fill_rule == NON_ZERO
        return winding != 0
    elseif fill_rule == POSITIVE
        return winding > 0
    elseif fill_rule == NEGATIVE
        return winding < 0
    end
    throw("Unknown fill_rule: $fill_rule")
end
