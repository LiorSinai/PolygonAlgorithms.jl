#========= windings for fill rules =========#

"""
Fill rules determine which regions of a self-intersecting polygon are considered filled.

EVEN_ODD rule: alternative regions are filled.

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

function calc_winding_top_to_bottom(start::Point2D, tail::Point2D ; atol::AbstractFloat=default_atol)
    dx = tail[1] - start[1]
    if abs(dx) < atol
        return Int8(0)
    end
    dx > Int8(0) ? Int8(1) : Int8(-1)
end

function calc_winding_left_to_right(start::Point2D, tail::Point2D; atol::AbstractFloat=default_atol)
    dy = tail[2] - start[2]
    if abs(dy) < atol
        return Int8(0)
    end
    dy > Int8(0) ? Int8(1) : Int8(-1)
end

function get_winding_top_to_bottom!(ev::SegmentEvent; atol::AbstractFloat=default_atol)
    if isnothing(ev.winding_top_to_bottom)
        pt1 = ev.forward < Int(0) ? ev.segment[1] : ev.segment[2]
        pt2 = ev.forward < Int(0) ? ev.segment[2] : ev.segment[1]
        ev.winding_top_to_bottom = calc_winding_top_to_bottom(pt1, pt2; atol=atol)
        if !isnothing(ev.other)
            ev.other.winding_top_to_bottom = ev.winding_top_to_bottom
        end
    end
    ev.winding_top_to_bottom
end

function get_winding_left_to_right!(ev::SegmentEvent; atol::AbstractFloat=default_atol)
    if isnothing(ev.winding_left_to_right)
        pt1 = ev.forward < Int(0) ? ev.segment[1] : ev.segment[2]
        pt2 = ev.forward < Int(0) ? ev.segment[2] : ev.segment[1]
        ev.winding_left_to_right = calc_winding_left_to_right(pt1, pt2; atol=atol)
        if !isnothing(ev.other)
            ev.other.winding_left_to_right = ev.winding_left_to_right
        end
    end
    ev.winding_left_to_right
end
