"""
    contains(vertices, point;
        on_border_is_inside=true,
        atol=default_atol,
    )

Runs in `O(n)` time where `n=length(vertices)`.

This algorithm is an an extension of the odd-even ray algorithm.
It is based on "A Simple and Correct Even-Odd Algorithm for the Point-in-Polygon Problem for Complex Polygons" 
by Michael Galetzka and Patrick Glauner (2017). 
It skips vertices that are on the ray. To compensate, the ray is projected backwards (to the left) so that an 
intersection can be found for a skipped vertex if needed.
"""
function contains(
    vertices::Path2D, point::Point2D{T}
    ; 
    on_border_is_inside::Bool=true,
    atol::AbstractFloat=default_atol,
    ) where T
    # First do a bounding box check: is it near the polygon?
    xmin, ymin, xmax, ymax = bounds(vertices)
    if point[1] < xmin - atol || point[1] > xmax + atol || point[2] < ymin - atol || point[2] > ymax + atol
        return false
    end

    # step 1: point intersects a vertex or edge
    if on_border(vertices, point; atol=atol)
        return on_border_is_inside
    end
    
    # step 3: check intersections with vertices
    s = 1
    n = length(vertices)
    num_intersections = 0
    extreme_left = (xmin - one(T), point[2])
    extreme_right = (xmax + one(T), point[2])
    while s <= n
        # step 2: find a vertex not on the same horizontal ray as the point
        while (s <= n) && (abs(vertices[s][2] - point[2]) <= atol)
            s += 1
        end
        if s > n
            break
        end
        # step 3a: find the next vertex not on the horizontal ray
        next_s = s
        skipped_right = false
        for i in 0:(n-1)
            next_s = (next_s) % n + 1
            if abs(vertices[next_s][2] - point[2]) > atol && 
                !is_same_point(vertices[s], vertices[next_s]; atol=atol)
                break
            end
            skipped_right = skipped_right || (vertices[next_s][1] > point[1])
        end
        if next_s == s
            # gone in a full loop
            break
        end
        # step 3b: edge intersect with the ray
        edge = (vertices[s], vertices[next_s])
        intersect = false
        if skipped_right # 3b.ii
            intersect = do_intersect(edge, (extreme_left, extreme_right); atol=atol)
        elseif (next_s - s) == 1 || (s == n && next_s == 1) # 3b.i
            intersect = do_intersect(edge, (point, extreme_right); atol=atol)
        end
        num_intersections += intersect
        if next_s <= s  # gone in a full loop
            break
        end
        s = next_s
    end
    return (num_intersections % 2) == 1
end

"""
    on_border(vertices, point; atol=default_atol)

Return true if `point` lies on any segment between the `vertices`.

Runs in `O(n)` time where `n=length(vertices)`.
"""
function on_border(
    vertices::Path2D, point::Point2D{T}
    ; atol::AbstractFloat=default_atol
    ) where T
    n = length(vertices)
    for i in 1:n
        next_i = (i % n) + 1
        segment = (vertices[i], vertices[next_i])
        if on_segment(point, segment; atol=atol)
            return true
        end
    end
    false
end
