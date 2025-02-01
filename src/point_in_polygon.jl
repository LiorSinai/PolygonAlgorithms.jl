"""
    contains(vertices, point; on_border_is_inside=true; rtol=1e-4, atol=1e-6)

Runs in `O(n)` time where `n=length(vertices)`.

This algorithm is an an extension of the odd-even ray algorithm.
It is based on "A Simple and Correct Even-Odd Algorithm for the Point-in-Polygon Problem for Complex Polygons" 
by Michael Galetzka and Patrick Glauner (2017). 
It skips vertices that are on the ray. To compensate, the ray is projected backwards (to the left) so that an 
intersection can be found for a skipped vertex if needed.
"""
function contains(
    vertices::Polygon2D, point::Point2D{T}
    ; on_border_is_inside::Bool=true, rtol::AbstractFloat=1e-4, atol::AbstractFloat=1e-6
    ) where T
    n = length(vertices)
    num_intersections = 0

    x = x_coords(vertices)
    extreme_left =  (minimum(x) - one(T), point[2])
    extreme_right = (maximum(x) + one(T), point[2])

    # step 1: point intersects a vertex or edge
    for i in 1:n
        next_i = (i % n) + 1
        segment = (vertices[i], vertices[next_i])
        if (point == segment[1]) || on_segment(point, segment)
            return on_border_is_inside
        end
    end

    # step 3: check intersections with vertices
    s = 1
    while s <= n
        # step 2: find a vertex not on the same horizontal ray as the point
        while (s <= n) && (vertices[s][2] == point[2])
            s += 1
        end
        if s > n
            break
        end
        # step 3a: find the next vertex not on the horizontal ray
        next_s = s
        skipped_right = false
        for i in 0:n
            next_s = (next_s) % n + 1
            if abs(vertices[next_s][2] - point[2]) > atol && 
                !is_same_point(vertices[s], vertices[next_s]; atol=atol)
                break
            end
            skipped_right = skipped_right || (vertices[next_s][1] > point[1])
        end
        # step 3b: edge intersect with the ray
        edge = (vertices[s], vertices[next_s])
        intersect = false
        if (next_s - s) == 1 || (s == n && next_s ==1) # 3b.i
            intersect = do_intersect(edge, (point, extreme_right); rtol=rtol)
        elseif skipped_right # 3b.ii
            intersect = do_intersect(edge, (extreme_left, extreme_right); rtol=rtol)
        end
        num_intersections += intersect
        if next_s <= s  # gone in a full loop
            break
        end
        s = next_s
    end
    return (num_intersections % 2) == 1
end
