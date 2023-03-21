"""
    point_in_polygon(point, vertices)

Runs in `O(n)` time when `n=length(vertices)`.

This algorithm is an an extension of the odd-even ray algorithm.
It is based on "A Simple and Correct Even-Odd Algorithm for the Point-in-Polygon Problem for Complex Polygons" 
by Michael Galetzka and Patrick Glauner (2017). 
It skips vertices that are on the ray. To compensate, the ray is projected backwards (to the left) so that an 
intersection can be found for a skipped vertix if needed.
"""
function point_in_polygon(point::Point2D{T}, vertices::Polygon2D; on_border_is_inside=true) where T
    n = length(vertices)
    num_intersections = 0

    x = x_coords(vertices)
    extreme_left =  (minimum(x) - T(1e6), point[2])
    extreme_right = (maximum(x) + T(1e6), point[2])

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
        # step 2: find a vertix not on the same horizontal ray as the point
        while (s <= n) && (vertices[s][2] == point[2])
            s += 1
        end
        if s > n
            break
        end
        # step 3a: find the next vertix not on the horizontal ray
        next_s = s
        skipped_right = false
        for i in 0:n
            next_s = (next_s) % n + 1
            if vertices[next_s][2] != point[2]
                break
            end
            skipped_right = skipped_right || (vertices[next_s][1] > point[1])
        end
        # step 3b: edge intersect with the ray
        edge = (vertices[s], vertices[next_s])
        intersect = 0
        if (next_s - s) == 1 || (s == n && next_s ==1) # 3b.i
            intersect = do_intersect(edge, (point, extreme_right))
        elseif skipped_right # 3b.ii
            intersect = do_intersect(edge, (extreme_left, extreme_right))
        end
        num_intersections += intersect
        if next_s <= s  # gone in a full loop
            break
        end
        s = next_s
    end
    return (num_intersections % 2) == 1
end
