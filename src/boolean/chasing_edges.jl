"""
    chasing_edges_algorithm(polygon1, poly2; atol, rtol):

As fast algorithm for finding the intersection, if it exists, between two
convex polygons.

- Time complexity: `O(n+m)` where `n` and `m` vertices on polygon 1 and 2 respectively.
- Description: this algorithm rotates two pointers, one around each polygon. Each iteration it only
    advances one pointer based on a set of "advance rules". After two cycles around the polygons 
    it is guaranteed to have found all (zero, one or both) intersection points and 
    all the points in between.
- Algorithm is from "A New Linear Algorithm for Intersecting Convex Polygons" (1981) by Joseph O'Rourke et. al.
- Reference: https://www.cs.jhu.edu/~misha/Spring16/ORourke82.pdf
"""
function chasing_edges_algorithm(
    polygon1::Path2D{T}, polygon2::Path2D{T}
    ; atol::AbstractFloat=default_atol, rtol::AbstractFloat=default_rtol
    ) where T
    n = length(polygon1)
    m = length(polygon2)
    points = Point2D{T}[]
    if is_clockwise(polygon1)
        polygon1 = reverse(polygon1)
    end
    if is_clockwise(polygon2)
        polygon2 = reverse(polygon2)
    end
    poly1_in_2 = false
    poly2_in_1 = false
    i = 1
    j = 1
    for k in 1:(2 * (m + n))
        i_prev = i == 1 ? n : i - 1
        j_prev = j == 1 ? m : j - 1
        edge1 = (polygon1[i_prev], polygon1[i])
        edge2 = (polygon2[j_prev], polygon2[j])
        inter = intersect_geometry(edge1, edge2; atol=atol, rtol=rtol)
        is_colinear = isapprox(cross_product(edge1, edge2), 0.0; atol=atol)
        if !isnothing(inter) && !is_colinear
            is_second_iter = k > (m + n)
            if length(points) > 1 && (all(inter .≈ points[1]) && is_second_iter)
                poly1_in_2 = false
                poly2_in_1 = false
                break
            end
            push!(points, inter)
            poly1_in_2 = in_half_plane(edge2, polygon1[i])
            poly2_in_1 = !poly1_in_2
        end
        advance_1 = false 
        if cross_product(edge2, edge1) >= 0
            advance_1 = !(in_half_plane(edge2, polygon1[i], ))
        else
            advance_1 = in_half_plane(edge1, polygon2[j])
        end
        if advance_1
            if poly1_in_2
                push!(points, polygon1[i])
            end
            i = i % n + 1
        else # advance_2
            if poly2_in_1
                push!(points, polygon2[j])
            end
            j = j % m + 1
        end
    end
    if isempty(points)
        if contains(polygon2, polygon1[1]; atol=atol, rtol=rtol)
            return polygon1
        elseif contains(polygon1, polygon2[1]; atol=atol, rtol=rtol)
            return polygon2
        end
    end
    points
end
