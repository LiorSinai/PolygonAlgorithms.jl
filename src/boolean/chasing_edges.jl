function chasing_edges_algorithm(polygon1::Polygon2D{T}, polygon2::Polygon2D{T}) where T
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
        inter = intersect_geometry(edge1, edge2)
        is_colinear = cross_product(edge1, edge2) ≈ 0.0
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
        if contains(polygon2, polygon1[1])
            return polygon1
        elseif contains(polygon1, polygon2[1])
            return polygon2
        end
    end
    points
end
