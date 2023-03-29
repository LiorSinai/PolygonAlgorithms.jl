"""
    intersect_geometry(polygon1::Vector{<:Point2D}, polygon2::Vector{<:Point2D})

Returns multiple regions, edges and single points of intersection. 
Only returns the larger type if one is within another e.g. an edge is also part of a region.

This uses the Weiler-Atherton algorithm.
It runs in `O(nm)` time where `n` and `m` are the number of vertices of polygon1 and polygon2 respectively.
Use `intersect_convex` for convex polygons for an `O(n+m)` algorithm.

Limitations
1. This version does not cater for holes.
2. It partially fails for self-intersecting areas. For example, a shared edge that connects to a region of intersection.
3. It can fail completely for self-intersecting polygons.

For a more general algorithm see the Martinez-Rueda polygon clipping algorithm.
"""
function intersect_geometry(polygon1::Polygon2D{T}, polygon2::Polygon2D{T}) where T
    if polygon1 == polygon2
        return [polygon1]
    end
    n = length(polygon1)
    m = length(polygon2)
    @assert n > 2 &&  m > 2 "require at least 3 points per polygons; lengths are $n and $m"
    if is_counter_clockwise(polygon1)
        polygon1 = reverse(polygon1)
    end
    if is_counter_clockwise(polygon2)
        polygon2 = reverse(polygon2)
    end
    # convert polygons to list so can follow loops betwween them
    list1 = generate_list(polygon1)
    list2 = generate_list(polygon2)
    # insert intersections
    find_and_insert_intersections!(list1, list2)
    # walk loops
    regions, visited = walk_linked_lists(list2)
    # special cases
    if isempty(regions)
        p1 = get_first_non_intersection_point(list1)
        p2 = get_first_non_intersection_point(list2)
        if point_in_polygon(p1, polygon2)
            return [polygon1]
        elseif point_in_polygon(p2, polygon1)
            return [polygon2]
        end
    end
    not_visited = get_unvisited_intercepts(list2, visited)
    push!(regions, not_visited...)
    regions
end

@enum IntersectionType NONE ENTRY EXIT VERTIX

mutable struct PointInfo{T}
    point::Point2D{T}
    intersection::Bool
    link::Union{Nothing,Node{PointInfo{T}}}
    type::IntersectionType
end

function generate_list(polygon::Polygon2D{T}) where T
    data = PointInfo{T}(polygon[1], false, nothing, NONE)
    head = Node(data)
    list = DoublyLinkedList(head)
    for point in polygon[2:end]
        data = PointInfo{T}(point, false, nothing, NONE)
        push!(list, data)
    end
    list
end

function get_first_non_intersection_point(polygon::DoublyLinkedList{<:PointInfo})
    start = polygon.head
    node = start
    while node.data.intersection && node.next != start
        node = node.next
    end
    node.data.point
end

function find_and_insert_intersections!(
        polygon1::DoublyLinkedList{PointInfo{T}}, 
        polygon2::DoublyLinkedList{PointInfo{T}};
        atol::Float64=1e-6
        ) where T
    ## collect original nodes before mutating in place
    vec1 = collect_nodes(polygon1)
    vec2 = collect_nodes(polygon2)
    ## find itersections
    for (node1, next1) in zip(vec1, vcat(vec1[2:end], vec1[1]))
        edge1 = (node1.data.point, next1.data.point) 
        for (node2, next2) in zip(vec2, vcat(vec2[2:end], vec2[1]))
            edge2 = (node2.data.point, next2.data.point)
            p = intersect_geometry(edge1, edge2)
            if !isnothing(p)
                i1 = insert_intersection_in_order!(p, node1, next1; atol=atol)
                i2 = insert_intersection_in_order!(p, node2, next2; atol=atol)
                link_intersections!(i1, i2, edge1, edge2; atol=atol)
            end
        end
    end
end

function insert_intersection_in_order!(
    point::Point2D, tail::Node{<:PointInfo}, head::Node{<:PointInfo}
    ; atol::Float64=1e-6
    )
    node = tail
    while node != head.next
        d1 = norm(point, node.data.point)
        d2 = norm(node.next.data.point, node.data.point)
        if d1 < atol # either on tail or intersection here already
            node.data.intersection = true
            return node
        elseif abs(d1 - d2) < atol # either on head or intersection here already
            node.next.data.intersection = true
            return node.next
        elseif d1 < d2
            break
        end
        node = node.next
    end
    info = PointInfo(point, true, nothing, NONE)
    insert!(node, info)
end

function link_intersections!(
        inter1::Node{<:PointInfo}, 
        inter2::Node{<:PointInfo}, 
        edge1::Segment2D, 
        edge2::Segment2D; 
        atol::Float64=1e-6
    )
    point = inter1.data.point
    on_edge1 = (is_same_point(point, edge1[1]; atol=atol) || is_same_point(point, edge1[2]; atol=atol))
    head2_on_edge = is_same_point(point, edge2[2]; atol=atol)
    tail2_on_edge = is_same_point(point, edge2[1]; atol=atol)
    on_edge2 = head2_on_edge || tail2_on_edge
    if on_edge1 && on_edge2
        # intersect at common vertix
        prev1 = inter1.prev.data.point
        next1 = inter1.next.data.point
        edge1_prev = (prev1, point)
        edge1_next = (point, next1)
        next2 = inter2.next.data.point
        prev2 = inter2.prev.data.point
        # case where only one is true: \|/__   
        next2_in_1 = in_half_plane(next2, edge1_prev, false; on_border_is_inside=false) ||
                     in_half_plane(next2, edge1_next, false; on_border_is_inside=false)
        prev2_in_1 = in_half_plane(prev2, edge1_prev, false; on_border_is_inside=false) ||
                     in_half_plane(prev2, edge1_next, false; on_border_is_inside=false)
        next2_on_edge1, prev2_on_edge1 = has_edge_overlap(point, prev1, next1, prev2, next2)
        share_plane = has_plane_overlap(point, prev1, next1, prev2, next2)
        share_edge = next2_on_edge1 || prev2_on_edge1
        if (!share_plane && !share_edge)
            set_vertix_intercept!(inter1, inter2) # lone outer vertix 
        elseif xor(next2_in_1, prev2_in_1) # entry/exit point of edges/regions
            set_exit!(inter1, inter2, !next2_in_1 && prev2_in_1)
        elseif xor(next2_on_edge1, prev2_on_edge1) # share one edge
            set_exit!(inter1, inter2, !next2_on_edge1 && prev2_on_edge1)
        else # vertix between edges or lone inner/outer vertix 
            set_vertix_intercept!(inter1, inter2)
        end
    elseif head2_on_edge # edge2 hitting edge
        tail_in_1 = in_half_plane(edge2[1], edge1, false)
        set_exit!(inter1, inter2, tail_in_1)
    elseif tail2_on_edge # edge2 leaving edge
        head_in_1 = in_half_plane(edge2[2], edge1, false)
        set_exit!(inter1, inter2, !head_in_1)  
    else # cross or edge1 hitting/leaving edge
        exiting_1_to_2 = in_half_plane(edge2[1], edge1, false)
        set_exit!(inter1, inter2, exiting_1_to_2)
    end
    if is_vertix_intercept(inter2) # bounces off (case 2+3) or cycles back (case 2/3+4)
        inter1.data.type = VERTIX
        inter2.data.type = VERTIX
    end
end

function set_exit!(inter1::Node{<:PointInfo}, inter2::Node{<:PointInfo}, exiting_1_to_2::Bool)
    if exiting_1_to_2 # edge of polygon 2 is exiting poylgon 1 (intersection region)
        inter2.data.link = inter1 # so move from polygon2 to polygon1
        inter2.data.type = EXIT
    else
        inter1.data.link = inter2
        inter2.data.type = ENTRY
    end 
end

function set_vertix_intercept!(inter1::Node{<:PointInfo}, inter2::Node{<:PointInfo})
    inter1.data.link = inter2
    inter2.data.link = inter1
    inter2.data.type = VERTIX
    inter1.data.type = VERTIX
end

function walk_linked_lists(polygon::DoublyLinkedList{PointInfo{T}}) where T
    regions = Vector{Point2D{T}}[]
    node = polygon.head
    visited = Set{Point2D{T}}()
    while !isnothing(node)
        if !(node.data.point in visited) && node.data.intersection
            if node.data.type == ENTRY
                loop, visited_in_loop = walk_loop(node)
                if !isempty(visited_in_loop)
                    push!(visited, visited_in_loop...)
                end
                push!(regions, loop)
            end
        end
        node = node.next == polygon.head ? nothing : node.next
    end
    regions, visited
end

function walk_loop(start::Node{PointInfo{T}}) where T
    loop = Point2D{T}[]
    visited = Set{Point2D{T}}()
    push!(loop, start.data.point)
    push!(visited, start.data.point)
    node = start.next
    from_link = false # for debugging purposes
    while (node != start) && (node.data.point != start.data.point)
        push!(loop, node.data.point)
        if (node.data.point in visited) && 
            (node.data.type != VERTIX) &&  # many edges can hit the same vertix.
            (node.prev.data.point != node.data.point) # exception for repeated nodes
            @warn "Cycle detected: start node: $start; repeated node: $node"
            break
        end
        push!(visited, node.data.point)
        if !isnothing(node.data.link) && node.data.type != VERTIX 
            from_link = true
            node = node.data.link.next
        else
            from_link = false
            node = node.next
        end
    end
    loop, visited
end

function is_vertix_intercept(node::Node{<:PointInfo})
    node2 = node.data.link
    !isnothing(node2) && !isnothing(node2.data.link) && node2.data.link == node
end

function get_unvisited_intercepts(polygon::DoublyLinkedList{PointInfo{T}}, visited::Set{Point2D{T}}) where T
    regions = Vector{Point2D{T}}[]
    node = polygon.head
    while !isnothing(node)
        if node.data.intersection && !(node.data.point in visited) 
            push!(regions, [node.data.point])
        end
        node = node.next == polygon.head ? nothing : node.next
    end
    regions
end

function has_edge_overlap(vertix::Point2D, prev1::Point2D, next1::Point2D, prev2::Point2D, next2::Point2D)
    tail_edge1 = (prev1, vertix)
    head_edge1 = (vertix, next1)
    tail_edge2 = (prev2, vertix)
    head_edge2 = (vertix, next2)
    mid_prev1 = segment_midpoint(tail_edge1)
    mid_head1 = segment_midpoint(head_edge1)
    mid_prev2 = segment_midpoint(tail_edge2)
    mid_next2 = segment_midpoint(head_edge2)
    # check both because one edge might be shorter
    prev2_on_head1 = on_segment(mid_prev2, head_edge1) || on_segment(mid_head1, tail_edge2)
    next2_on_tail1 = on_segment(mid_next2, tail_edge1) || on_segment(mid_prev1, head_edge2)
    prev2_on_tail1 = on_segment(mid_prev2, tail_edge1) || on_segment(mid_prev1, tail_edge2)
    next2_on_head1 = on_segment(mid_head1, head_edge2) || on_segment(mid_next2, head_edge1)
    prev2_on_edge1 = prev2_on_head1 || prev2_on_tail1
    next2_on_edge1 = next2_on_tail1 || next2_on_head1
    next2_on_edge1, prev2_on_edge1
end

function segment_midpoint(segment::Segment2D)
    x = (segment[1][1] + segment[2][1])/2
    y = (segment[1][2] + segment[2][2])/2
    (x, y)
end

function has_plane_overlap(vertix::Point2D, prev1::Point2D, next1::Point2D, prev2::Point2D, next2::Point2D)
    # share a common plane: >> or << might only intersect at the vertix or share edges or share a whole region
    # do not share a common plane: >< might intersect at a vertix or share edges 
    edge1_prev = (prev1, vertix)
    mid1 = ((next1[1] + prev1[1])/2, (next1[2] + prev1[2])/2)
    mid2 = ((next2[1] + prev2[1])/2, (next2[2] + prev2[2])/2)
    in_plane12 = in_half_plane(mid2, edge1_prev, false; on_border_is_inside=false)
    in_plane11 = in_half_plane(mid1, edge1_prev, false; on_border_is_inside=false)
    in_plane11 == in_plane12
end