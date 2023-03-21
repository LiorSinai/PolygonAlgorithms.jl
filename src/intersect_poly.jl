"""
    intersect_geometry(polygon1, polygon2)

Returns the points of regions of intersections of two general polygons. 
Uses the Weiler-Atherton algorithm which fails for self-intersecting polygons. 
This version also does not cater for holes.
For a more general algorithm, see the Martinez-Rueda polygon clipping algorithm.

Runs in `O(nm)` time where `n` and `m` are the number of vertices of polygon1 and polygon2 respectively.

Returns:
    1. Regions of intersection.
    2. Edges of intersection.
    3. Single points of intersection.

If one type is within another e.g. an edge is also part of a region, only returns the larger type.

For convex polygons use `intertect_convex` for an `O(n+m)` algorithm.
"""
function intersect_geometry(polygon1::Polygon{T}, polygon2::Polygon{T}) where T
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
    vertix_intercepts = get_lone_vertix_intercepts(list2, visited)
    push!(regions, vertix_intercepts...)
    regions
end

@enum IntersectionType NONE ENTRY EXIT VERTIX

mutable struct PointInfo{T}
    point::Point2D{T}
    intersection::Bool
    link::Union{Nothing,Node{PointInfo{T}}}
    type::IntersectionType
end

function generate_list(polygon::Polygon{T}) where T
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
        polygon2::DoublyLinkedList{PointInfo{T}}
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
                i1 = insert_intersection!(p, node1, next1; id=1)
                i2 = insert_intersection!(p, node2, next2; id=2)
                link_intersections!(i1, i2, edge1, edge2)
            end
        end
    end
end

function insert_intersection!(point::Point2D, node::Node{<:PointInfo}, next::Node{<:PointInfo}; atol=1e-6, id=1)
    if is_same_point(point, node.data.point)
        i1 = node
        i1.data.intersection = true
    elseif is_same_point(point, next.data.point)
        i1 = next
        i1.data.intersection = true
    else
        i1 = insert_intersection_in_order!(node, point; id=id, atol=atol)
    end
    i1
end

function insert_intersection_in_order!(node::Node{<:PointInfo}, point::Point2D; atol=1e-6, id=1)
    start = node
    while node.next.data.intersection && node.next != start
        d1 = norm2(point, node.data.point)
        d2 = norm2(node.next.data.point, node.data.point)
        if d1 < d2
            break
        end
        node = node.next
    end
    if is_same_point(point, node.data.point; atol=atol) # (vertix) intersection here already
        node.data.intersection = true
        return node
    end
    info = PointInfo(point, true, nothing, NONE)
    insert!(node, info)
end

function link_intersections!(
        inter1::Node{<:PointInfo}, 
        inter2::Node{<:PointInfo}, 
        edge1::Segment2D, 
        edge2::Segment2D; 
        atol=1e-6
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
        head_in_1 = in_half_plane(next2, edge1_prev; on_border_is_inside=false) ||
                    in_half_plane(next2, edge1_next; on_border_is_inside=false)
        tail_in_1 = in_half_plane(prev2, edge1_prev; on_border_is_inside=false) ||
                    in_half_plane(prev2, edge1_next; on_border_is_inside=false)
        if head_in_1 != tail_in_1 # entry/exit point
            set_exit!(inter1, inter2, !head_in_1 && tail_in_1)
        else 
            # check if there is an segment overlap
            share_head_edge, share_tail_edge = share_edges(point, prev1, next1, prev2, next2)
            if share_head_edge && !share_tail_edge
                set_exit!(inter1, inter2, false)
                return
            elseif !share_head_edge && share_tail_edge
                set_exit!(inter1, inter2, true)
                return
            end
            inter1.data.link = inter2
            inter2.data.link = inter1
            inter2.data.type = VERTIX
            inter1.data.type = VERTIX
        end
        return
    elseif head2_on_edge # edge2 hitting edge
        tail_in_1 = in_half_plane(edge2[1], edge1)
        set_exit!(inter1, inter2, tail_in_1) # note: if true it might bounce back in and will only be a pseudo-exit
    elseif tail2_on_edge # edge2 leaving edge
        head_in_1 = in_half_plane(edge2[2], edge1)
        set_exit!(inter1, inter2, !head_in_1)  
    else # cross or edge1 hitting/leaving edge
        exiting_1_to_2 = in_half_plane(edge2[1], edge1; on_border_is_inside=false)
        set_exit!(inter1, inter2, exiting_1_to_2)
    end
    if is_vertix_intercept(inter2) # bounces off (case 2+3) or cycles back (case 2+3+4)
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

function walk_linked_lists(polygon::DoublyLinkedList{PointInfo{T}}) where T
    regions = Vector{Point2D{T}}[]
    node = polygon.head
    visited = Set{Node{PointInfo{T}}}()
    while !isnothing(node)
        if !(node in visited) && node.data.intersection
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
    visited = Set{Node{PointInfo{T}}}()
    push!(loop, start.data.point)
    push!(visited, start)
    node = start.next
    from_link = false
    while (node != start) && (node.data.point != start.data.point)
        if (node in visited) && (node.data.type != VERTIX) # many edges can hit the same vertix
            @warn "Cycle detected: start node: $start; repeated node: $node"
            break
        end
        push!(visited, node)
        if !from_link
            push!(loop, node.data.point)
        end
        if node.data.type == VERTIX 
            push!(visited, node.data.link)
            from_link = false
            node = node.next
        elseif !isnothing(node.data.link)
            from_link = true
            node = node.data.link
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

function get_lone_vertix_intercepts(polygon::DoublyLinkedList{PointInfo{T}}, visited::Set{Node{PointInfo{T}}}) where T
    regions = Vector{Point2D{T}}[]
    node = polygon.head
    while !isnothing(node)
        if node.data.intersection && !(node in visited) && node.data.type == VERTIX
            push!(regions, [node.data.point])
        end
        node = node.next == polygon.head ? nothing : node.next
    end
    regions
end

function share_edges(vertix::Point2D, prev1::Point2D, next1::Point2D, prev2::Point2D, next2::Point2D)
    orientation1 = get_orientation(prev1, vertix, next1)
    orientation2 = get_orientation(prev2, vertix, next2)
    tail_edge1 = orientation1 == orientation2 ? (prev1, vertix) : (vertix, next1)
    head_edge1 = orientation1 == orientation2 ? (vertix, next1) : (prev1, vertix)
    mid_prev1 = (((tail_edge1[1][1] + tail_edge1[2][1])/2), ((tail_edge1[1][2] + tail_edge1[2][2])/2))
    mid_next1 = (((head_edge1[1][1] + head_edge1[2][1])/2), ((head_edge1[1][2] + head_edge1[2][2])/2))
    mid_prev2 = (((prev2[1] + vertix[1])/2), ((prev2[2] + vertix[2])/2))
    mid_next2 = (((next2[1] + vertix[1])/2), ((next2[2] + vertix[2])/2))
    # check both because one edge might be shorter
    share_tail_edge = on_segment(mid_prev2, tail_edge1) || on_segment(mid_prev1, (prev2, vertix))
    share_head_edge = on_segment(mid_next2, head_edge1) || on_segment(mid_next1, (vertix, next2))
    share_head_edge, share_tail_edge
end
