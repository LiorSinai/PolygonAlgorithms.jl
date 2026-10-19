"""
    weiler_atherton_algorithm(polygon1, polygon2; atol=default_atol)

The Weiler-Atherton polygon clipping algorithm.
Returns regions, edges and single points of intersection. 
Only returns the larger type if one is within another e.g. an edge is also part of a region.

It runs in `O(nm)` time where `n` and `m` are the number of vertices of `polygon1` and `polygon2` respectively.
Use `intersect_convex` for convex polygons for an `O(n+m)` algorithm.

Description: operates at a point level. Walks from point to point along `polygon2`. 
It starts recording loops at "entry points" - crossings from `polygon2` to `polygon1` - and 
stops recording when it gets back to the same entry point. Then it continues walking along `polygon2`
until it reaches the start point.

Limitations
1. This version does not cater for holes.
2. It can fail completely for self-intersecting polygons.

For a more general algorithm see the Martinez-Rueda polygon clipping algorithm.
"""
function weiler_atherton_algorithm(
    polygon1::Path2D{T}, polygon2::Path2D{T}
    ; atol::AbstractFloat=default_atol
    ) where T
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
    list1 = convert_to_linked_list(polygon1)
    list2 = convert_to_linked_list(polygon2)
    find_and_insert_intersections!(list1, list2; atol=atol)
    regions, visited = walk_linked_lists(list2)
    # special cases
    if isempty(regions)
        p1 = get_first_non_intersection_point(list1)
        p2 = get_first_non_intersection_point(list2)
        if contains(polygon2, p1; atol=atol)
            return [polygon1]
        elseif contains(polygon1, p2; atol=atol)
            return [polygon2]
        end
    end
    not_visited = get_unvisited_intercepts(list2, visited)
    push!(regions, not_visited...)
    regions
end

"""
    IntersectionType

Classification of a point in the Weiler-Atherton algorithm, from the perspective of `polygon2`:
- `NONE`: not an intersection.
- `ENTRY`: `polygon2` enters `polygon1` here.
- `EXIT`: `polygon2` exits `polygon1` here.
- `VERTEX`: the polygons touch at this point but `polygon2` does not cross into or out of `polygon1`.
"""
@enum IntersectionType NONE ENTRY EXIT VERTEX

"""
    PointEvent(point, intersection, link, type)

The data for a node in the Weiler-Atherton linked lists.
`intersection` is true if the point is an intersection between the polygons.
`link` is the corresponding node on the other polygon's linked list, which is followed when walking a region.
`type` is an `IntersectionType`.
"""
mutable struct PointEvent{T}
    point::Point2D{T}
    intersection::Bool
    link::Union{Nothing,Node{PointEvent{T}}}
    type::IntersectionType
end

"""
    convert_to_linked_list(polygon)

Convert a polygon to a circular `DoublyLinkedList` of `PointEvent`s with no intersections.
"""
function convert_to_linked_list(polygon::Path2D{T}) where T
    data = PointEvent{T}(polygon[1], false, nothing, NONE)
    head = Node(data)
    list = DoublyLinkedList(head)
    for point in polygon[2:end]
        data = PointEvent{T}(point, false, nothing, NONE)
        push!(list, data)
    end
    list
end

"""
    get_first_non_intersection_point(polygon::DoublyLinkedList)

The first point in the linked list which is not an intersection.
If all points are intersections, the last point is returned.
"""
function get_first_non_intersection_point(polygon::DoublyLinkedList{<:PointEvent})
    start = polygon.head
    node = start
    while node.data.intersection && node.next != start
        node = node.next
    end
    node.data.point
end

"""
    find_and_insert_intersections!(polygon1, polygon2; atol=default_atol)

Intersect every edge of `polygon1` with every edge of `polygon2` in `O(nm)` time.
Each intersection point is inserted into both linked lists in order along the edge
"""
function find_and_insert_intersections!(
        polygon1::DoublyLinkedList{PointEvent{T}}, 
        polygon2::DoublyLinkedList{PointEvent{T}};
        atol::AbstractFloat=default_atol
        ) where T
    ## collect original nodes before mutating in place
    vec1 = collect_nodes(polygon1)
    vec2 = collect_nodes(polygon2)
    n1 = length(vec1)
    n2 = length(vec2)
    ## find intersections
    for (j, node2) in enumerate(vec2)
        next2 = vec2[j % n2 + 1]
        edge2 = (node2.data.point, next2.data.point)
        for (i, node1) in enumerate(vec1)
            next1 = vec1[i % n1 + 1]
            edge1 = (node1.data.point, next1.data.point) 
            p = intersect_geometry(edge1, edge2; atol=atol)
            if !isnothing(p)
                i1 = insert_intersection_in_order!(p, node1, next1; atol=atol)
                i2 = insert_intersection_in_order!(p, node2, next2; atol=atol)
                link_intersections!(i1, i2, edge1, edge2; atol=atol)
            end
        end
    end
end

"""
    insert_intersection_in_order!(point, tail, head; atol=default_atol)

Insert an intersection `point` between the nodes `tail` and `head`, which were originally adjacent
but may now have other intersections between them.
The point is inserted in order of distance from `tail`.

If the point is within `atol` of an existing node, that node is marked as an intersection and returned
instead of inserting a new node.
"""
function insert_intersection_in_order!(
    point::Point2D, tail::Node{<:PointEvent}, head::Node{<:PointEvent}
    ; atol::AbstractFloat=default_atol
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
    info = PointEvent(point, true, nothing, NONE)
    insert!(node, info)
end

"""
    link_intersections!(inter1, inter2, edge1, edge2; atol=default_atol)

Link the intersection nodes `inter1` (on `polygon1`) and `inter2` (on `polygon2`)
and classify the intersection as an `ENTRY`, `EXIT` or `VERTEX`.

Cases:
1. The intersection is at a vertex of both polygons. Use the angles of the adjacent edges to determine if
    `polygon2` crosses `polygon1` or only touches it.
2. `edge2` ends on `edge1`.
3. `edge2` starts on `edge1`.
4. The edges cross, or `edge1` starts or ends on `edge2`.

Assumes both polygons are clockwise.
"""
function link_intersections!(
        inter1::Node{<:PointEvent}, 
        inter2::Node{<:PointEvent}, 
        edge1::Segment2D, 
        edge2::Segment2D; 
        atol::AbstractFloat=default_atol
    )
    point = inter1.data.point
    on_edge1 = (is_same_point(point, edge1[1]; atol=atol) || is_same_point(point, edge1[2]; atol=atol))
    head2_on_edge = is_same_point(point, edge2[2]; atol=atol)
    tail2_on_edge = is_same_point(point, edge2[1]; atol=atol)
    on_edge2 = head2_on_edge || tail2_on_edge
    if on_edge1 && on_edge2 # case 1: intersect at common vertex
        prev1 = inter1.prev.data.point
        next1 = inter1.next.data.point
        next2 = inter2.next.data.point
        prev2 = inter2.prev.data.point
        next2_on_edge1, prev2_on_edge1 = has_edge_overlap(point, prev1, next1, prev2, next2)
        tail2_in_1, head2_in_1 = in_plane((prev1, point, next1), prev2, next2)
        tail2_in_1 = tail2_in_1 || prev2_on_edge1
        head2_in_1 = head2_in_1 || next2_on_edge1
        if head2_in_1 == tail2_in_1
            # vertex between edges or lone inner/outer vertex 
            set_vertex_intercept!(inter1, inter2)  
        else 
            set_link!(inter1, inter2, head2_in_1)
        end
    elseif head2_on_edge # case 2: edge2 hitting edge
        tail_in_1 = in_half_plane(edge1, edge2[1], false)
        set_link!(inter1, inter2, !tail_in_1)
    elseif tail2_on_edge # case 3: edge2 leaving edge
        head_in_1 = in_half_plane(edge1, edge2[2], false)
        set_link!(inter1, inter2, head_in_1)  
    else # case 4: cross or edge1 hitting/leaving edge
        entering_1_from_2 = in_half_plane(edge1, edge2[2], false)
        set_link!(inter1, inter2, entering_1_from_2)
    end
    if is_vertex_intercept(inter2) # bounces off (case 2+3) or cycles back (case 2/3+4)
        inter1.data.type = VERTEX
        inter2.data.type = VERTEX
    end
    inter1, inter2
end

function set_link!(inter1::Node{<:PointEvent}, inter2::Node{<:PointEvent}, entering_1_from_2::Bool)
    if entering_1_from_2 # edge of polygon 2 is entering polygon 1
        inter1.data.link = inter2
        inter2.data.type = ENTRY
    else 
        inter2.data.link = inter1
        inter2.data.type = EXIT
    end 
end

function set_vertex_intercept!(inter1::Node{<:PointEvent}, inter2::Node{<:PointEvent})
    inter1.data.link = inter2
    inter2.data.link = inter1
    inter2.data.type = VERTEX
    inter1.data.type = VERTEX
end

"""
    walk_linked_lists(polygon2)

Walk along `polygon2` and record a region with `walk_loop` from each unvisited `ENTRY` point.

Returns the regions and a `PointSet` of all visited points.
"""
function walk_linked_lists(polygon::DoublyLinkedList{PointEvent{T}}) where T
    regions = Vector{Point2D{T}}[]
    node = polygon.head
    visited = PointSet(Dict{Point2D{T}, Nothing}(), 6)
    while !isnothing(node)
        if !(node.data.point in visited) && node.data.type == ENTRY
            loop, visited_in_loop = walk_loop(node)
            if !isempty(visited_in_loop)
                visited_points = [p.data.point for p in visited_in_loop]
                push!(visited, visited_points...)
            end
            push!(regions, loop)
        end
        node = node.next == polygon.head ? nothing : node.next
    end
    regions, visited
end

"""
    walk_loop(start)

Record a single region of intersection, starting from the `ENTRY` node `start`.
At each linked intersection (excluding `VERTEX`s), jump to the other polygon's list.
Stops when it returns to `start`.

Returns the points of the region and the set of visited nodes.
"""
function walk_loop(start::Node{PointEvent{T}}) where T
    loop = Point2D{T}[]
    visited = Set{Node{PointEvent{T}}}()
    push!(loop, start.data.point)
    push!(visited, start)
    node = start.next
    from_link = false # for debugging purposes
    while (node != start) && (node.data.point != start.data.point)
        push!(loop, node.data.point)
        if (node in visited) && # can go to the same point on different nodes!
            (node.data.type != VERTEX) &&  # many edges can hit the same vertex
            (node.prev.data.point != node.data.point) # exception for repeated nodes
            @warn "Cycle detected: start node: $start; repeated node: $node"
            break
        end
        push!(visited, node)
        if !isnothing(node.data.link) && node.data.type != VERTEX 
            if (node.data.link == start)
                break
            end
            from_link = true
            node = node.data.link.next
        else
            from_link = false
            node = node.next
        end
    end
    loop, visited
end

function is_vertex_intercept(node::Node{<:PointEvent})
    node2 = node.data.link
    !isnothing(node2) && !isnothing(node2.data.link) && node2.data.link == node
end

"""
    get_unvisited_intercepts(polygon, visited)

Return each intersection point on `polygon` that was not visited while walking the regions, as a single point region.
These are points where the polygons touch but do not overlap.
"""
function get_unvisited_intercepts(polygon::DoublyLinkedList{PointEvent{T}}, visited::PointSet{T}) where T
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

"""
    has_edge_overlap(vertex, prev1, next1, prev2, next2)

At a common `vertex`, check if the edges of `polygon2` to `next2` and from `prev2` overlap with
either edge of `polygon1`. Midpoints are used to compare edges of different lengths.

Returns `(next2_on_edge1, prev2_on_edge1)`.
"""
function has_edge_overlap(vertex::Point2D, prev1::Point2D, next1::Point2D, prev2::Point2D, next2::Point2D)
    tail_edge1 = (prev1, vertex)
    head_edge1 = (vertex, next1)
    tail_edge2 = (prev2, vertex)
    head_edge2 = (vertex, next2)
    mid_tail1 = segment_midpoint(tail_edge1)
    mid_head1 = segment_midpoint(head_edge1)
    mid_tail2 = segment_midpoint(tail_edge2)
    mid_head2 = segment_midpoint(head_edge2)
    # check both because one edge might be shorter
    prev2_on_tail1 = on_segment(mid_tail2, tail_edge1) || on_segment(mid_tail1, tail_edge2)
    prev2_on_head1 = on_segment(mid_tail2, head_edge1) || on_segment(mid_head1, tail_edge2)
    next2_on_tail1 = on_segment(mid_head2, tail_edge1) || on_segment(mid_tail1, head_edge2)
    next2_on_head1 = on_segment(mid_head1, head_edge2) || on_segment(mid_head2, head_edge1)
    # check both because unsure of directions
    prev2_on_edge1 = prev2_on_head1 || prev2_on_tail1
    next2_on_edge1 = next2_on_tail1 || next2_on_head1
    next2_on_edge1, prev2_on_edge1
end

function segment_midpoint(segment::Segment2D)
    x = (segment[1][1] + segment[2][1])/2
    y = (segment[1][2] + segment[2][2])/2
    (x, y)
end

"""
    in_plane(edge, point1, point2)

For the two edges `edge[1] → edge[2] → edge[3]` meeting at the vertex `edge[2]`, check if `point1` and `point2`
lie in the interior angle between the edges. Assumes the interior is clockwise from the edges.

Angles are compared instead of half-planes because a point can be in one half-plane but outside the polygon.

Returns `(point1_inside, point2_inside)`.
"""
function in_plane(edge::NTuple{3, Point2D}, point1::Point2D, point2::Point2D)
    # assume half-plane is clockwise of the edge
    # if in half-planes of edge[1:2] && edge[2:3] then definitely in the plane.
    # But if in one half-plane it can be either in or out. 
    # A better strategy is to check angles instead.
    vertex = edge[2]
    angle_tail = atan_pos(edge[1][2] - vertex[2], edge[1][1] - vertex[1])
    angle_head = atan_pos(edge[3][2] - vertex[2], edge[3][1] - vertex[1])
    angle1 = atan_pos(point1[2] - vertex[2], point1[1] - vertex[1])
    angle2 = atan_pos(point2[2] - vertex[2], point2[1] - vertex[1])
    if angle_tail < angle_head
        tail2_in_1 = (angle_tail < angle1 < angle_head) 
        head2_in_1 = (angle_tail < angle2 < angle_head) 
    else
        tail2_in_1 = (angle1 > angle_tail) || (angle1 < angle_head)
        head2_in_1 = (angle2 > angle_tail) || (angle2 < angle_head)
    end
    tail2_in_1, head2_in_1
end

"""
    atan_pos(y, x)

`atan(y, x)` mapped to the positive range `[0, 2π)`.
"""
function atan_pos(y::Real, x::Real)
    angle = atan(y, x)
    angle < 0 ? angle + 2π : angle
end
