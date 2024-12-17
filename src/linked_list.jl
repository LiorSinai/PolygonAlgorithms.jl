# see DataStrctures.jl for a full implementation

mutable struct Node{T}
    data::T
    prev::Node{T}
    next::Node{T}
    function Node(data::T) where T
        node = new{T}(data)
        node.next = node
        node.prev = node
        return node
    end
end

function Node(data::T, prev::Node{T}, next::Node{T}) where T
    node = Node(data)
    node.prev = prev
    node.next = next
    node
end

"""
    insert!(node, data)

Insert a new node `Node(data)` after `node`.
"""
function Base.insert!(node::Node{T}, data::T) where {T}
    middle = Node(data, node, node.next)
    node.next.prev = middle
    node.next = middle
    middle
end

function Base.show(io::IO, node::Node)
    print(io, typeof(node), '(')
    print(io, node.data)
    print(io, ')')
end

mutable struct DoublyLinkedList{T}
    head::Union{Node{T}, Nothing}
end

function Base.push!(list::DoublyLinkedList{T}, data::T) where {T}
    if isnothing(list.head)
        list.head = Node(data)
        return list.head
    end
    head = list.head
    tail = head.prev
    node = Node(data, tail, head)
    head.prev = node
    tail.next = node
    node
end

Base.iterate(list::DoublyLinkedList) = (list.head.data, list.head.next)
Base.iterate(list::DoublyLinkedList, state::Node) = state == list.head ? nothing : (state.data, state.next)

function Base.show(io::IO, list::DoublyLinkedList)
    print(io, typeof(list), "(")
    print(io, isnothing(list.head) ? nothing : join(list, ", "))
    print(io, ")")
end

function Base.show(io::IO, m::MIME"text/plain", list::DoublyLinkedList)
    if isnothing(list.head)
        print(io, typeof(list), "(", nothing, ")")
    else
        print(io, typeof(list), "(\n   ")
        print(io, join(list, ",\n   "))
        print(io, "\n)")
    end
end

function Base.length(list::DoublyLinkedList)
    count = 1
    next = list.head.next
    while next != list.head
        count += 1
        next = next.next
    end
    count
end

function collect_nodes(list::DoublyLinkedList{T}) where T
    # A standard collect(list) will return the data, not the nodes themselves
    vec = Node{T}[]
    node = list.head
    push!(vec, node)
    while node.next != list.head
        node = node.next
        push!(vec, node)
    end
    vec
end
