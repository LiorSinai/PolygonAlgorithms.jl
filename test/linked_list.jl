@testset "DoublyLinkedList" begin

using PolygonAlgorithms: DoublyLinkedList, Node, collect_nodes

@testset "list basic" begin
    head = Node(1.0)
    list = DoublyLinkedList(head)
    push!(list, 2.0)
    push!(list, 3.0)

    node2 = head.next
    @test !isnothing(node2)
    @test node2.data == 2.0
    @test node2.prev == head
    
    node3 = node2.next
    @test !isnothing(node3)
    @test node3.data == 3.0
    @test node3.prev == node2
    @test node3.next == head

    nodes = collect_nodes(list)
    @test nodes == [head, node2, node3]

    expected = "1.0, 2.0, 3.0"
    result = join(list, ", ")
    @test result == expected
end

@testset "list insert" begin
    head = Node(1.0)
    list = DoublyLinkedList(head)
    node2 = push!(list, 2.0)
    node3 = push!(list, 3.0)

    node1_5 = insert!(head, 1.5)
    @test node1_5.data == 1.5
    @test node1_5 == head.next
    @test node1_5.prev == head
    @test node1_5 != node2
    @test node1_5 == node2.prev

    node2_5 = insert!(node2, 2.5)
    @test node2_5.data == 2.5
    @test node2_5 == node2.next
    @test node2_5 == node3.prev

    node3_5 = insert!(node3, 3.5)
    @test node3_5 == node3.next
    @test node3_5.prev == node3
    @test node3_5.next == head

    expected = "1.0, 1.5, 2.0, 2.5, 3.0, 3.5"
    result = join(list, ", ")
    @test result == expected
end

end