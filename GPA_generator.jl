"""
'generalized_preferential_attachment_graph'
===========================================
Generate an instance of a generalized preferential attachment graph which
follows the Avin,Lotker,Nahum,Peleg description. This is an undirected graph
that is generated as follows:
- Start with a k0-node clique
- Add n - k0 vertices where at each time step one of three events occurs: A new
node is added with probability p, a new edge between two existing nodes is added
with probability r, two new nodes with an edge between them is added with
probability 1 - p - r
Functions
---------
The following functions are synonyms
- 'generalized_preferential_attachment_graph'
- 'gpa_graph'
and
- 'generalized_preferential_attachment_edges!'
- 'gpa_edges!'
The computational functions are
- 'gpa_graph(n,p,r,k0)' Generate a GPA graph with a k0 clique and n total nodes.
    This returns a MatrixNetwork type
- 'gpa_graph(n,p,r,k0,Val{true})' Generate a GPA graph with a k0 clique and
    n total nodes, allowing self-loops. This returns a MatrixNetwork type
The edge functions are
-   'gpa_edges!(n,p,r,edges,n0)' Add new edges to an existing set, by taking
    n0 time steps. Edges are added in one of three ways: From a new node to
    an existing node with probability p, between two existing nodes with
    probability r, between two new nodes with probability 1-p-r
-   'gpa_edges!(n,p,r,edges,n0,Val{true})' Add new edges to an existing set, by
    taking n0 time steps. Edges are added in one of three ways: From a new node
    to an existing node with probability p, between two existing nodes with
    probability r (allowing self-loops), between two new nodes with probability
    1-p-r
Input
-----
- 'n': the number of nodes in the final graph.
- 'p': The probability of a node event, p must be a constant.
- 'r': The probability of an edge event, r must be a constant. p+r <=1
- 'k0': the number of nodes in the starting clique.
- 'Val{true}': Include this parameter if self-loops are allowed. Default is false
- 'edges': A list of edges to be manipulated in the process of generating
          new edges.
Output
------
- A matrix network type for the generalized preferential attachment graph.
- 'edges': An updated list of edges.
Example:
generalized_preferential_attachment_graph(100,1/3,1/2,2)
"""
:generalized_preferential_attachment_graph, :gpa_graph,
:generalized_preferential_attachment_edges!, :gpa_edges!

using IterTools

generalized_preferential_attachment_graph(n::Int,p::Float64,r::Float64,k0::Int) =
    generalized_preferential_attachment_graph(n,p,r,k0,Val{false})
generalized_preferential_attachment_edges!(n::Int,p::Float64,r::Float64,edges::Vector{Tuple{Int,Int}},n0::Int) =
    generalized_preferential_attachment_edges!(n,p,r,edges,n0,Val{false})

function generalized_preferential_attachment_graph(
    n::Int,p::Float64,r::Float64,k0::Int,::Type{Val{true}})
    k0 >= 0 || throw(ArgumentError(@sprintf("k0=%i must be non-negative",k0)))
    n >= k0 || throw(ArgumentError(@sprintf("n=%i must be >= k0=%i",n,k0)))
    0<=p<=1 || throw(ArgumentError(@sprintf("p=%0.3f must be between 0 and 1",p)))
    0<=r<=1 || throw(ArgumentError(@sprintf("r=%0.3f must be between 0 and 1",r)))
    p+r <= 1 || throw(ArgumentError(@sprintf("(p=%0.3f)+(r=%0.3f) must be <=1",p,r)))
    edges = Vector{Tuple{Int,Int}}()
    #add the clique
    for i = 1:k0
        for j = 1:i-1
            push!(edges,(i,j))
            push!(edges, (j,i))
        end
    end
    return MatrixNetwork(generalized_preferential_attachment_edges!(n,p,r,edges,k0,Val{true}),n)
end

function generalized_preferential_attachment_graph(
    n::Int,p::Float64,r::Float64,k0::Int,::Type{Val{false}})
    k0 >= 0 || throw(ArgumentError(@sprintf("k0=%i must be non-negative",k0)))
    n >= k0 || throw(ArgumentError(@sprintf("n=%i must be >= k0=%i",n,k0)))
    0<=p<=1 || throw(ArgumentError(@sprintf("p=%0.3f must be between 0 and 1",p)))
    0<=r<=1 || throw(ArgumentError(@sprintf("r=%0.3f must be between 0 and 1",r)))
    p+r <= 1 || throw(ArgumentError(@sprintf("(p=%0.3f)+(r=%0.3f) must be <=1",p,r)))
    edges = Vector{Tuple{Int,Int}}()
    #add the clique
    for i = 1:k0
        for j = 1:i-1
            push!(edges,(i,j))
            push!(edges, (j,i))
        end
    end
    return MatrixNetwork(generalized_preferential_attachment_edges!(n,p,r,edges,k0,Val{false}),n)
end

function generalized_preferential_attachment_edges!(
    n::Int,p::Float64,r::Float64,edges::Vector{Tuple{Int,Int}},n0::Int,::Type{Val{true}})
    i = n0
    while i < n
        #generate a random value between 0 and 1
        x = rand()
        if x < p #node event
            v = rand(edges)[1]
            push!(edges, (i+1,v[1]))
            push!(edges, (v[1], i+1))
            i = i+1;
        elseif x < p+r #edge event, self loops permitted
            v1 = rand(edges)[1]
            v2 = rand(edges)[1]
            push!(edges, (v1, v2))
            push!(edges, (v2, v1))
        else #component event
            if i+2 <= n #only allow this step if there is room for two more nodes
                push!(edges, (i+1, i+2))
                push!(edges, (i+2, i+1))
                i = i+2;
            end
        end
    end
    return edges
end

function _check_for_two_distinct_nodes(edges::Vector{Tuple{Int,Int}})
    length(edges) > 0 || throw(ArgumentError("requires at least one edge"))
    firstnode = edges[1][1]
    return any(IterTools.imap(x -> firstnode != x[1] || firstnode != x[2], edges))
end

function generalized_preferential_attachment_edges!(
    n::Int,p::Float64,r::Float64,edges::Vector{Tuple{Int,Int}},n0::Int,::Type{Val{false}})
    i = n0
    if i >= n
        return edges
    end

    if !_check_for_two_distinct_nodes(edges::Vector{Tuple{Int,Int}})
        throw(ArgumentError("The starting graph must have at least two distinct nodes"))
    end

    while i < n
        #generate a random value between 0 and 1
        x = rand()
        if x < p #node event
            v = rand(edges)[1]
            push!(edges, (i+1,v[1]))
            push!(edges, (v[1], i+1))
            i = i+1;
        elseif x < p+r #edge event, no self-loops permitted
            v1 = rand(edges)[1]
            v2 = rand(edges)[1]
            while (v1 == v2 && i != 1) #i != 1 because we want more than 1 node for this to work
                v1 = rand(edges)[1]
                v2 = rand(edges)[1]
            end
            push!(edges, (v1, v2))
            push!(edges, (v2, v1))
        else #component event
            if i+2 <= n #only allow this step if there is room for two more nodes
                push!(edges, (i+1, i+2))
                push!(edges, (i+2, i+1))
                i = i+2;
            end
        end
    end
    return edges
end

gpa_graph = generalized_preferential_attachment_graph
gpa_edges! = generalized_preferential_attachment_edges!
