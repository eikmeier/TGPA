"""
'triangle generalized_preferential_attachment_graph'
===========================================
Generate an instance of a triangle generalized preferential attachment graph which
follows the Eikmeier, Gleich description. This is an undirected graph
that is generated as follows:
- Start with an empty graph (or a clique of size k0, if you insist)
- Add n vertices where at each time step one of three events occurs:
A new node is added with probability p, along with two edges to existing nodes;
A new wedge between two existing nodes is added with probability r,
Three new nodes with two edges between them is added with probability 1 - p - r

Functions
---------
The following functions are synonyms
- 'triangle_generalized_preferential_attachment_graph'
- 'tgpa_graph'
and
- 'triangle_generalized_preferential_attachment_edges!'
- 'tgpa_edges!'
The computational functions are
- 'tgpa_graph(n,p,r,k0)' Generate a TGPA graph with n total nodes.
    This returns a MatrixNetwork type. Note that this process allows self loops.
The edge functions are
-   'tgpa_edges!(n,p,r,edges,n0)' Add new edges to an existing set, by taking
    n0 time steps. Edges are added in one of three ways: From a new node to
    an existing node with probability p, between three existing nodes with
    probability r, between three new nodes with probability 1-p-r

Input
-----
- 'n': the number of nodes in the final graph.
- 'p': The probability of a node event, p must be a constant.
- 'r': The probability of a wedge event, r must be a constant. p+r <=1
- 'k0': The number of nodes in the starting clique.
- 'edges': A list of edges to be manipulated in the process of generating
          new edges.
Output
------
- A matrix network type for the triangle generalized preferential attachment graph.
- 'edges': An updated list of edges.
Example:
triangle_generalized_preferential_attachment_graph(100,1/3,1/2)
"""
:triangle_generalized_preferential_attachment_graph, :tgpa_graph,
:triangle_generalized_preferential_attachment_edges!, :tgpa_edges!

using Printf #For Julia 1.0
using MatrixNetworks

function triangle_generalized_preferential_attachment_graph(
    n::Int,p::Float64,r::Float64,k0::Int)
    k0 >= 0 || throw(ArgumentError(@sprintf("k0=%i must be >= 0",n)))
    n >= k0 || throw(ArgumentError(@sprintf("n=%i must be >= k0=%i",n,k0)))
    0<=p<=1 || throw(ArgumentError(@sprintf("p=%0.3f must be between 0 and 1",p)))
    0<=r<=1 || throw(ArgumentError(@sprintf("r=%0.3f must be between 0 and 1",r)))
    p+r <= 1 || throw(ArgumentError(@sprintf("(p=%0.3f)+(r=%0.3f) must be <=1",p,r)))
    edges = Vector{Tuple{Int,Int}}()
    # add the clique
    for i=1:k0
        for j=1:i-1
            push!(edges, (i,j))
            push!(edges, (j,i))
        end
    end
    return MatrixNetwork(triangle_generalized_preferential_attachment_edges!(n,p,r,edges,k0),n)
end

function _edges_to_neighbors(edges, k0)
    neighbors = Vector{Vector{Int}}(undef, k0)
    for i=1:k0
        neighbors[i] = Vector{Int}()
    end
    for e in edges
        push!(neighbors[e[1]], e[2])
    end
    return neighbors
end

function triangle_generalized_preferential_attachment_edges!(
    n::Int,p::Float64,r::Float64,edges::Vector{Tuple{Int,Int}},k0::Int)
    i = k0
    # for this code, we need to have random neighbor access, so we'll keep
    neighbors = _edges_to_neighbors(edges, k0)
    while i < n
        #generate a random value between 0 and 1
        x = rand()
        if x < p #node event, self loops permitted
            if length(edges) > 0
                v = rand(edges)[1]
                v_neighbors = neighbors[v]
            else
                v = i+1
                v_neighbors = [];
            end
            #edge between new node and v
            push!(edges, (i+1, v))
            push!(edges, (v, i+1))
            #update neighbor vector
            push!(neighbors,Vector{Int}())
            push!(neighbors[i+1],v)
            push!(neighbors[v],i+1)

            #add another edge between new node and a neighbor of v
            if length(v_neighbors) > 0
                v2 = rand(v_neighbors)
                push!(edges, (i+1, v2))
                push!(edges, (v2, i+1))
                push!(neighbors[i+1],v2)
                push!(neighbors[v2],i+1)
            end
            i = i+1;
        elseif x < p+r #wedge event, self loops permitted
            if length(edges) > 0 #only allow this event if the graph is non-empty
                v1 = rand(edges)[1]
                v1_neighbors = neighbors[v1]
                v2 = rand(edges)[1]
                push!(edges, (v1, v2))
                push!(edges, (v2, v1))
                push!(neighbors[v1],v2)
                push!(neighbors[v2],v1)

                #find a neighbor of v1, and add another edge between v1 and that node
                v3 = rand(v1_neighbors)
                push!(edges, (v1, v3))
                push!(edges, (v3, v1))
                push!(neighbors[v1],v3)
                push!(neighbors[v3],v1)
            end
        else #component event
            if i+3 <= n #only allow this step if there is room for three more nodes
                push!(edges, (i+1, i+2))
                push!(edges, (i+2, i+1))
                push!(edges, (i+1, i+3))
                push!(edges, (i+3, i+1))
                push!(neighbors,Vector{Int}()) #three new nodes
                push!(neighbors,Vector{Int}())
                push!(neighbors,Vector{Int}())
                push!(neighbors[i+1],i+2)
                push!(neighbors[i+1],i+3)
                push!(neighbors[i+2],i+1)
                push!(neighbors[i+3],i+1)
                i = i+3;
            end
        end
    end
    return edges
end



tgpa_graph = triangle_generalized_preferential_attachment_graph
tgpa_edges! = triangle_generalized_preferential_attachment_edges!
