"""
'Holme_Graph'
===========================================
Generate an instance of a graph which
follows the Holme description. This is an undirected graph
that is generated as follows:
- Start with a clique of size k0
- Add n vertices where at each time step a new node is added:
with m new triangles
Functions
---------
The computational functions are
- 'holme_graph(n,k0,m)' Generate a TGPA graph with n total nodes.
    This returns a MatrixNetwork type. Note that this process allows self loops.
The edge functions are
-   'holme_edges!(n,edges,n0)' Add new edges to an existing set, by taking
    n0 time steps. Edges are added in one of two ways: From a new node to
    an existing node with probability p, between three existing nodes with
    probability 1-p

Input
-----
- 'n': the number of nodes in the final graph.
- 'k0': The number of nodes in the starting clique.
- 'm': The number of new triangles to add at each time step
- 'edges': A list of edges to be manipulated in the process of generating
          new edges.
Output
------
- A matrix network type for the triangle generalized preferential attachment graph.
- 'edges': An updated list of edges.
Example:
triangle_generalized_preferential_attachment_graph(100,1/2,3)
"""
#:triangle_generalized_preferential_attachment_graph, :tgpa_graph,
#:triangle_generalized_preferential_attachment_edges!, :tgpa_edges!

using Printf #For Julia 1.0
using MatrixNetworks
using Distributions

function holme_graph(n::Int,k0::Int,m::Int,p::Float64)
    k0 >= 0 || throw(ArgumentError(@sprintf("k0=%i must be >= 0",n)))
    n >= k0 || throw(ArgumentError(@sprintf("n=%i must be >= k0=%i",n,k0)))
    m >= 1 || throw(ArgumentError(@sprintf("m=%i must be >= 1",m)))
    m <= max(n,1) || throw(ArgumentError(@sprintf("m=%i must be <= n=%i",m,n)))
    0<= p <= 1 || throw(ArgumentError(@sprintf("p=%i must be between 0 and 1",p)))
    edges = Vector{Tuple{Int,Int}}()
    # add the clique
    for i=1:k0
        for j=1:i-1
            push!(edges, (i,j))
            push!(edges, (j,i))
        end
    end
    return MatrixNetwork(holme_edges!(n,edges,k0,m,p),n)
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


function holme_edges!(
    n::Int,edges::Vector{Tuple{Int,Int}},k0::Int,m::Int,p::Float64)
    i = k0

    # for this code, we need to have random neighbor access, so we'll keep
    neighbors = _edges_to_neighbors(edges, k0)

    while i < n
        push!(neighbors,Vector{Int}()) #add a new node to the graph
        for j = 1:m #do m PA steps
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
            push!(neighbors[i+1],v)
            push!(neighbors[v],i+1)

            #add another edge between new node and a neighbor of v, with probability $p$
            x = rand(1)[1]
            if x > p
                if length(v_neighbors) > 0
                    v2 = rand(v_neighbors)
                    push!(edges, (i+1, v2))
                    push!(edges, (v2, i+1))
                    push!(neighbors[i+1],v2)
                    push!(neighbors[v2],i+1)
                end
            end #triad closure
        end #end of PA steps
        i = i+1;
    end

    return edges
end
