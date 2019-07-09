# graph plots
##
using Plots
#font = Plots.font("Helvetica", 18)
using SparseArrays
using LinearAlgebra
using MatrixNetworks
using PyCall
using ColorTypes
@pyimport igraph
include("TGPA_generator.jl")
include("TGPA_generator_v2.jl")
include("Holme_generator.jl")

##

r = RGB(36/255,155/255,223/255)
#bg=RGBA(r,0.1)
# This function uses tricks we have learned so far to plot a graph
function graphplot(A,xy,loc,short_name)
    f = plot(leg=false, axis = false,grid=false)
    ei,ej,w = findnz(triu(A))
    lx = [xy[ei,1]';xy[ej,1]';NaN*ones(1,length(ei))]
    ly = [xy[ei,2]';xy[ej,2]';NaN*ones(1,length(ei))]
    for i = 1:length(w)
        plot!(f,lx[:,i],ly[:,i],color = :grey, linewidth = 0.8)
    end
    scatter!(f,xy[:,1],xy[:,2],color = 1,marker = (8,stroke(1,:gray)))
    title!(short_name,titlefont=font)
    save_to = loc*".pdf"
    savefig(save_to)
    return f
end

function igraph_layout(A::SparseMatrixCSC{T}, layoutname::AbstractString="lgl") where T
    (ei,ej) = findnz(A)
    edgelist = [(ei[i]-1,ej[i]-1) for i = 1:length(ei)]
    nverts = size(A)
    G = igraph.Graph(nverts, edges=edgelist, directed=true)
    layoutname = "fr"
    xy = G[:layout](layoutname)
    xy = [ Float64(xy[i][j]) for i in 1:length(xy),  j in 1:length(xy[1])]
end





## TGPA
A = tgpa_graph(50,0.8,0.1,0);
A = sparse(A);
A = A+A'; #symmetric
A = min.(A, 1); #simple graph
A = convert(SparseMatrixCSC{Int64,Int64},A);
xy = igraph_layout(A)
f = graphplot(A,xy,"/Users/eikmeier/Dropbox/TGPA/WWW/images/TGPA_example_network","TGPA(p,r,q)")
#f = graphplot(A,xy,"testing","TGPA")

## TGPA_v2
A = tgpa_graph_v2(50,0.85,0,2);
A = sparse(A);
A = A+A'; #symmetric
A = min.(A, 1); #simple graph
A = convert(SparseMatrixCSC{Int64,Int64},A);
xy = igraph_layout(A)
f = graphplot(A,xy,"/Users/eikmeier/Dropbox/TGPA/WWW/images/TGPA2_example_network","TGPA(p,q)")
#f = graphplot(A,xy,"testing","TGPA")

## GPA
A = gpa_graph(50,0.8,0.1,2);
A = sparse(A)
A = A+A';
A = min.(A,1);
A = convert(SparseMatrixCSC{Int64,Int64},A);

xy = igraph_layout(A)
f = graphplot(A,xy, "/Users/eikmeier/Dropbox/TGPA/WWW/images/GPA_example_network","GPA")


## Holme TGPA
A = holme_graph(50,2,2,0.8);
A = sparse(A);
A = A+A'; #symmetric
A = min.(A, 1); #simple graph
A = convert(SparseMatrixCSC{Int64,Int64},A);
xy = igraph_layout(A)
#f = graphplot(A,xy,"/Users/eikmeier/Dropbox/TGPA/WWW/images/Holme_example_network","Holme")
