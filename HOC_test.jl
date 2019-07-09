#This code will test the HigherOrderClustering Coefficients of the TGPA graph
#
#Note that the HigherOrderClustering package has not been updated for Julia 1.0
#I tested this code on Julia 0.6 and it works

using MatrixNetworks
using HigherOrderClustering
include("TGPA_generator.jl")
include("TGPA_generator_v2.jl")


p = [0.9,0.8,0.7,0.6,0.5,0.4,0.3]
r = [0.09, 0.15,0.25,0.35,0.45,0.55,0.65]
for i = 1:length(p)
    pii = p[i]
    rii = r[i]
    n = 1000
    A = triangle_generalized_preferential_attachment_graph(n,pii,rii,0);
    #A = triangle_generalized_preferential_attachment_graph_v2(n,pii,0,4);
    A = sparse(A);
    A = A+A'; #symmetric
    A = min.(A, 1); #simple graph
    A = convert(SparseMatrixCSC{Int64,Int64},A);
    ccfs = higher_order_ccfs(A,2);
    @show ccfs.global_hoccf
    @show ccfs.avg_hoccf
end
