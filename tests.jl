using Base.Test
include("TGPA_generator.jl")
include("TGPA_generator_v2.jl")

@test typeof(tgpa_graph(10,0.5,0.3,0)) == MatrixNetwork{Bool}
@test typeof(tgpa_edges!(4,.75,.25,[(1,2)],2)) == Vector{Tuple{Int,Int}}
@test is_empty(tgpa_graph(0,0.0,0.0,0)) == true
@test is_undirected(tgpa_graph(10, .2, .8,0))
@test_throws ArgumentError tgpa_graph(-1,.3,.4,0)
@test_throws ArgumentError tgpa_graph(10,.3,.4,12)
@test_throws ArgumentError tgpa_graph(8,.8,.4,0)
@test size(tgpa_graph(12,.5,.4,0),1) == 12
@test size(tgpa_graph(21,.1,.1,0),1) == 21


@test typeof(tgpa_graph_v2(10,0.5,3,3)) == MatrixNetwork{Bool}
@test typeof(tgpa_edges_v2!(4,.75,2,[(1,2)],1)) == Vector{Tuple{Int,Int}}
@test is_empty(tgpa_graph_v2(0,0.0,0,1)) == true
@test is_undirected(tgpa_graph_v2(10, .2, 0,2))
@test_throws ArgumentError tgpa_graph_v2(-1,.3,0,1)
@test_throws ArgumentError tgpa_graph_v2(10,.3,12,3)
@test_throws ArgumentError tgpa_graph_v2(8,1.1,0,3)
@test size(tgpa_graph_v2(12,.5,0,2),1) == 12
@test size(tgpa_graph_v2(21,.1,0,3),1) == 21
