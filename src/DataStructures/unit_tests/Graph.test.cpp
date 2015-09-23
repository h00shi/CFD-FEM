#include "DataStructures/Graph.h"
#include "gtest/gtest.h"
#include "my_incl.h"
TEST(Graph, Constructor) {
  //---> Setup
  const intT N = 10;
  List2D<intT> adj(N,3*(N-2) + 2*2);
  Array2D<intT> edge2node(N-1,2);
  Array1D<intT> ncol(N);
  Graph* graph;
  ncol(0) = 2;
  for(intT i = 1; i < N-1; i++){
    ncol(i) = 3;
  }
  ncol(N - 1) = 2;
  adj.set_ncol(ncol);
 
  adj(0,0) = 0;
  adj(0,1) = 1;
  for(intT i = 1; i < N-1; i++){
    adj(i,0) = i - 1;
    adj(i,1) = i;
    adj(i,2) = i + 1;
  }

  adj(N-1,0) = N-2;
  adj(N-1,1) = N-1;
  
  for(intT i = 0; i < N-1; i++){
    edge2node(i,0) = i;
    edge2node(i,1) = i + 1;
  }

  //---> Allocate Graph and transfer adj and edge2node
  graph = new Graph(adj,edge2node);
  //---> Check that old structures are gone
  EXPECT_EQ(0,adj.get_lead_size());
  EXPECT_EQ(0,adj.get_total_size());
  EXPECT_EQ(0,edge2node.get_size(0));
  EXPECT_EQ(0,edge2node.get_size(1));
  
  const List2D<intT>& adj2 = graph->get_GraphAdj();
  const Array2D<intT>& e2n = graph->get_GraphEdge2Node();

  EXPECT_EQ(N, graph->get_nnode());
  EXPECT_EQ(N-1,graph->get_nedge());

  for(intT i = 0; i < graph->get_nnode(); i++){
    for(intT j = graph->NeighborBegin(i); 
	j < graph->NeighborEnd(i); j++){
      SystemModule::cout << graph->GetNeighbor(i,j) << " ";
    }
    SystemModule::cout << std::endl;
  }
  
  
}
