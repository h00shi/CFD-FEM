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
  
  EXPECT_EQ(N, graph->get_nnode());
  EXPECT_EQ(N-1,graph->get_nedge());

  for(intT i = 0; i < graph->get_nnode(); i++){
    for(intT j = graph->NeighborBegin(i); 
	j < graph->NeighborEnd(i); j++){
      EXPECT_EQ(adj2(i,j), graph->GetNeighbor(i,j));
    }
  }
  
}

TEST(Graph, Constructor2) {
  //---> Setup
  const intT N = 10;
  List2D<intT> adj(N,3*(N-2) + 2*2);
  List2D<intT> adj2(N,3*(N-2) + 2*2);
  Array2D<intT> edge2node(N-1,2);
  Array1D<intT> ncol(N);
  Graph* graph;
  ncol(0) = 2;
  for(intT i = 1; i < N-1; i++){
    ncol(i) = 3;
  }
  ncol(N - 1) = 2;
  adj.set_ncol(ncol);
  adj2.set_ncol(ncol);
 
  adj(0,0) = 0;
  adj(0,1) = 1;
  adj2(0,0) = 0;
  adj2(0,1) = 1;
  for(intT i = 1; i < N-1; i++){
    adj(i,0) = i - 1;
    adj(i,1) = i;
    adj(i,2) = i + 1;
    adj2(i,0) = i - 1;
    adj2(i,1) = i;
    adj2(i,2) = i + 1;
  }
  adj(N-1,0) = N-2;
  adj(N-1,1) = N-1;
  adj2(N-1,0) = N-2;
  adj2(N-1,1) = N-1;

  for(intT i = 0; i < N-1; i++){
    edge2node(i,0) = i;
    edge2node(i,1) = i + 1;
  }

  //---> Allocate Graph and transfer adj and edge2node
  graph = new Graph(adj);
  
  //---> Check that old structures are gone
  EXPECT_EQ(0,adj.get_lead_size());
  EXPECT_EQ(0,adj.get_total_size());
  
  EXPECT_EQ(N, graph->get_nnode());
  EXPECT_EQ(N-1,graph->get_nedge());

  for(intT i = 0; i < graph->get_nnode(); i++){
    for(intT j = graph->NeighborBegin(i); j < graph->NeighborEnd(i); j++){
      EXPECT_EQ(adj2(i,j), graph->GetNeighbor(i,j));
    }
  }
  
  for(intT e = 0; e < graph->get_nedge(); e++){
    EXPECT_EQ(edge2node(e,0), graph->get_GraphEdge2Node()(e,1));
    EXPECT_EQ(edge2node(e,1), graph->get_GraphEdge2Node()(e,0));
  }
  
}
TEST(Graph, ExtraPointers)
{
  //---> Setup
  const intT N = 10;
  List2D<intT> adj(N,3*(N-2) + 2*2);
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
  //---> Allocate Graph and transfer adj and edge2node
  graph = new Graph(adj);
  
  for(intT i = 0; i < graph->get_nnode(); i++){
    for(intT j = 0; j < graph->get_GraphAdj().get_ncol(i); j++){
      if(graph->get_GraphAdj()(i,j) == i){
        EXPECT_EQ(j, graph->get_GraphNode2SelfIndex()(i));
      }
    }
  }
  
  for(intT i = 0; i < graph->get_nedge(); i++){
    intT nl = graph->get_GraphEdge2Node()(i,0);
    intT nr = graph->get_GraphEdge2Node()(i,1);
    
    //---> Left node check
    for(intT j = 0; j < graph->get_GraphAdj().get_ncol(nl); j++){
      if(graph->get_GraphAdj()(nl,j) == nr){ 
        EXPECT_EQ(j,graph->get_GraphEdge2AdjIndex()(i,0));
      }
    }
    for(intT j = 0; j < graph->get_GraphAdj().get_ncol(nr); j++){
      if(graph->get_GraphAdj()(nr,j) == nl){ 
        EXPECT_EQ(j,graph->get_GraphEdge2AdjIndex()(i,1));
      } 
    }

  }

  std::cout << graph->get_GraphAdj() << std::endl;
}
