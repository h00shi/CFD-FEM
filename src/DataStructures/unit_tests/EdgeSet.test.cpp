#include "gtest/gtest.h"
#include "DataStructures/EdgeSet.h"

TEST(EdgeSet, Initialize)
{
  intT nnode = 5;
  intT nedge = 7;
  //---> Now make up edges from nodes
  EdgeSet edges(nnode,nedge);
  edges.InsertEdge(0,1);
  edges.InsertEdge(1,2);
  edges.InsertEdge(2,0);
  
  edges.InsertEdge(2,1);
  edges.InsertEdge(1,3);
  edges.InsertEdge(3,2);
  
  edges.InsertEdge(1,4);
  edges.InsertEdge(4,3);
  edges.InsertEdge(3,1);

  Array2D<intT> edge2node = edges.RemoveEdge2Node();
  EXPECT_EQ(nedge, edge2node.get_size(0));

  EXPECT_EQ(0,edge2node(0,1));
  EXPECT_EQ(1,edge2node(0,0));
  
  EXPECT_EQ(3,edge2node(3,0));
  EXPECT_EQ(1,edge2node(3,1));
}
TEST(EdgeSet, Node2Edge)
{
  intT nnode = 5;
  intT nedge = 7;
  //---> Now make up edges from nodes
  EdgeSet edges(nnode,nedge);
  edges.InsertEdge(0,1);
  edges.InsertEdge(1,2);
  edges.InsertEdge(2,0);

  edges.InsertEdge(2,1);
  edges.InsertEdge(1,3);
  edges.InsertEdge(3,2);

  edges.InsertEdge(1,4);
  edges.InsertEdge(4,3);
  edges.InsertEdge(3,1);

  List2D<intT> node2edge = edges.FormNode2Edge();

  EXPECT_EQ(0,node2edge(0,0));
  EXPECT_EQ(2,node2edge(0,1));

  EXPECT_EQ(0,node2edge(1,0));
  EXPECT_EQ(1,node2edge(1,1));
  EXPECT_EQ(3,node2edge(1,2));
  EXPECT_EQ(5,node2edge(1,3));

  EXPECT_EQ(1,node2edge(2,0));
  EXPECT_EQ(2,node2edge(2,1));
  EXPECT_EQ(4,node2edge(2,2));


  EXPECT_EQ(3,node2edge(3,0));
  EXPECT_EQ(4,node2edge(3,1));
  EXPECT_EQ(6,node2edge(3,2));

  EXPECT_EQ(5,node2edge(4,0));
  EXPECT_EQ(6,node2edge(4,1));
}
TEST(EdgeSet, Adj)
{
  intT nnode = 5;
  intT nedge = 7;
  //---> Now make up edges from nodes
  EdgeSet edges(nnode,nedge);
  edges.InsertEdge(0,1);
  edges.InsertEdge(1,2);
  edges.InsertEdge(2,0);

  edges.InsertEdge(2,1);
  edges.InsertEdge(1,3);
  edges.InsertEdge(3,2);

  edges.InsertEdge(1,4);
  edges.InsertEdge(4,3);
  edges.InsertEdge(3,1);

  List2D<intT> adj = edges.FormNodeAdj();

  EXPECT_EQ(0,adj(0,0));
  EXPECT_EQ(1,adj(0,1));
  EXPECT_EQ(2,adj(0,2));

  EXPECT_EQ(0,adj(1,0));
  EXPECT_EQ(1,adj(1,1));
  EXPECT_EQ(2,adj(1,2));
  EXPECT_EQ(3,adj(1,3));
  EXPECT_EQ(4,adj(1,4));

  EXPECT_EQ(0,adj(2,0));
  EXPECT_EQ(1,adj(2,1));
  EXPECT_EQ(2,adj(2,2));
  EXPECT_EQ(3,adj(2,3));

  EXPECT_EQ(1,adj(3,0));
  EXPECT_EQ(2,adj(3,1));
  EXPECT_EQ(3,adj(3,2));
  EXPECT_EQ(4,adj(3,3));

  EXPECT_EQ(1,adj(4,0));
  EXPECT_EQ(3,adj(4,1));
  EXPECT_EQ(4,adj(4,2));
}
