#include "Mesh/UnstMeshElements.h"
#include "Mesh/UnstMeshEdges.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/List2D.h"
#include "gtest/gtest.h"
#include <set>
TEST(UnstMeshEdges, Construct)
{
  UnstMeshReader* mesh_reader = new UnstMeshReaderNKBGrid("Square.grid");
  UnstMeshElements mesh_elem(*mesh_reader);

  UnstMeshEdges mesh_edges(mesh_elem);
  
  const Array2D<intT>& edge2node = mesh_edges.get_edge2node();
  EXPECT_EQ(7, edge2node(9,0));
  EXPECT_EQ(3, edge2node(9,1));

  EXPECT_EQ(17, edge2node(37,0));
  EXPECT_EQ(13, edge2node(37,1));

  const List2D<intT>& elem2edge = mesh_edges.get_element2edge();
  EXPECT_EQ(32, elem2edge.get_lead_size());
  EXPECT_EQ(96, elem2edge.get_total_size());

  EXPECT_EQ(3,  elem2edge.get_ncol(6));
  EXPECT_EQ(12, elem2edge(6,0));
  EXPECT_EQ(13, elem2edge(6,1));
  EXPECT_EQ(14, elem2edge(6,2));

  EXPECT_EQ(3,  elem2edge.get_ncol(12));
  EXPECT_EQ(23, elem2edge(12,0));
  EXPECT_EQ(24, elem2edge(12,1));
  EXPECT_EQ(11, elem2edge(12,2));

  const List2D<intT>& edge2elem = mesh_edges.get_edge2element();
  EXPECT_EQ(56, edge2elem.get_lead_size());
  EXPECT_EQ(96, edge2elem.get_total_size());
  EXPECT_EQ(3,edge2elem(8,0));
  EXPECT_EQ(4,edge2elem(8,1));
  EXPECT_EQ(1,edge2elem.get_ncol(6));
  EXPECT_EQ(2,edge2elem(6,0));
  EXPECT_EQ(10,edge2elem(21,0));
  EXPECT_EQ(11,edge2elem(21,1));

  delete mesh_reader;
}

