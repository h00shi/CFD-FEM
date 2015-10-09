#include "Mesh/UnstMeshElements.h"
#include "DataStructures/List2D.h"
#include "gtest/gtest.h"
TEST(UnstMeshGoem, Init) 
{
  UnstMeshReader* mesh_reader = new UnstMeshReaderNKBGrid("Square.grid");
  UnstMeshElements mesh_elem(*mesh_reader);
  
  const List2D<intT>& element2node = mesh_elem.get_element2node();
  const List2D<intT>& node2element = mesh_elem.get_node2element();

  EXPECT_EQ(5, element2node(8,0));
  EXPECT_EQ(10, element2node(8,1));
  EXPECT_EQ(6, element2node(8,2));

  EXPECT_EQ(12, element2node(20, 0));
  EXPECT_EQ(17, element2node(20,1));
  EXPECT_EQ(13, element2node(20,2));

  EXPECT_EQ(1, node2element(6,0));
  EXPECT_EQ(10, node2element(6,5));

  EXPECT_EQ(13, node2element(12, 2));

  delete mesh_reader;
}
