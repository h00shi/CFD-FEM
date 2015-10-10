#include "Mesh/CGMesh.h"
#include "gtest/gtest.h"
#include "IO/UnstMeshReaderNKBGrid.h"
TEST(UnstMesh, InitfromFile) 
{
  UnstMeshReaderNKBGrid grid_reader("Square.grid");
  CGMesh mesh(grid_reader);
 
  mesh.Diagnostic(std::cout);

}

TEST(UnstMesh, Coords) {
  UnstMeshReaderNKBGrid grid_reader("Square.grid");
  CGMesh mesh(grid_reader);
  const Array2D<realT>& x = mesh.get_MeshGeom().get_x();

  EXPECT_DOUBLE_EQ(-1.0, x(0,0));
  EXPECT_DOUBLE_EQ(-0.5, x(0,1));
  EXPECT_DOUBLE_EQ(0.0, x(13,0));
  EXPECT_DOUBLE_EQ(.25, x(13,1));
  EXPECT_DOUBLE_EQ(1.0, x(24,0));
  EXPECT_DOUBLE_EQ(0.5, x(24,1));
}

TEST(UnstMesh, Element2node) {
  UnstMeshReaderNKBGrid grid_reader("Square.grid");
  CGMesh mesh(grid_reader);
  const List2D<intT>& element2node =
      mesh.get_MeshElements().get_element2node();

  EXPECT_EQ(5, element2node(8,0));
  EXPECT_EQ(10, element2node(8,1));
  EXPECT_EQ(6, element2node(8,2));

  EXPECT_EQ(12, element2node(20, 0));
  EXPECT_EQ(17, element2node(20,1));
  EXPECT_EQ(13, element2node(20,2));

}

TEST(UnstMesh, Node2Element){
  UnstMeshReaderNKBGrid grid_reader("Square.grid");
  CGMesh mesh(grid_reader);
  const List2D<intT>& node2element =
      mesh.get_MeshElements().get_node2element();

  EXPECT_EQ(1, node2element(6,0));
  EXPECT_EQ(10, node2element(6,5));

  EXPECT_EQ(13, node2element(12, 2));

}

TEST(UnstMesh, BcConnectivity){
  UnstMeshReaderNKBGrid grid_reader("Square.grid");
  CGMesh mesh(grid_reader);
  const Array1D<intT>& bc_face2elem =
      mesh.get_MeshBcFaces().get_bc_face2element();
  const Array1D<intT>& bc_face_id =
      mesh.get_MeshBcFaces().get_bc_face_id();

  EXPECT_EQ(0, bc_face2elem(0));
  EXPECT_EQ(25, bc_face2elem(4));
  EXPECT_EQ(23, bc_face2elem(9));

  EXPECT_EQ(0, bc_face_id(0));
  EXPECT_EQ(1, bc_face_id(4));
  EXPECT_EQ(2, bc_face_id(9));


}

