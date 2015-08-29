#include "UnstMesh.h"
#include "UnstMeshWriter.h"
#include "UnstMeshWriterCGNS.h"
#include "Array1D.h"
#include "Array2D.h"
#include "List2D.h"

#include "gtest/gtest.h"

static const realT small = 5.0e-15;
TEST(UnstMesh, InitfromFile) {
  UnstMesh mesh("Square.grid", "Grid-NKB");
  mesh.Diagnostic();
  
  //---> Test writing
  UnstMeshWriter* writer = new UnstMeshWriterCGNS(mesh);
  writer->Write("Square");
  
}

TEST(UnstMesh, Coords) {
  UnstMesh mesh("Square.grid", "Grid-NKB");
  const Array2D<realT>& x = mesh.get_x();
  
  EXPECT_NEAR(-1.0, x(0,0), small);
  EXPECT_NEAR(-0.5, x(0,1), small);
  EXPECT_NEAR(0.0, x(13,0), small);
  EXPECT_NEAR(.25, x(13,1), small);
  EXPECT_NEAR(1.0, x(24,0), small);
  EXPECT_NEAR(0.5, x(24,1), small); 
}

TEST(UnstMesh, Element2node) {
  UnstMesh mesh("Square.grid", "Grid-NKB");
  const List2D<intT>& element2node = mesh.get_element2node();
  
  EXPECT_EQ(5, element2node(8,0));
  EXPECT_EQ(10, element2node(8,1));
  EXPECT_EQ(6, element2node(8,2));

  EXPECT_EQ(12, element2node(20, 0));
  EXPECT_EQ(17, element2node(20,1));
  EXPECT_EQ(13, element2node(20,2));

}

TEST(UnstMesh, Node2Element){
  UnstMesh mesh("Square.grid", "Grid-NKB");
  const List2D<intT>& node2element = mesh.get_node2element();

  EXPECT_EQ(1, node2element(6,0));
  EXPECT_EQ(10, node2element(6,5));

  EXPECT_EQ(13, node2element(12, 2));

}

TEST(UnstMesh, BcConnectivity){
  UnstMesh mesh("Square.grid", "Grid-NKB");
  const Array1D<intT>& bc_face2elem = mesh.get_bc_face2elem();
  const Array1D<intT>& bc_face_id = mesh.get_bc_face_id();
  
  EXPECT_EQ(0, bc_face2elem(0));
  EXPECT_EQ(25, bc_face2elem(4));
  EXPECT_EQ(23, bc_face2elem(9));

  EXPECT_EQ(0, bc_face_id(0));
  EXPECT_EQ(1, bc_face_id(4));
  EXPECT_EQ(2, bc_face_id(9));
  
  
}

TEST(UnstMesh, Adjacency){
  UnstMesh mesh("Square.grid", "Grid-NKB");
  const List2D<intT>& adj = mesh.get_adj();
}
