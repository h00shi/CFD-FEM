/*
 * UnstMeshBcFaces.test.cpp
 *
 *  Created on: Oct 9, 2015
 *      Author: rabbit
 */

#include "Mesh/UnstMeshBcFaces.h"
#include "IO/UnstMeshReaderNKBGrid.h"
#include "Mesh/UnstMeshElements.h"
#include "gtest/gtest.h"

TEST(BcFaces, Read){

  UnstMeshReader* mesh_reader = new UnstMeshReaderNKBGrid("Square.grid");
  UnstMeshBcFaces mesh_bc_faces(*mesh_reader);
  UnstMeshElements mesh_elements(*mesh_reader);
  mesh_bc_faces.FormBcFace2Element(mesh_elements);

  EXPECT_EQ(16, mesh_bc_faces.get_nbc_face());
  EXPECT_EQ(4, mesh_bc_faces.get_nbc_id());

  const Array1D<intT>& bc_face2elem = mesh_bc_faces.get_bc_face2element();
  const Array1D<intT>& bc_face_id = mesh_bc_faces.get_bc_face_id();
  const Array1D<ElementTopology::face_types>& bc_face_type = mesh_bc_faces.get_bc_face_type();
  EXPECT_EQ(0, bc_face2elem(0));
  EXPECT_EQ(25, bc_face2elem(4));
  EXPECT_EQ(23, bc_face2elem(9));

  EXPECT_EQ(0, bc_face_id(0));
  EXPECT_EQ(1, bc_face_id(4));
  EXPECT_EQ(2, bc_face_id(9));

  const Array1D<intT>& bc_local_face = mesh_bc_faces.get_bc_local_face();

  EXPECT_EQ(0, bc_local_face(0));
  EXPECT_EQ(1, bc_local_face(9));
  EXPECT_EQ(2, bc_local_face(15));

 for(intT i = 0; i < mesh_bc_faces.get_nbc_face(); i++){
    EXPECT_EQ(ElementTopology::face_types::FACE_BAR, bc_face_type(i));
  }
}
TEST(BcFaces, Read2){

  UnstMeshReader* mesh_reader = new UnstMeshReaderNKBGrid("Square.grid");
  UnstMeshElements mesh_elements(*mesh_reader);

  UnstMeshBcFaces mesh_bc_faces(*mesh_reader, mesh_elements);
  EXPECT_EQ(16, mesh_bc_faces.get_nbc_face());
  EXPECT_EQ(4, mesh_bc_faces.get_nbc_id());

  const Array1D<intT>& bc_face2elem = mesh_bc_faces.get_bc_face2element();
  const Array1D<intT>& bc_face_id = mesh_bc_faces.get_bc_face_id();
  const Array1D<ElementTopology::face_types>& bc_face_type = mesh_bc_faces.get_bc_face_type();
  EXPECT_EQ(0, bc_face2elem(0));
  EXPECT_EQ(25, bc_face2elem(4));
  EXPECT_EQ(23, bc_face2elem(9));

  EXPECT_EQ(0, bc_face_id(0));
  EXPECT_EQ(1, bc_face_id(4));
  EXPECT_EQ(2, bc_face_id(9));

  const Array1D<intT>& bc_local_face = mesh_bc_faces.get_bc_local_face();

  EXPECT_EQ(0, bc_local_face(0));
  EXPECT_EQ(1, bc_local_face(9));
  EXPECT_EQ(2, bc_local_face(15));

 for(intT i = 0; i < mesh_bc_faces.get_nbc_face(); i++){
    EXPECT_EQ(ElementTopology::face_types::FACE_BAR, bc_face_type(i));
  }
}
