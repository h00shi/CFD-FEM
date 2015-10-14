/*
 * UnstMeshReaderGMSH.test.cpp
 *
 *  Created on: Oct 11, 2015
 *      Author: rabbit
 */
#include "gtest/gtest.h"
#include "IO/UnstMeshReaderGMSH.h"

TEST(GMSH_Reader, Read2D_ASCII)
{
  UnstMeshReaderGMSH mesh_reader("Circle.msh", "Circle.idmap");
  EXPECT_EQ(694, mesh_reader.ReadNnode());
  EXPECT_EQ(1287, mesh_reader.ReadNelement());
  EXPECT_EQ(100, mesh_reader.ReadNbcFace());
  EXPECT_EQ(2, mesh_reader.ReadNbcID());

  Array2D<realT> x = mesh_reader.ReadNodes();
  EXPECT_EQ(694, x.get_size(0));
  EXPECT_EQ(3, x.get_size(1));

  EXPECT_DOUBLE_EQ(-10.0, x(0,0));
  EXPECT_DOUBLE_EQ(-10.0, x(0,1));
  EXPECT_DOUBLE_EQ(0.0, x(0,2));

  EXPECT_DOUBLE_EQ(-7.99999999999753, x(6,0));
  EXPECT_DOUBLE_EQ(-10.0, x(6,1));
  EXPECT_DOUBLE_EQ(0.0, x(6,2));

  EXPECT_DOUBLE_EQ(0.8116612181008351, x(693,0));
  EXPECT_DOUBLE_EQ(-0.241928798846167, x(693,1));
  EXPECT_DOUBLE_EQ(0.0, x(693,2));

  List2D<intT> elem2node = mesh_reader.ReadElement2Node();
  EXPECT_EQ(1287, elem2node.get_lead_size());
  EXPECT_EQ(1287*3, elem2node.get_total_size());

  EXPECT_EQ(3, elem2node.get_ncol(0));
  EXPECT_EQ(167, elem2node(0,0));
  EXPECT_EQ(138, elem2node(0,1));
  EXPECT_EQ(683, elem2node(0,2));

  EXPECT_EQ(3, elem2node.get_ncol(757));
  EXPECT_EQ(58, elem2node(757,0));
  EXPECT_EQ(515, elem2node(757,1));
  EXPECT_EQ(59, elem2node(757,2));

  EXPECT_EQ(3, elem2node.get_ncol(1286));
  EXPECT_EQ(418, elem2node(1286,0));
  EXPECT_EQ(669, elem2node(1286,1));
  EXPECT_EQ(664, elem2node(1286,2));

  Array1D<ElementTopology::element_types> elem_type =
      mesh_reader.ReadElementType();

  EXPECT_EQ(mesh_reader.ReadNelement(), elem_type.get_size(0));
  for(intT i = 0; i < elem_type.get_size(0); i++){
    EXPECT_EQ(ElementTopology::element_types::TRI, elem_type(i));
  }

  Array1D<intT> elem_region = mesh_reader.ReadElementRegion();

  EXPECT_EQ(mesh_reader.ReadNelement(), elem_region.get_size(0));
  for(intT i = 0; i < elem_type.get_size(0); i++){
    EXPECT_EQ(0, elem_region(i));
  }
  List2D<intT> bc_face2node = mesh_reader.ReadBcFace2Node();
  EXPECT_EQ(100, bc_face2node.get_lead_size());
  EXPECT_EQ(100*2, bc_face2node.get_total_size());

  EXPECT_EQ(2, bc_face2node.get_ncol(0));
  EXPECT_EQ(0, bc_face2node(0,0));
  EXPECT_EQ(6, bc_face2node(0,1));

  EXPECT_EQ(2, bc_face2node.get_ncol(99));
  EXPECT_EQ(99, bc_face2node(99,0));
  EXPECT_EQ(4, bc_face2node(99,1));

  Array1D<intT> bc_id = mesh_reader.ReadBcID();
  EXPECT_EQ(0, bc_id(0));
  EXPECT_EQ(1, bc_id(40));

  Array1D<ElementTopology::face_types> bc_face_type = mesh_reader.ReadBcFaceType();
  for(intT i = 0; i < bc_face_type.get_size(0); i++){
    EXPECT_EQ(ElementTopology::face_types::FACE_BAR, bc_face_type(i));
  }
}
TEST(GMSH_Reader, Read3D_ASCII)
{
  UnstMeshReaderGMSH mesh_reader("Box.msh", "Box.idmap");
  EXPECT_EQ(292, mesh_reader.ReadNnode());
  EXPECT_EQ(1129, mesh_reader.ReadNelement());
  EXPECT_EQ(396, mesh_reader.ReadNbcFace());
  EXPECT_EQ(1, mesh_reader.ReadNbcID());

  Array2D<realT> x = mesh_reader.ReadNodes();
  EXPECT_EQ(292, x.get_size(0));
  EXPECT_EQ(3, x.get_size(1));

  EXPECT_DOUBLE_EQ(0.0, x(0,0));
  EXPECT_DOUBLE_EQ(0.0, x(0,1));
  EXPECT_DOUBLE_EQ(0.0, x(0,2));

  EXPECT_DOUBLE_EQ(1.0, x(6,0));
  EXPECT_DOUBLE_EQ(1.0, x(6,1));
  EXPECT_DOUBLE_EQ(1.0, x(6,2));

  EXPECT_DOUBLE_EQ(0.508718729019165, x(291,0));
  EXPECT_DOUBLE_EQ(0.2503562569618225, x(291,1));
  EXPECT_DOUBLE_EQ(0.6738438010215759, x(291,2));

  List2D<intT> elem2node = mesh_reader.ReadElement2Node();
  EXPECT_EQ(1129, elem2node.get_lead_size());
  EXPECT_EQ(1129*4, elem2node.get_total_size());

  EXPECT_EQ(4, elem2node.get_ncol(0));
  EXPECT_EQ(285, elem2node(0,0));
  EXPECT_EQ(268, elem2node(0,1));
  EXPECT_EQ(267, elem2node(0,2));
  EXPECT_EQ(207, elem2node(0,3));

  EXPECT_EQ(4, elem2node.get_ncol(360));
  EXPECT_EQ(88, elem2node(360,0));
  EXPECT_EQ(234, elem2node(360,1));
  EXPECT_EQ(83, elem2node(360,2));
  EXPECT_EQ(96, elem2node(360,3));

  EXPECT_EQ(4, elem2node.get_ncol(1128));
  EXPECT_EQ(186, elem2node(1128,0));
  EXPECT_EQ(220, elem2node(1128,1));
  EXPECT_EQ(189, elem2node(1128,2));
  EXPECT_EQ(194, elem2node(1128,3));
  Array1D<ElementTopology::element_types> elem_type =
      mesh_reader.ReadElementType();

  EXPECT_EQ(mesh_reader.ReadNelement(), elem_type.get_size(0));
  for(intT i = 0; i < elem_type.get_size(0); i++){
    EXPECT_EQ(ElementTopology::element_types::TET, elem_type(i));
  }

  Array1D<intT> elem_region = mesh_reader.ReadElementRegion();

  EXPECT_EQ(mesh_reader.ReadNelement(), elem_region.get_size(0));
  for(intT i = 0; i < elem_type.get_size(0); i++){
    EXPECT_EQ(0, elem_region(i));
  }
  List2D<intT> bc_face2node = mesh_reader.ReadBcFace2Node();
  EXPECT_EQ(396, bc_face2node.get_lead_size());
  EXPECT_EQ(396*3, bc_face2node.get_total_size());

  EXPECT_EQ(3, bc_face2node.get_ncol(0));
  EXPECT_EQ(58, bc_face2node(0,0));
  EXPECT_EQ(72, bc_face2node(0,1));
  EXPECT_EQ(73, bc_face2node(0,2));

  EXPECT_EQ(3, bc_face2node.get_ncol(395));
  EXPECT_EQ(179, bc_face2node(395,0));
  EXPECT_EQ(199, bc_face2node(395,1));
  EXPECT_EQ(192, bc_face2node(395,2));

  Array1D<intT> bc_id = mesh_reader.ReadBcID();
  EXPECT_EQ(0, bc_id(0));
  EXPECT_EQ(0, bc_id(40));

  Array1D<ElementTopology::face_types> bc_face_type = mesh_reader.ReadBcFaceType();
  for(intT i = 0; i < bc_face_type.get_size(0); i++){
    EXPECT_EQ(ElementTopology::face_types::FACE_TRI, bc_face_type(i));
  }

}
