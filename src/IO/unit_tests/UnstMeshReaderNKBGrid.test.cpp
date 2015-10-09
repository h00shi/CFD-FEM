#include "IO/UnstMeshReaderNKBGrid.h"
#include "gtest/gtest.h"
TEST(UnstMeshReaderNKBGrid, ReadIn)
{
  UnstMeshReaderNKBGrid grid_reader("Square.grid");

  Array2D<realT> x = grid_reader.ReadNodes();
  
  //---> Check mesh sizes
  EXPECT_EQ(25, grid_reader.ReadNnode());
  EXPECT_EQ(32, grid_reader.ReadNelement());
  EXPECT_EQ(16, grid_reader.ReadNbcFace());
  EXPECT_EQ(4, grid_reader.ReadNbcID());

  //---> Check nodal coordinates
  EXPECT_DOUBLE_EQ(-1.0, x(0,0));
  EXPECT_DOUBLE_EQ(-0.5, x(0,1));
  EXPECT_DOUBLE_EQ(-0.5, x(5,0));
  EXPECT_DOUBLE_EQ(-0.5, x(5,1));
  EXPECT_DOUBLE_EQ(1.0, x(23,0));
  EXPECT_DOUBLE_EQ(.25, x(23,1));
  
  //---> Check element2node
  List2D<intT> e2n(grid_reader.ReadElement2Node());
  for(intT i = 0; i < e2n.get_lead_size(); i++){
    EXPECT_EQ(3, e2n.get_ncol(i));
  }
 
  EXPECT_EQ(5, e2n(0,1));
  EXPECT_EQ(6, e2n(8,2));
  EXPECT_EQ(6, e2n(9,2));
  
  //---> Check Element Type
  Array1D<ElementTopology::element_types> etype = grid_reader.ReadElementType();
  for(intT i = 0; i < etype.get_size(0); i++){
    EXPECT_EQ(ElementTopology::element_types::TRI, etype(i));
  }

  Array1D<intT> region = grid_reader.ReadElementRegion();
  for(intT i = 0; i < region.get_size(0); i++){
    EXPECT_EQ(0, region(i));
  }

  List2D<intT> bc_face2node = grid_reader.ReadBcFace2Node();
  for(intT i = 0; i < bc_face2node.get_lead_size(); i++){
    EXPECT_EQ(2, bc_face2node.get_ncol(i));
  }
  EXPECT_EQ(0,bc_face2node(0,0));
  EXPECT_EQ(5,bc_face2node(0,1));
  EXPECT_EQ(2,bc_face2node(14,0));
  EXPECT_EQ(1,bc_face2node(14,1));

  Array1D<intT> bc_id = grid_reader.ReadBcID();
  EXPECT_EQ(3, bc_id(15));
  EXPECT_EQ(0, bc_id(0));
  EXPECT_EQ(0, bc_id(2));

  Array1D<ElementTopology::face_types> bc_face_type = grid_reader.ReadBcFaceType();
  for(intT i = 0; i < bc_face_type.get_size(0); i++){
    EXPECT_EQ(ElementTopology::face_types::FACE_BAR, bc_face_type(i));
  }

}
