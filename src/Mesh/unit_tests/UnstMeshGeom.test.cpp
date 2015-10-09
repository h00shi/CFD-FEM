#include "Mesh/UnstMeshGeom.h"
#include "DataStructures/Array2D.h"
#include "gtest/gtest.h"

TEST(UnstMeshGoem, Init) 
{
  UnstMeshReader* mesh_reader = new UnstMeshReaderNKBGrid("Square.grid");
  UnstMeshGeom mesh_geom(*mesh_reader);
  
  EXPECT_EQ(25, mesh_geom.get_nnode());
  EXPECT_EQ(2, mesh_geom.get_ndim());
  const Array2D<realT>& x = mesh_geom.get_x();
  
  //---> Check nodal coordinates
  EXPECT_DOUBLE_EQ(-1.0, x(0,0));
  EXPECT_DOUBLE_EQ(-0.5, x(0,1));
  EXPECT_DOUBLE_EQ(-0.5, x(5,0));
  EXPECT_DOUBLE_EQ(-0.5, x(5,1));
  EXPECT_DOUBLE_EQ(1.0, x(23,0));
  EXPECT_DOUBLE_EQ(.25, x(23,1));
  
  
}
