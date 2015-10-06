#include "Mesh/UnstMeshGeom.h"
#include "DataStructures/Array2D.h"
#include "gtest/gtest.h"

TEST(UnstMeshGoem, Init) 
{
  const intT N = 11;
  Array2D<realT> x(N,3);

  for(intT i = 0; i < N; i++){
    x(i,0) = (realT)i + 0;
    x(i,1) = (realT)i + 1;
    x(i,2) = (realT)i + 2;
  }
  
  // UnstMeshGeom mesh_geom(std::move(x));

  // EXPECT_EQ(0,x.get_size(0));
  // EXPECT_EQ(0,x.get_size(1));
  
  // for(intT i = 0; i < N; i++){
  //   EXPECT_DOUBLE_EQ((realT)i + 0, mesh_geom.get_x()(i,0));
  //   EXPECT_DOUBLE_EQ((realT)i + 1, mesh_geom.get_x()(i,1));
  //   EXPECT_DOUBLE_EQ((realT)i + 2, mesh_geom.get_x()(i,2));
  // }
  
}
