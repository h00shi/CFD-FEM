#include "Mesh/DummyMesh.h"
#include "gtest/gtest.h"
#include "Include/my_incl.h"

static const realT small = 5.0e-15;
TEST(DummyMesh,CubicHexMesh) 
{
  std::string mesh = DummyMesh::SetupCubicHexMesh(3.0,3.0,3.0,4,4,4);
  std::cout << mesh << std::endl;
}
