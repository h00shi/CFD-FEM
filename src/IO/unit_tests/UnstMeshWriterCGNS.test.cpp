#include "gtest/gtest.h"
#include "Mesh/UnstMesh.h"
#include "IO/UnstMeshWriterCGNS.h"

TEST(CGNS, Write){
  UnstMesh mesh("Square.grid","GRID-NKB");
  UnstMeshWriterCGNS writer(mesh);
  
  writer.Write("Test");

}
