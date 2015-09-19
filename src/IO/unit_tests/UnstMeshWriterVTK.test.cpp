#include "gtest/gtest.h"
#include "Mesh/UnstMesh.h"
#include "IO/UnstMeshWriterVTK.h"

TEST(CGNS, Write){
  UnstMesh mesh("Square.grid","Grid-NKB");
  mesh.Diagnostic(std::cout);
  UnstMeshWriterVTK writer(mesh);
  
  writer.Write("Test");

}
