#include "gtest/gtest.h"
#include "Mesh/CGMesh.h"
#include "IO/UnstMeshWriterVTK.h"

TEST(VTK, Write){
  UnstMeshReaderNKBGrid grid_reader("Square.grid");
  CGMesh mesh(grid_reader);
  mesh.Diagnostic(std::cout);
  UnstMeshWriterVTK writer(mesh);
  
  writer.Write("Test");

}
