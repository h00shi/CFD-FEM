#include "gtest/gtest.h"
#include "Mesh/CGMesh.h"
#include "IO/UnstMeshWriterVTK.h"
#include "IO/UnstMeshReaderGMSH.h"
#include "ParallelComm/Communication.h"
TEST(VTK, Write){
  Communication::Initialize();
  UnstMeshReaderNKBGrid grid_reader("Square.grid");
  CGMesh mesh(grid_reader);
  mesh.Diagnostic(std::cout);
  UnstMeshWriterVTK writer(mesh);
  
  writer.Write("Test");

}

TEST(VTK, Write_2D_Mixed){
  UnstMeshReaderGMSH grid_reader("Mixed_2D.msh-bin", "Mixed_2D.idmap");
  CGMesh mesh(grid_reader);
  mesh.Diagnostic(std::cout);
  UnstMeshWriterVTK writer(mesh);

  writer.Write("Test_Mixed_2D");

}

TEST(VTK, Write_3D){
  UnstMeshReaderGMSH grid_reader("Mixed_3D.msh-bin", "Mixed_3D.idmap");
  CGMesh mesh(grid_reader);
  mesh.Diagnostic(std::cout);
  UnstMeshWriterVTK writer(mesh);

  writer.Write("Test_Mixed_3D");
  Communication::Finalize();
}
