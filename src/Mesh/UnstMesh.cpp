#include "UnstMesh.h"

//****************************************************************************80
UnstMesh::UnstMesh(UnstMeshReader& mesh_reader) :
mesh_geometry_(mesh_reader),
mesh_elements_(mesh_reader),
mesh_bc_faces_(mesh_reader)
{
  grid_mem_ = 0.0;
}

//****************************************************************************80
UnstMesh::~UnstMesh() {}

//****************************************************************************80
void UnstMesh::Diagnostic(std::ostream& out_stream) 
{
  out_stream << std::endl;
  out_stream << "------------------------- Mesh Diagnostics --------------------"
	    << std::endl;

  out_stream << "Total Grid Memory: " << grid_mem_ << " MB" << std::endl;
  out_stream << "---------------------- End  Mesh Diagnostics -----------------"
	    << std::endl << std::endl;
}


