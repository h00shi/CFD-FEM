#include "Mesh/UnstMeshGeom.h"
//****************************************************************************80
UnstMeshGeom::UnstMeshGeom(UnstMeshReader& mesh_reader)
{
  //---> Read nodes from mesh reader
  x_ = mesh_reader.ReadNodes();
  //---> Query coordinate array for number of nodes 
  nnode_ = x_.get_size(0);
  //---> Query coordinate array for number of dimensions
  ndim_  = x_.get_size(1);
}// End UnstMeshGeom::UnstMeshGeom
//****************************************************************************80
UnstMeshGeom::~UnstMeshGeom(){}
