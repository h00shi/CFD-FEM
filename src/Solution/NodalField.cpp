/*
 * NodalField.cpp
 *
 *  Created on: Oct 20, 2015
 *      Author: rabbit
 */

#include "Solution/NodalField.h"

//****************************************************************************80
NodalField::NodalField(const UnstMesh& mesh, const Array1D<intT>& nvar) :
mesh_(mesh)
{

  data_index_.initialize(mesh.get_MeshGeom().get_nnode() + 1);
  nvar_.initialize(mesh.get_MeshGeom().get_nnode());

  intT field_size = 0;
  data_index_(0) = 0;
  for(intT i = 0; i < nvar.get_size(0); i++){
    field_size += nvar(i);
    nvar_(i) = nvar(i);
    data_index_(i + 1) = data_index_(i) + nvar_(i);
  }

  MeshField::Initialize(field_size);
}
//****************************************************************************80
NodalField::~NodalField(){}
