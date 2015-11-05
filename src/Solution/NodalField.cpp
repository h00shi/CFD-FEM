/*
 * NodalField.cpp
 *
 *  Created on: Oct 20, 2015
 *      Author: rabbit
 */

#include "Solution/NodalField.h"

//****************************************************************************80
NodalField::NodalField(const UnstMesh& mesh, const Array1D<intT>& nvar) :
mesh_(mesh), nvar_(nvar), elem2dof_(mesh.get_MeshElements().get_element2node())
{

}
//****************************************************************************80
NodalField::~NodalField(){}
