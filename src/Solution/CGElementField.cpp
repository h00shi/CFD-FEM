/*
 * CGElementField.cpp
 *
 *  Created on: Oct 20, 2015
 *      Author: rabbit
 */

#include "Solution/CGElementField.h"

//****************************************************************************80
CGElementField::CGElementField(const CGMesh& mesh, const intT& nvar) :
ElementalField(mesh, mesh.get_MeshGeom().get_nnode(), nvar,
               mesh.get_MeshElements().get_element2node())
{


}
//****************************************************************************80
CGElementField::~CGElementField(){}
