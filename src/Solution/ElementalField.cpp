/*
 * ElementalField.cpp
 *
 *  Created on: Oct 21, 2015
 *      Author: rabbit
 */

#include "Solution/ElementalField.h"

//****************************************************************************80
ElementalField::ElementalField(const UnstMesh& mesh,
                               const intT& ndof,
                               const intT& nvar,
                               const List2D<intT>& elem2dof) : mesh_(mesh)

{
  nvar_.initialize(ndof);
  for(intT i = 0; i < nvar_.get_size(0); i++){
    nvar_(i) = nvar;
  }

  Field::data_.initialize(nvar_);

  Field::elem2dof_.initialize_copy_pattern(elem2dof);
  for(intT i = 0; i < elem2dof.get_total_size(); i++){
    elem2dof_(i) = elem2dof(i);
  }

  Array1D<intT> count_type(7);
  for(intT e = 0; e < mesh_.get_MeshElements().get_nelement(); e++)
  {
    count_type(mesh_.get_MeshElements().get_element_type()(e))++;
  }

  Field::elem_groups_.initialize(count_type);

  count_type.set_value(0);
  for(intT e = 0; e < mesh_.get_MeshElements().get_nelement(); e++)
  {
    ElementTopology::element_types etype =
        mesh_.get_MeshElements().get_element_type()(e);

    intT i = count_type(etype);
    Field::elem_groups_(etype,i) = e;

    count_type(etype)++;
  }

}
//****************************************************************************80
ElementalField::~ElementalField(){}
