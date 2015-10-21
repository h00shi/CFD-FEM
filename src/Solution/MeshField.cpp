/*
 * MeshField.cpp
 *
 *  Created on: Oct 18, 2015
 *      Author: rabbit
 */
#include "Solution/MeshField.h"
//****************************************************************************80
MeshField::MeshField()
{

}
//****************************************************************************80
MeshField::MeshField(const intT field_size)
{
  Initialize(field_size);
}
//****************************************************************************80
MeshField::~MeshField(){}

//****************************************************************************80
void MeshField::Initialize(const intT field_size)
{
  data_.initialize(field_size);
  return;
}
