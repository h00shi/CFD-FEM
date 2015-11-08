/*
 * CGElementField.h
 *
 *  Created on: Oct 16, 2015
 *      Author: rabbit
 */

#ifndef CGELEMENTFIELD_H_
#define CGELEMENTFIELD_H_
#include "Mesh/ElementTopology.h"
#include "Mesh/CGMesh.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "Solution/ElementalField.h"
//****************************************************************************80
//! \class CGElementField
//! \brief And interpretation of a MeshField as a field over the nodes of the
//!        mesh.  This field is naturally continuous at all element boundaries
//! \nick
//! \version $Rev$
//! \date $Date$
//!
//****************************************************************************80
class CGElementField : public ElementalField
{
public:
//****************************************************************************80
//! \brief Constructor
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[in] mesh The mesh on which we are defining the field
//****************************************************************************80
  CGElementField(const CGMesh& mesh, const intT& nvar);
//****************************************************************************80
//! \brief Destructor
//! \nick
//! \version $Rev$
//! \date $Date$
//!
//****************************************************************************80
  ~CGElementField();

private:
  CGElementField() = delete;
};



#endif /* CGELEMENTFIELD_H_ */
