/*
 * Element.h
 *
 *  Created on: Nov 8, 2015
 *      Author: rabbit
 */

#ifndef ELEMENT_H_
#define ELEMENT_H_

//****************************************************************************80
#include "my_incl.h"
#include "consts.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/Array3D.h"
#include "DataStructures/Array4D.h"
class Element {
//****************************************************************************80
//!
//! \brief Element: The constructor, it's empty ...but we declare anyway
//! \details
//! \nick
//! \version
//! \date
//!
//****************************************************************************80
  Element(){}

//****************************************************************************80
//!
//! \brief ~Element : The destructor but we declare it anyway
//! \details
//! \nick
//! \version
//! \date
//!
//****************************************************************************80
  virtual ~Element(){}

private:

  Array2D<realT> dof_; //!< Tabulation of degrees of freedom
  Array2D<realT> basis_; /*!< Tabulation of basis functions at specified
                                physical points*/
  Array3D<realT> grad_basis_;/*!< Tabulation of the gradient of the the basis
                                  functions.*/
};
#endif /* ELEMENT_H_ */
