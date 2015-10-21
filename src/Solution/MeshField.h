/*
 * MeshField.h
 *
 *  Created on: Oct 16, 2015
 *      Author: rabbit
 */

#ifndef MESHFIELD_H_
#define MESHFIELD_H_
#include "my_incl.h"
#include "DataStructures/Array1D.h"
//****************************************************************************80
//! \brief An abstract class of a field that exists on a mesh.
//! \nick
//! \version $Rev$
//! \date $Date$
//****************************************************************************80
class MeshField
{
public:
//****************************************************************************80
//! \brief Constructor for Mesh Field, requires that one know the size of
//!        of the field
//! \nick
//! \version $Rev$
//! \date $Date$
//****************************************************************************80
  MeshField(const intT field_size);
//****************************************************************************80
//! \brief Destructor for the mesh field...pure virtual.
//! \nick
//! \version $Rev$
//! \date $Date$
//****************************************************************************80
  virtual ~MeshField()=0;

  inline Array1D<realT>& get_Data() {return data_;}
  inline const Array1D<realT>& get_Data()const{return data_;}

protected:
  Array1D<realT> data_;//!< Array1D representing all data in the field over
  //****************************************************************************80
  //! \brief Method for initializing the class...allows use of default
  //!        constructor in children.
  //! \nick
  //! \version $Rev$
  //! \date $Date$
  //****************************************************************************80
  void Initialize(const intT field_size);

  MeshField();
private:

};
#endif /* FIELD_H_ */
