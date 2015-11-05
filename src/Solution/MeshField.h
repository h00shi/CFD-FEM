/*
 * MeshField.h
 *
 *  Created on: Oct 16, 2015
 *      Author: rabbit
 */

#ifndef MESHFIELD_H_
#define MESHFIELD_H_
#include "my_incl.h"
#include "DataStructures/List2D.h"
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

  inline void set_Data(List2D<realT>* data){data_ = data;};
  inline List2D<realT>& get_Data() {return *data_;}
  inline const List2D<realT>& get_Data()const{return *data_;}

protected:
  List2D<realT>* data_;//!< Array1D representing all data in the field over

  MeshField();
private:

};
#endif /* FIELD_H_ */
