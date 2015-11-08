/*
 * Field.h
 *
 *  Created on: Oct 16, 2015
 *      Author: rabbit
 */

#ifndef FIELD_H_
#define FIELD_H_
#include "my_incl.h"
#include "DataStructures/List2D.h"

//****************************************************************************80
//! \brief An abstract class of a field that exists on a mesh.
//! \nick
//! \version $Rev$
//! \date $Date$
//****************************************************************************80
class Field
{
public:
//****************************************************************************80
//! \brief Destructor for the mesh field...pure virtual.
//! \nick
//! \version $Rev$
//! \date $Date$
//****************************************************************************80
  virtual ~Field()=0;

  inline intT get_nDof(){return data_.get_lead_size();}
  inline List2D<realT>& get_Data() {return data_;}
  inline const List2D<realT>& get_Data()const {return data_;}
  inline const List2D<intT>& get_Elem2Dof(){return elem2dof_;}
  inline const intT get_nVar(const intT& dof){return data_.get_ncol(dof);}
//****************************************************************************80
//! \brief Operator for getting data for node n and variable j
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[in] idof Degree of freedom for which you want data
//! \param[in] ivar The index for the variable you want
//****************************************************************************80
  inline realT& operator()(const intT& idof, const intT& ivar)
  {
    return data_(idof, ivar);
  }

//****************************************************************************80
//! \brief Operator for getting data for node n and variable j const version
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[in] idof Degree of freedom for which you want data
//! \param[in] ivar The index for the variable you want
//****************************************************************************80
  inline const realT& operator()(const intT& idof, const intT& ivar) const
  {
    return data_(idof, ivar);
  }
//****************************************************************************80
//! \brief Provides the beginning pointer to node data
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[in] idof The degree of freedom at which you want the pointer
//! \return Pointer to data beginning at node n
//****************************************************************************80
  inline realT* DataPtr(const intT& n)
  {
    return Field::data_.get_ptr(n,0);
  }
//****************************************************************************80
//! \brief Provides the beginning pointer to node data const version.
//! \nick
//! \version $Rev$
//! \date $Date$
//! \return Pointer to data beginning at node n
//****************************************************************************80
  inline const realT* DataPtr(const intT& n) const
  {
    return Field::data_.get_ptr(n,0);
  }

protected:
  List2D<realT> data_;//!< Array1D representing all data in the field over
  List2D<intT> elem2dof_;//!< A list that maps each element's degrees of freedom
  List2D<intT> elem_groups_;//!< Lists elements by region number
  Field();

private:

};
#endif /* FIELD_H_ */
