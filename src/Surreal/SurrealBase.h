//-*-c++-*-
#ifndef SURREALBASE_H
#define SURREALBASE_H
#include <cmath>
#include <iostream>
//****************************************************************************80
//! \brief This is the base class template for all Surreal data types and
//!         surreal operation types.
//! \details  This class template is the key component of the curiously
//!  recurring template pattern that we will use to formulate automatic
//!  differentiation operator classes and operators to construct these classes.
//! \nick
//! \version $Rev$
//! \tparam DerivedType A type that indicates the derived type for which we are
//!                       creating a base class.
//! \tparam N           The number of derivatives you want to take
//****************************************************************************80
template<class DerivedType, int N>
class SurrealBase
{
  //---> There are no private data members

public:
//****************************************************************************80
//! \brief CastToDerived : A method to cast the current instance of
//!                        SurrealBase to a reference to a constant
//!                        instance of DerivedType type.
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//!
//****************************************************************************80
  const DerivedType& CastToDerived () const
  {
    return static_cast<DerivedType const &> (*this);
  } // End CastToDerived
}; // End class SurrealBase
#endif
