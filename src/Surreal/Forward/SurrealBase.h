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
template<class DerivedType, class realT, int N>
class SurrealBase
{
public:
  typedef realT realT_;
  static const int N_=N;

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
//****************************************************************************80
//!
//! \brief Value : Get the reference to the value...modifiable
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
  inline realT& Value() {
    //---> Return value to user
    return(value_);
  }// End Value

//****************************************************************************80
//!
//! \brief Value : Get the value const correct version to return to const
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
  inline realT Value() const {
    //---> Return value to user
    return(value_);
  }// End Value
protected:
  realT value_;   /*!< The value of the number. */
}; // End class SurrealBase
#endif
