//-*-c++-*-
#ifndef SURREALRBASE_H
#define SURREALRBASE_H

template<class DerivedType, class realT, int N>
class SurrealRBase
{
public:
  typedef realT realT_;
  static const int N_=N;

  SurrealRBase() : value_(0.0) {}

  inline SurrealRBase& operator =
  (SurrealRBase<DerivedType,realT,N> const &) = default;

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
  DerivedType const & CastToDerived () const
  {
    return static_cast<DerivedType const &> (*this);
  } // End CastToDerived

//****************************************************************************80
//! \brief Value : Get the reference to the value...modifiable
//****************************************************************************80
  inline realT& Value() {
    //---> Return value to user
    return(value_);
  }// End Value

//****************************************************************************80
//!
//! \brief Value : Get the value const correct version to return to const
//****************************************************************************80
  inline realT Value() const {
    //---> Return value to user
    return(value_);
  }// End Value

//****************************************************************************80
//! \brief GetValue : Returns value
//! \nick
//****************************************************************************80
  inline realT GetValue() const {
    return(value_);
  }
//****************************************************************************80
//! \brief SetValue : Sets the value
//! \nick
//! \param[in] val - value to set to
//****************************************************************************80
  void SetValue(const realT val) {
    value_ = val;
  }// End SetValue
protected:
  realT value_;   /*!< The value of the number. */
}; // End class SurrealRBase
#endif
