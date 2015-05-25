// -*-c++-*-
#ifndef SURREALMULTIPLY_H
#define SURREALMULTIPLY_H

#include "SurrealBase.h"
//----------------------------- Surreal * Surreal ------------------------------
//****************************************************************************80
//! \brief SurrealMultiply : A class template to represent automatic
//!                     differentiation of the multiplication operator
//!                     on two arbitrary (Surreal) types
//! \details
//! \nick
//! \tparam LHSType Type used for the arguments on the left hand side of *
//! \tparam RHSType Type used for the arguments on the right hand side of *
//****************************************************************************80
template<class LHSType, class RHSType, int N>
class SurrealMultiply :
  public SurrealBase<SurrealMultiply<LHSType, RHSType, N>, N>
{
private:
  const LHSType& lhs_; //!< Reference to constant object on left hand side
  const RHSType& rhs_; //!< Reference to constant object on right hand side
public:
  typedef typename LHSType::realT_ realT; //get fundamental real type from left
  typedef realT realT_; //store real type

//****************************************************************************80
//! \brief Constructor for constructing multiply operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to constant object on left side of *
//! \param[in] rhs_in The reference to constant object on right side of *
//****************************************************************************80
  SurrealMultiply(const LHSType& lhs_in, const RHSType& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {}

//****************************************************************************80
//! \brief Value : Returns the value of the multiplication between left and
//!                right hand sides
//! \details Value of (Surreal1 * Surreal2) =
//!                   Surreal1.Value() * Surreal2.Value()
//! \nick
//! \return Multiplication of lhs.Value() and rhs.Value()
//****************************************************************************80
  inline realT Value() const {
    return(lhs_.Value() * rhs_.Value());
  }

//****************************************************************************80
//! \brief Deriv : Returns the derivative of the multiplication between left and
//!                right hand sides
//! \details Derivative of (Surreal1 * Surreal2) =
//!                         Surreal1.Deriv() * Surreal2.Value() +
//!                         Surreal1.Value() * Surreal2.Deriv()
//! \nick
//! \param[in] i The index of the derivative you wish to compute
//! \return ith Derivative of (Surreal1 * Surreal2)
//****************************************************************************80
  inline realT Deriv(const int& i) const {
    return(lhs_.Deriv(i)*rhs_.Value() + lhs_.Value()*rhs_.Deriv(i));
  }
}; // End class SurrealMultiply

//----------------------------- Real * Surreal -------------------------------
//****************************************************************************80
//! \brief SurrealMultiply : A class template to represent automatic
//!                     differentiation of the multiplication operator
//!                     on a real type and an arbitrary (Surreal) type
//! \details
//! \nick
//! \tparam RHSType Type used for the argument on the right hand side of *
//****************************************************************************80
template<class RHSType, int N>
class SurrealMultiply<typename RHSType::realT_, RHSType, N> :
  public SurrealBase<SurrealMultiply<typename RHSType::realT_, RHSType, N>, N>
{
public:
  typedef typename RHSType::realT_ realT;  //get fundamental real type from rght
  typedef realT realT_;  //store real type

//****************************************************************************80
//! \brief Constructor for constructing multiply operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to constant real   on left side of *
//! \param[in] rhs_in The reference to constant object on right side of *
//****************************************************************************80
  SurrealMultiply(const typename RHSType::realT_& lhs_in, const RHSType& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {}

//****************************************************************************80
//! \brief Value : Returns the value of the multiplication between left and
//!                right hand sides
//! \details Value of (real * Surreal) = real * Surreal.Value()
//! \nick
//! \return Multiplication of lhs and rhs.Value()
//****************************************************************************80
  inline realT Value() const {
    return(lhs_*rhs_.Value());
  }

//****************************************************************************80
//! \brief Deriv : Returns the derivative of the multiplication between left and
//!                right hand sides
//! \details Derivative of (real * Surreal) = real * Surreal2.Deriv()
//! \nick
//! \param[in] i The index of the derivative you wish to compute
//! \return ith Derivative of (real * Surreal)
//****************************************************************************80
  inline realT Deriv(const int& i) const {
    return(lhs_*rhs_.Deriv(i));
  }
private:
  const typename RHSType::realT_ & lhs_; //left  hand side
  const RHSType& rhs_;                   //right hand side
}; // End class SurrealMultiply

//-------------------------- MULTIPLICATION OPERATORS --------------------------
//****************************************************************************80
//! \brief Operator * : This operator declaration declares how to multiply two
//!                     SurrealBase types
//! \details  SurrealMultiply = SurrealBase * SurrealBase
//! \nick
//! \param[in] lhs The object on left side of * sign
//! \param[in] rhs The object on right side of * sign
//! \return SurrealMultiply object to represent multiplication of lhs and rhs
//****************************************************************************80
template<class LHSType, class RHSType, int N>
inline SurrealMultiply<LHSType, RHSType, N>
operator*(const SurrealBase<LHSType, N>& lhs,
          const SurrealBase<RHSType, N>& rhs);

//****************************************************************************80
//! \brief Operator * : This operator declaration declares how to multiply a
//!                     real number with a SurrealBase type
//! \details SurrealMultiply = real * SurrealBase
//! \nick
//! \param[in] lhs The real number on left side of * sign
//! \param[in] rhs The object on right side of * sign
//! \return SurrealMultiply object to represent multiplication of lhs and rhs
//****************************************************************************80
template<class RHSType, int N>
inline SurrealMultiply<typename RHSType::realT_, RHSType, N>
operator*
(const typename RHSType::realT_& lhs, const SurrealBase<RHSType, N>& rhs);

//****************************************************************************80
//! \brief Operator * : This operator declaration declares how to multiply a
//!                     SurrealBase type with a real number
//! \details  SurrealMultiply = SurrealBase * real
//!           Since multiplication is commutative, this function will simply
//!           flip the order of the arguments and use the SurrealMultiply class
//!           specialization with the real type in the left hand side.
//! \nick
//! \param[in] lhs The object on left side of * sign
//! \param[in] rhs The real number on right side of * sign
//! \return SurrealMultiply object to represent multiplication of lhs and rhs
//****************************************************************************80
template<class LHSType, int N>
inline SurrealMultiply<typename LHSType::realT_, LHSType, N>
operator*
(const SurrealBase<LHSType, N>& lhs, const typename LHSType::realT_& rhs);

#include "SurrealMultiply_Imple.h"
#endif
