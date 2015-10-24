#ifndef SURREALMULTIPLY_H
#define SURREALMULTIPLY_H

#include "Surreal/Forward/SurrealBase.h"
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
template<class LHSType, class RHSType>
class SurrealMultiply :
  public SurrealBase<SurrealMultiply<LHSType, RHSType>,
                     typename RHSType::realT_, RHSType::N_>
{
private:
  typedef typename RHSType::realT_ realT;  //get fundamental real type from rght
  const LHSType& lhs_; //!< Reference to constant object on left hand side
  const RHSType& rhs_; //!< Reference to constant object on right hand side
public:
//****************************************************************************80
//! \brief Constructor for constructing multiply operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to constant object on left side of *
//! \param[in] rhs_in The reference to constant object on right side of *
//****************************************************************************80
  SurrealMultiply(const LHSType& lhs_in, const RHSType& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {
    static_assert(std::is_same<typename LHSType::realT_,
                  typename RHSType::realT_>::value,
                  "Surreal binary operations require the same floating-point "
                  "data type on left and right sides");
    static_assert(LHSType::N_ == RHSType::N_,
                  "Surreal binary operations require the same number of "
                  "derivatives on left and right sides");
    this->value_ = lhs_.Value() * rhs_.Value();
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
  inline realT Deriv(const int i) const {
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
template<class RHSType>
class SurrealMultiply<typename RHSType::realT_, RHSType> :
  public SurrealBase<SurrealMultiply<typename RHSType::realT_, RHSType>,
                     typename RHSType::realT_, RHSType::N_>
{
private:
  typedef typename RHSType::realT_ realT;  //get fundamental real type from rght
  const realT lhs_; //left  hand side
  const RHSType& rhs_;                   //right hand side
public:
//****************************************************************************80
//! \brief Constructor for constructing multiply operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to constant real   on left side of *
//! \param[in] rhs_in The reference to constant object on right side of *
//****************************************************************************80
  SurrealMultiply(const realT lhs_in, const RHSType& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {
    this->value_ = lhs_*rhs_.Value();
  }
//****************************************************************************80
//! \brief Deriv : Returns the derivative of the multiplication between left and
//!                right hand sides
//! \details Derivative of (real * Surreal) = real * Surreal2.Deriv()
//! \nick
//! \param[in] i The index of the derivative you wish to compute
//! \return ith Derivative of (real * Surreal)
//****************************************************************************80
  inline realT Deriv(const int i) const {
    return(lhs_*rhs_.Deriv(i));
  }
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
template<class LHSType, class RHSType, class realT, int N>
inline SurrealMultiply<LHSType, RHSType>
operator*(const SurrealBase<LHSType, realT, N>& lhs,
          const SurrealBase<RHSType, realT, N>& rhs)
{
  return(SurrealMultiply<LHSType,RHSType>
         (lhs.CastToDerived(), rhs.CastToDerived()));
}

//****************************************************************************80
//! \brief Operator * : This operator declaration declares how to multiply a
//!                     real number with a SurrealBase type
//! \details SurrealMultiply = real * SurrealBase
//! \nick
//! \param[in] lhs The real number on left side of * sign
//! \param[in] rhs The object on right side of * sign
//! \return SurrealMultiply object to represent multiplication of lhs and rhs
//****************************************************************************80
template<class RHSType>
inline SurrealMultiply<typename RHSType::realT_, RHSType>
operator*
(const typename RHSType::realT_ lhs,
 const SurrealBase<RHSType,typename RHSType::realT_, RHSType::N_>& rhs)
{
  return(SurrealMultiply<typename RHSType::realT_,RHSType>
         (lhs, rhs.CastToDerived()));
}

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
template<class LHSType>
inline SurrealMultiply<typename LHSType::realT_, LHSType>
operator*
(const SurrealBase<LHSType, typename LHSType::realT_, LHSType::N_>& lhs,
 const typename LHSType::realT_ rhs)
{
  return(SurrealMultiply<typename LHSType::realT_,LHSType>
         (rhs, lhs.CastToDerived()));
}

#endif
