// -*-c++-*-
#ifndef SURREALDIVIDE_H
#define SURREALDIVIDE_H

#include "Surreal/SurrealBase.h"
//----------------------------- Surreal / Surreal ------------------------------
//****************************************************************************80
//! \brief SurrealDivide : A class template to represent automatic
//!                     differentiation of the division operator
//!                     on two arbitrary (Surreal) types
//! \details
//! \nick
//! \tparam LHSType Type used for the arguments on the left hand side of /
//! \tparam RHSType Type used for the arguments on the right hand side of /
//****************************************************************************80
template<class LHSType, class RHSType, int N>
class SurrealDivide : public SurrealBase< SurrealDivide<LHSType, RHSType, N>, N >
{
public:
  typedef typename LHSType::realT_ realT; //get fundamental real type from left
  typedef realT realT_; //store real type

//****************************************************************************80
//! \brief Constructor for constructing divide operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to constant object on left side of /
//! \param[in] rhs_in The reference to constant object on right side of /
//****************************************************************************80
  SurrealDivide(const LHSType& lhs_in, const RHSType& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {}

//****************************************************************************80
//! \brief Value : Returns the value of the left hand side divided by the
//!                right hand side
//! \details Value of (Surreal1 / Surreal2) =
//!                   Surreal1.Value() / Surreal2.Value()
//! \nick
//! \return Division of lhs.Value() by rhs.Value()
//****************************************************************************80
  inline realT Value() const {
    return(lhs_.Value() / rhs_.Value());
  }

//****************************************************************************80
//! \brief Deriv : Returns the derivative of the left hand side divided by the
//!                right hand side
//! \details Derivative of (Surreal1 / Surreal2) =
//! (Surreal1.Deriv() * Surreal2.Value() - Surreal1.Value() * Surreal2.Deriv())/
//! (Surreal2.Value)^2
//! \nick
//! \param[in] i The index you of the Derivative you wish to compute
//! \return ith Derivative of (Surreal1 / Surreal2)
//****************************************************************************80
  inline realT Deriv(const int& i) const {
    return((lhs_.Deriv(i)*rhs_.Value() - lhs_.Value()*rhs_.Deriv(i))
           / (rhs_.Value() * rhs_.Value()));
  }
private:
  const LHSType& lhs_; //!< Reference to constant object on left hand side
  const RHSType& rhs_; //!< Reference to constant object on right hand side
}; // End class SurrealDivide

//----------------------------- Real / Surreal -------------------------------
//****************************************************************************80
//! \brief SurrealDivide : A class template to represent automatic
//!                     differentiation of the division of a real type by
//!                     an arbitrary (Surreal) type
//! \details
//! \nick
//! \tparam RHSType Type used for the argument on the right hand side of /
//****************************************************************************80
template<class RHSType, int N>
class SurrealDivide<typename RHSType::realT_, RHSType, N> :
  public SurrealBase<SurrealDivide<typename RHSType::realT_, RHSType, N>, N>
{
public:
  typedef typename RHSType::realT_ realT; //get fundamental real type from right
  typedef realT realT_; //store real type

//****************************************************************************80
//! \brief Constructor for constructing divide operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to constant real   on left side of /
//! \param[in] rhs_in The reference to constant object on right side of /
//****************************************************************************80
  SurrealDivide(const typename RHSType::realT_& lhs_in, const RHSType& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {}

//****************************************************************************80
//! \brief Value : Returns the value of the left hand side divided by the
//!                right hand side
//! \details Value of (real / Surreal) = real / Surreal.Value()
//! \nick
//! \return Division of lhs by rhs.Value()
//****************************************************************************80
  inline realT Value() const {
    return(lhs_/rhs_.Value());
  }

//****************************************************************************80
//! \brief Deriv : Returns the derivative of the left hand side divided by the
//!                right hand side
//! \details Derivative of (real / Surreal) =
//!          - real / Surreal.Value()^2 * Surreal.Deriv()
//! \nick
//! \param[in] i The index you of the Derivative you wish to compute
//! \return ith Derivative of (real / Surreal)
//****************************************************************************80
  inline realT Deriv(const int& i) const {
    return(-lhs_ / (rhs_.Value() * rhs_.Value()) * rhs_.Deriv(i));
  }
private:
  const typename RHSType::realT_ & lhs_; //left  hand side
  const RHSType& rhs_;                   //right hand side
}; // End class SurrealDivide

//----------------------------- Surreal / real -------------------------------
//****************************************************************************80
//! \brief SurrealDivide : A class template to represent automatic
//!                     differentiation of the division of an arbitrary
//!                     (Surreal) type by a real type
//! \details
//! \nick
//! \tparam LHSType Type used for the argument on the left hand side of /
//****************************************************************************80
template<class LHSType, int N>
class SurrealDivide<LHSType, typename LHSType::realT_, N> :
  public SurrealBase<SurrealDivide<LHSType, typename LHSType::realT_,N>, N>
{
public:
  typedef typename LHSType::realT_ realT; //get fundamental real type from left
  typedef realT realT_; //store real type

//****************************************************************************80
//! \brief Constructor for constructing divide operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to constant object on left side of /
//! \param[in] rhs_in The reference to constant real   on right side of /
//****************************************************************************80
  SurrealDivide(const LHSType& lhs_in, const typename LHSType::realT_& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {}

//****************************************************************************80
//! \brief Value : Returns the value of the left hand side divided by the
//!                right hand side
//! \details Value of (Surreal / real) = Surreal.Value() / real
//! \nick
//! \return Division of lhs.Value() by rhs
//****************************************************************************80
  inline realT Value() const {
    return(lhs_.Value() /rhs_);
  }

//****************************************************************************80
//! \brief Deriv : Returns the derivative of the left hand side divided by the
//!                right hand side
//! \details Derivative of (Surreal / real) = Surreal.Deriv() / real
//! \nick
//! \param[in] i The index you of the derivative you wish to compute
//! \return ith Derivative of (Surreal / real)
//****************************************************************************80
  inline realT Deriv(const int& i) const {
    return( lhs_.Deriv(i) / rhs_);
  }
private:
  const LHSType& lhs_;                   //left  hand side
  const typename LHSType::realT_ & rhs_; //right hand side
}; // End class SurrealDivide

//-------------------------- DIVISION OPERATORS --------------------------------
//****************************************************************************80
//! \brief Operator / : This operator declaration declares how to divide two
//!                     SurrealBase types
//! \details  SurrealDivide = SurrealBase / SurrealBase
//! \nick
//! \param[in] lhs The object on left side of / sign
//! \param[in] rhs The object on right side of / sign
//! \return SurrealDivide object to represent division of lhs by rhs
//****************************************************************************80
template<class LHSType, class RHSType, int N>
inline SurrealDivide<LHSType, RHSType, N>
operator /(const SurrealBase< LHSType, N>& lhs,
           const SurrealBase< RHSType, N>& rhs);

//****************************************************************************80
//! \brief Operator / : This operator declaration declares how to divide a
//!                     real number with a SurrealBase type
//! \details SurrealDivide = real / SurrealBase
//! \nick
//! \param[in] lhs The real number on left side of / sign
//! \param[in] rhs The object on right side of / sign
//! \return SurrealDivide object to represent division of lhs by rhs
//****************************************************************************80
template<class RHSType, int N>
inline SurrealDivide<typename RHSType::realT_, RHSType, N>
operator/
(const typename RHSType::realT_& lhs, const SurrealBase< RHSType,N >& rhs);

//****************************************************************************80
//! \brief Operator / : This operator declaration declares how to divide a
//!                     SurrealBase type with a real number
//! \details SurrealDivide = SurrealBase / real
//! \nick
//! \param[in] lhs The object on left side of / sign
//! \param[in] rhs The real number on right side of / sign
//! \return SurrealDivide object to represent division of lhs by rhs
//****************************************************************************80
template<class LHSType, int N>
inline SurrealDivide<LHSType, typename LHSType::realT_, N>
operator/
(const SurrealBase<LHSType, N>& lhs, const typename LHSType::realT_& rhs);

#include "SurrealDivide_Imple.h"
#endif
