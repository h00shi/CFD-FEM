//-*-c++-*-
#ifndef SURREALADD_H
#define SURREALADD_H

#include "Surreal/SurrealBase.h"

//----------------------------- Surreal + Surreal ------------------------------
//****************************************************************************80
//! \brief SurrealAdd : A class template to represent automatic
//!                     differentiation of the addition operator
//!                     on two arbitrary (Surreal) types
//! \details
//! \nick
//! \tparam LHSType Type used for the argument on the left hand side of +
//! \tparam RHSType Type used for the argument on the right hand side of +
//****************************************************************************80
template<class LHSType, class RHSType, int N>
class SurrealAdd : public SurrealBase<SurrealAdd<LHSType, RHSType, N>, N>
{
public:
  typedef typename LHSType::realT_ realT; //get fundamental real type from left
  typedef realT realT_; //store real type

//****************************************************************************80
//! \brief Constructor for constructing add operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to constant object on left side of +
//! \param[in] rhs_in The reference to constant object on right side of +
//****************************************************************************80
  SurrealAdd(const LHSType& lhs_in, const RHSType& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {}

//****************************************************************************80
//! \brief Value : Returns the value of the addition between left and right hand
//!                sides
//! \details Value of (Surreal1 + Surreal2) =
//!                    Surreal1.Value() + Surreal2.Value()
//! \nick
//! \return Addition of lhs.Value() and rhs.Value()
//****************************************************************************80
  inline realT Value() const {
    return(lhs_.Value() + rhs_.Value());
  }

//****************************************************************************80
//! \brief Deriv : Returns the Derivative of lhs + rhs
//! \details Derivative of (Surreal1 + Surreal2) =
//!                         Surreal1.Deriv() + Surreal2.Deriv()
//! \nick
//! \param[in] i The index of the derivative you wish to compute
//! \return Addition of lhs.Deriv(i) and rhs.Deriv(i)
//****************************************************************************80
  inline realT Deriv(const int& i) const {
    return(lhs_.Deriv(i) + rhs_.Deriv(i));
  }
private:
  const LHSType& lhs_;//!< Reference to constant object on left hand side
  const RHSType& rhs_;//!< Reference to constant object on right hand side
}; // End class SurrealAdd

//----------------------------- Real + Surreal ---------------------------------
//****************************************************************************80
//! \brief SurrealAdd : A class template to represent automatic
//!                     differentiation of the addition operator
//!                     on a real type + arbitrary (Surreal) type
//! \details
//! \nick
//! \tparam RHSType Type used for the argument on the right hand side of +
//****************************************************************************80
template<class RHSType, int N>
class SurrealAdd<typename RHSType::realT_, RHSType, N> :
  public SurrealBase<SurrealAdd<typename RHSType::realT_, RHSType, N>, N>
{
public:
  typedef typename RHSType::realT_ realT; //get fundamental real type from right
  typedef realT realT_; //store real type

//****************************************************************************80
//! \brief Constructor for constructing add operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to constant real   on left side of +
//! \param[in] rhs_in The reference to constant object on right side of +
//****************************************************************************80
  SurrealAdd(const typename RHSType::realT_& lhs_in, const RHSType& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {}

//****************************************************************************80
//! \brief Value : Returns the value of the addition between left and right hand
//!                sides
//! \details  Value of (real + Surreal) = real + Surreal.Value()
//! \nick
//! \return Addition of lhs and rhs.Value()
//****************************************************************************80
  inline realT Value() const {
    return(lhs_ + rhs_.Value());
  }

//****************************************************************************80
//! \brief Deriv : Returns the Derivative of lhs + rhs
//! \details  Derivative of (real + Surreal) = Surreal.Deriv()
//! \nick
//! \param[in] i The index of the derivative you wish to compute
//! \return rhs.Deriv(i)
//****************************************************************************80
  inline realT Deriv(const int& i) const {
    return(rhs_.Deriv(i));
  }
private:
  const typename RHSType::realT_ & lhs_; //left  hand side
  const RHSType& rhs_;                   //right hand side
}; // End class SurrealAdd

//----------------------------- ADDITION OPERATORS -----------------------------
//****************************************************************************80
//! \brief Operator + : This operator declaration declares how to add two
//!                     SurrealBase types
//! \details  SurrealAdd = SurrealBase + SurrealBase
//! \nick
//! \param[in] lhs The object on left side of + sign
//! \param[in] rhs The object on right side of + sign
//! \return SurrealAdd object to represent addition of lhs + rhs
//****************************************************************************80
template<class LHSType, class RHSType, int N>
inline SurrealAdd<LHSType, RHSType, N>
operator+(const SurrealBase< LHSType, N>& lhs,
          const SurrealBase< RHSType, N>& rhs);

//****************************************************************************80
//! \brief Operator + : This operator declaration declares how to add a real
//!                     number with a SurrealBase type
//! \details SurrealAdd = real + SurrealBase
//! \nick
//! \param[in] lhs The real number on left side of + sign
//! \param[in] rhs The object on right side of + sign
//! \return SurrealAdd object to represent addition of lhs + rhs
//****************************************************************************80
template<class RHSType, int N>
inline SurrealAdd<typename RHSType::realT_, RHSType, N>
operator+
(const typename RHSType::realT_& lhs, const SurrealBase<RHSType, N>& rhs);

//****************************************************************************80
//! \brief Operator + : This operator declaration declares how to add a
//!                     SurrealBase type with a real number
//! \details  SurrealAdd = SurrealBase + real
//!           Since addition is commutative, this function will simply flip the
//!           order of the arguments and use the SurrealAdd class specialization
//!           with the real type in the left hand side.
//! \nick
//! \param[in] lhs The object on left side of + sign
//! \param[in] rhs The real number on right side of + sign
//! \return SurrealAdd object to represent addition of lhs + rhs
//****************************************************************************80
template<class LHSType, int N>
inline SurrealAdd<typename LHSType::realT_, LHSType, N>
operator+
(const SurrealBase<LHSType, N>& lhs, const typename LHSType::realT_& rhs);

#include "SurrealAdd_Imple.h"
#endif
