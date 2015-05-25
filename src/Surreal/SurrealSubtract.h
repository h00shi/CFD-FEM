//-*-c++-*-
#ifndef SURREALSUBTRACT_H
#define SURREALSUBTRACT_H

#include "SurrealBase.h"
//----------------------------- Surreal - Surreal ------------------------------
//****************************************************************************80
//! \brief SurrealSubtract : A class template to represent automatic
//!                     differentiation of the subtraction of an arbitrary type
//!                     from another arbitrary type
//! \details
//! \nick
//! \tparam LHSType Type used for the argument on the left hand side of -
//! \tparam RHSType Type used for the argument on the right hand side of -
//****************************************************************************80
template<class LHSType, class RHSType, int N>
class SurrealSubtract :
  public SurrealBase<SurrealSubtract<LHSType, RHSType, N>, N>
{
public:
  typedef typename LHSType::realT_ realT; //get fundamental real type from left
  typedef realT realT_; //store real type

//****************************************************************************80
//! \brief Constructor for constructing subtraction operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to constant object on left side of -
//! \param[in] rhs_in The reference to constant object on right side of -
//****************************************************************************80
  SurrealSubtract(const LHSType& lhs_in, const RHSType& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {}

//****************************************************************************80
//! \brief Value : Returns the value of the left hand side minus the right hand
//!                side
//! \details Value of (Surreal1 - Surreal2) =
//!                    Surreal1.Value() - Surreal2.Value()
//! \nick
//! \return Value of lhs minus rhs
//****************************************************************************80
  inline realT Value() const {
    return(lhs_.Value() - rhs_.Value());
  }

//****************************************************************************80
//! \brief Deriv : Returns the derivative of lhs - rhs
//! \details Derivative of (Surreal1 - Surreal2) =
//!                         Surreal1.Deriv() - Surreal2.Deriv()
//! \nick
//! \param[in] i The index of the derivative you wish to compute
//! \return ith derivative of lhs minus rhs
//****************************************************************************80
   inline realT Deriv(const int& i) const {
    return(lhs_.Deriv(i) - rhs_.Deriv(i));
  }
private:
  const LHSType& lhs_;//!< Reference to constant object on left hand side
  const RHSType& rhs_;//!< Reference to constant object on right hand side
}; // End class SurrealSubtract

//----------------------------- Real - Surreal -------------------------------
//****************************************************************************80
//! \brief SurrealSubtract : A class template to represent automatic
//!                     differentiation of the subtraction of an arbitrary type
//!                     from a real type
//! \details
//! \nick
//! \tparam RHSType Type used for the argument on the right hand side of -
//****************************************************************************80
template<class RHSType, int N>
class SurrealSubtract<typename RHSType::realT_, RHSType, N> :
  public SurrealBase<SurrealSubtract<typename RHSType::realT_, RHSType, N>, N>
{
public:
  typedef typename RHSType::realT_ realT;  //get fundamental real type from rght
  typedef realT realT_;  //store real type

//****************************************************************************80
//! \brief Constructor for constructing subtraction operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to constant real   on left side of -
//! \param[in] rhs_in The reference to constant object on right side of -
//****************************************************************************80
  SurrealSubtract
  (const typename RHSType::realT_& lhs_in, const RHSType& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {}

//****************************************************************************80
//! \brief Value : Returns the value of the left hand side minus the right hand
//!                side
//! \details Value of (real - Surreal) =
//!                    real - Surreal.Value()
//! \nick
//! \return Value of lhs minus rhs
//****************************************************************************80
  inline realT Value() const {
    return(lhs_ - rhs_.Value());
  }

//****************************************************************************80
//! \brief Deriv : Returns the derivative of lhs - rhs
//! \details Derivative of (real - Surreal) =
//!                         real - Surreal.Deriv()
//! \nick
//! \param[in] i The index of the derivative you wish to compute
//! \return -rhs.Deriv(i)
//****************************************************************************80
  inline realT Deriv(const int& i) const {
    return(-rhs_.Deriv(i));
  }
private:
  const typename RHSType::realT_& lhs_; //left hande side
  const RHSType& rhs_;                  //right hand side
}; // End class SurrealSubtract

//----------------------------- Surreal - Real ---------------------------------
//****************************************************************************80
//! \brief SurrealSubtract : A class template to represent automatic
//!                     differentiation of the subtraction of real type
//!                     from an arbitrary (Surreal) type
//! \details
//! \nick
//! \tparam LHSType Type used for the argument on the left hand side of -
//****************************************************************************80
template<class LHSType, int N>
class SurrealSubtract<LHSType, typename LHSType::realT_, N> :
  public SurrealBase<SurrealSubtract<LHSType, typename LHSType::realT_, N>, N>
{
public:
  typedef typename LHSType::realT_ realT; //get fundamental real type from left
  typedef realT realT_; //store real type

//****************************************************************************80
//! \brief Constructor for constructing subtraction operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to constant real   on left side of -
//! \param[in] rhs_in The reference to constant object on right side of -
//****************************************************************************80
  SurrealSubtract
  (const LHSType& lhs_in, const typename LHSType::realT_& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {}

//****************************************************************************80
//! \brief Value : Returns the value of the left hand side minus the right hand
//!                side
//! \details Value of (Surreal - real) =
//!                    Surreal.Value() - real
//! \nick
//! \return Value of lhs minus rhs
//****************************************************************************80
  inline realT Value() const {
    return(lhs_.Value() - rhs_);
  }

//****************************************************************************80
//! \brief Deriv : Returns the derivative of lhs - rhs
//! \details Derivative of (Surreal - real) =
//!                         Surreal.Deriv()
//! \nick
//! \param[in] i The index you of the derivative you wish to compute
//! \return lhs.Deriv(i)
//****************************************************************************80
  inline realT Deriv(const int& i) const {
    return(lhs_.Deriv(i));
  }
private:
  const LHSType& lhs_;                  //left  hand side
  const typename LHSType::realT_& rhs_; //right hand side
}; // End class SurrealSubtract

//-------------------------- SUBTRACTION OPERATORS -----------------------------
//****************************************************************************80
//! \brief Operator - : This operator declaration declares how to subtract a
//!                     SurrealBase type from another SurrealBase type
//! \details  SurrealSubtract = SurrealBase - SurrealBase
//! \nick
//! \param[in] lhs The object on left side of - sign
//! \param[in] rhs The object on right side of - sign
//! \return SurrealSubtract object to represent subtraction of lhs - rhs
//****************************************************************************80
template<class LHSType, class RHSType, int N>
inline SurrealSubtract<LHSType, RHSType, N>
operator-(const SurrealBase<LHSType, N>& lhs,
          const SurrealBase<RHSType, N>& rhs);

//****************************************************************************80
//! \brief Operator - : This operator declaration declares how to subtract a
//!                     SurrealBase type from a real
//! \details SurrealSubtract = real - SurrealBase
//! \nick
//! \param[in] lhs The real number on left side of - sign
//! \param[in] rhs The object on right side of - sign
//! \return SurrealSubtract object to represent subtraction of lhs - rhs
//****************************************************************************80
template<class RHSType, int N>
inline SurrealSubtract<typename RHSType::realT_, RHSType, N>
operator-
(const typename RHSType::realT_& lhs, const SurrealBase<RHSType, N>& rhs);

//****************************************************************************80
//! \brief Operator - : This operator declaration declares how to subtract a
//!                     real from a SurrealBase type
//! \details SurrealSubtract = real - SurrealBase
//! \nick
//! \param[in] lhs The object on left side of - sign
//! \param[in] rhs The real number on right side of - sign
//! \return SurrealSubtract object to represent subtraction of lhs - rhs
//****************************************************************************80
template<class LHSType, int N>
inline SurrealSubtract<LHSType, typename LHSType::realT_, N>
operator-
(const SurrealBase<LHSType, N>& lhs, const typename LHSType::realT_& rhs);
#include "SurrealSubtract_Imple.h"
#endif
