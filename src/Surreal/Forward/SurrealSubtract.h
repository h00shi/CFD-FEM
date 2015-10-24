#ifndef SURREALSUBTRACT_H
#define SURREALSUBTRACT_H

#include "Surreal/Forward/SurrealBase.h"
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
template<class LHSType, class RHSType>
class SurrealSubtract :
  public SurrealBase<SurrealSubtract<LHSType, RHSType>,
                     typename RHSType::realT_, RHSType::N_>
{
private:
  typedef typename RHSType::realT_ realT;  //get fundamental real type from rght
  const LHSType& lhs_;//!< Reference to constant object on left hand side
  const RHSType& rhs_;//!< Reference to constant object on right hand side
public:
//****************************************************************************80
//! \brief Constructor for constructing subtraction operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to constant object on left side of -
//! \param[in] rhs_in The reference to constant object on right side of -
//****************************************************************************80
  SurrealSubtract(const LHSType& lhs_in, const RHSType& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {
    static_assert(std::is_same<typename LHSType::realT_,
                  typename RHSType::realT_>::value,
                  "Surreal binary operations require the same floating-point "
                  "data type on left and right sides");
    static_assert(LHSType::N_ == RHSType::N_,
                  "Surreal binary operations require the same number of "
                  "derivatives on left and right sides");
    this->value_ = lhs_.Value() - rhs_.Value();
  }
//****************************************************************************80
//! \brief Deriv : Returns the derivative of lhs - rhs
//! \details Derivative of (Surreal1 - Surreal2) =
//!                         Surreal1.Deriv() - Surreal2.Deriv()
//! \nick
//! \param[in] i The index of the derivative you wish to compute
//! \return ith derivative of lhs minus rhs
//****************************************************************************80
   inline realT Deriv(const int i) const {
    return(lhs_.Deriv(i) - rhs_.Deriv(i));
  }
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
template<class RHSType>
class SurrealSubtract<typename RHSType::realT_, RHSType> :
  public SurrealBase<SurrealSubtract<typename RHSType::realT_, RHSType>,
                     typename RHSType::realT_, RHSType::N_>
{
private:
  typedef typename RHSType::realT_ realT;  //get fundamental real type from rght
  const realT lhs_;    //left hand side
  const RHSType& rhs_; //right hand side
public:
//****************************************************************************80
//! \brief Constructor for constructing subtraction operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to constant real   on left side of -
//! \param[in] rhs_in The reference to constant object on right side of -
//****************************************************************************80
  SurrealSubtract(const realT lhs_in, const RHSType& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {
    this->value_ = lhs_ - rhs_.Value();
  }
//****************************************************************************80
//! \brief Deriv : Returns the derivative of lhs - rhs
//! \details Derivative of (real - Surreal) =
//!                         real - Surreal.Deriv()
//! \nick
//! \param[in] i The index of the derivative you wish to compute
//! \return -rhs.Deriv(i)
//****************************************************************************80
  inline realT Deriv(const int i) const {
    return(-rhs_.Deriv(i));
  }
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
template<class LHSType>
class SurrealSubtract<LHSType, typename LHSType::realT_> :
  public SurrealBase<SurrealSubtract<LHSType, typename LHSType::realT_>,
                     typename LHSType::realT_, LHSType::N_>
{
private:
  typedef typename LHSType::realT_ realT; //get fundamental real type from left
  const LHSType& lhs_; //left  hand side
  const realT rhs_;    //right hand side
public:
//****************************************************************************80
//! \brief Constructor for constructing subtraction operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to constant real   on left side of -
//! \param[in] rhs_in The reference to constant object on right side of -
//****************************************************************************80
  SurrealSubtract
  (const LHSType& lhs_in, const realT rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {
    this->value_ = lhs_.Value() - rhs_;
  }
//****************************************************************************80
//! \brief Deriv : Returns the derivative of lhs - rhs
//! \details Derivative of (Surreal - real) =
//!                         Surreal.Deriv()
//! \nick
//! \param[in] i The index you of the derivative you wish to compute
//! \return lhs.Deriv(i)
//****************************************************************************80
  inline realT Deriv(const int i) const {
    return(lhs_.Deriv(i));
  }
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
template<class LHSType, class RHSType, class realT, int N>
inline SurrealSubtract<LHSType, RHSType>
operator-(const SurrealBase<LHSType, realT, N>& lhs,
          const SurrealBase<RHSType, realT, N>& rhs){
  return(SurrealSubtract<LHSType,RHSType>
         (lhs.CastToDerived(), rhs.CastToDerived()));
}

//****************************************************************************80
//! \brief Operator - : This operator declaration declares how to subtract a
//!                     SurrealBase type from a real
//! \details SurrealSubtract = real - SurrealBase
//! \nick
//! \param[in] lhs The real number on left side of - sign
//! \param[in] rhs The object on right side of - sign
//! \return SurrealSubtract object to represent subtraction of lhs - rhs
//****************************************************************************80
template<class RHSType>
inline SurrealSubtract<typename RHSType::realT_, RHSType>
operator-
(const typename RHSType::realT_ lhs,
 const SurrealBase<RHSType, typename RHSType::realT_, RHSType::N_>& rhs){
  return(SurrealSubtract<typename RHSType::realT_,RHSType>
         (lhs, rhs.CastToDerived()));
}

//****************************************************************************80
//! \brief Operator - : This operator declaration declares how to subtract a
//!                     real from a SurrealBase type
//! \details SurrealSubtract = SurrealBase - real
//! \nick
//! \param[in] lhs The object on left side of - sign
//! \param[in] rhs The real number on right side of - sign
//! \return SurrealSubtract object to represent subtraction of lhs - rhs
//****************************************************************************80
template<class LHSType>
inline SurrealSubtract<LHSType, typename LHSType::realT_>
operator-
(const SurrealBase<LHSType, typename LHSType::realT_, LHSType::N_>& lhs,
 const typename LHSType::realT_ rhs){
  return(SurrealSubtract<LHSType, typename LHSType::realT_>
         (lhs.CastToDerived(), rhs));
}
#endif
