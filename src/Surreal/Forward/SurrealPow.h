#ifndef SURREALPOW_H
#define SURREALPOW_H

#include "Surreal/Forward/SurrealBase.h"
#include <cmath>

//----------------------------- SurrealPow -------------------------------------
//---> pow(Surreal, Surreal)
//****************************************************************************80
//! \brief SurrealPow : A class template to represent automatic
//!                     differentiation of a surreal type raised to
//!                     a surreal type
//! \details
//! \nick
//! \tparam LHSType Type used for the argument on the left hand side of pow
//! \tparam RHSType Type used for the argument on the right hand side of pow
//****************************************************************************80
template<class LHSType, class RHSType>
class SurrealPow :  public SurrealBase<SurrealPow<LHSType, RHSType>,
                                       typename RHSType::realT_, RHSType::N_>
{
private:
  typedef typename RHSType::realT_ realT;  //get fundamental real type from rght
  const LHSType& lhs_;//!< Reference to constant object on left hand side
  const RHSType& rhs_;//!< Reference to constant object on right hand side
public:
//****************************************************************************80
//! \brief Constructor for constructing pow operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to const object on "left hand side" of pow
//! \param[in] rhs_in The reference to const object on "right hand side" of pow
//****************************************************************************80
  SurrealPow(const LHSType& lhs_in, const RHSType& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {
    static_assert(std::is_same<typename LHSType::realT_,
                  typename RHSType::realT_>::value,
                  "Surreal binary operations require the same floating-point "
                  "data type on left and right sides");
    static_assert(LHSType::N_ == RHSType::N_,
                  "Surreal binary operations require the same number of "
                  "derivatives on left and right sides");
    this->value_ = std::pow(lhs_.Value(), rhs_.Value());
  }
//****************************************************************************80
//! \brief Deriv : Returns the derivative of pow(lhs,rhs)
//! \details Derivative of pow(Surreal1,Surreal2) =
//! \nick
//! \param[in] i The index of the derivative you wish to compute
//! \return Derivative of pow(lhs,rhs)
//****************************************************************************80
 inline realT Deriv(const int i) const {
   return(std::pow(lhs_.Value(), rhs_.Value())*
          ( rhs_.Deriv(i) * std::log(lhs_.Value()) +
            rhs_.Value()/lhs_.Value() * lhs_.Deriv(i) ) );
 }
};

//---> pow(real, Surreal)
//****************************************************************************80
//! \brief SurrealPow : A class template to represent automatic
//!                     differentiation of a real raised to
//!                     a surreal type
//! \details pow(real,Surreal)
//! \nick
//! \tparam RHSType Type used for the argument on the "right hand side" of pow
//****************************************************************************80
template<class RHSType>
class SurrealPow<typename RHSType::realT_, RHSType> :
  public SurrealBase<SurrealPow<typename RHSType::realT_, RHSType>,
                     typename RHSType::realT_, RHSType::N_>
{
private:
  typedef typename RHSType::realT_ realT;  //get fundamental real type from rght
  const realT lhs_;    //left hand side
  const RHSType& rhs_; //right hand side
public:
//****************************************************************************80
//! \brief Constructor for constructing add operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to const real   on "left hand side" of pow
//! \param[in] rhs_in The reference to const object on "right hand side" of pow
//****************************************************************************80
  SurrealPow(const realT lhs_in, const RHSType& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {
    this->value_ = std::pow(lhs_,rhs_.Value());
  }
//****************************************************************************80
//! \brief Deriv : Returns the derivative of pow(lhs,rhs)
//! \details Derivative of pow(real,Surreal) =
//!                        Surreal.Value()*std::log(real)*Surreal.Deriv()
//! \nick
//! \param[in] i The index of the derivative you wish to compute
//! \return Derivative of pow(lhs,rhs)
//****************************************************************************80
 inline realT Deriv(const int i) const {
   return std::pow(lhs_,rhs_.Value()) *std::log(lhs_)*rhs_.Deriv(i);
 }
};

//---> pow(Surreal, real)
//****************************************************************************80
//! \brief SurrealPow : A class template to represent automatic
//!                     differentiation of a surreal raised to a real
//! \details pow(Surreal, real)
//! \nick
//! \tparam LHSType Type used for the argument on the "left hand side" of pow
//****************************************************************************80
template<class LHSType>
class SurrealPow<LHSType, typename LHSType::realT_> :
  public SurrealBase<SurrealPow<LHSType, typename LHSType::realT_>,
                     typename LHSType::realT_, LHSType::N_>

{
private:
  typedef typename LHSType::realT_ realT; //get fundamental real type from left
  const LHSType& lhs_; //left  hand side
  const realT rhs_;    //right hand side
public:
//****************************************************************************80
//! \brief Constructor for constructing add operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to const object on "left hand side" of pow
//! \param[in] rhs_in The reference to const real   on "right hand side" of pow
//****************************************************************************80
  SurrealPow(const LHSType& lhs_in, const realT rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {
    this->value_ = std::pow(lhs_.Value(), rhs_);
  }
//****************************************************************************80
//! \brief Deriv : Returns the derivative of pow(lhs,rhs)
//! \details Derivative of pow(Surreal, real) =
//!           real * pow(Surreal.Value, real-1) * Surreal.Deriv()
//! \nick
//! \param[in] i The index of the derivative you wish to compute
//! \return Derivative of pow(lhs,rhs)
//****************************************************************************80
 inline realT Deriv(const int i) const {
      return rhs_*std::pow(lhs_.Value(), rhs_-1)*lhs_.Deriv(i);
 }
};

//------------------------------ POW OPERATORS --------------------------------
//****************************************************************************80
//! \brief pow : This operator declaration declares how to raise a SurrealBase
//!              type to the power of a SurrealBase type
//! \details SurrealPow = pow(real, SurrealBase)
//! \nick
//! \param[in] lhs The object on the "left hand side" of pow
//! \param[in] rhs The object on the "right hand side" of pow
//! \return SurrealPow object to represent pow(lhs,rhs)
//****************************************************************************80
template<class LHSType, class RHSType, class realT, int N>
inline SurrealPow<LHSType, RHSType>
pow(const SurrealBase<LHSType, realT, N>& lhs,
    const SurrealBase<RHSType, realT, N>& rhs){
  return(SurrealPow<LHSType,RHSType>
         (lhs.CastToDerived(), rhs.CastToDerived()));
}

//****************************************************************************80
//! \brief pow : This operator declaration declares how to raise a real
//!              number to the power of a SurrealBase type
//! \details SurrealPow = pow(real, SurrealBase)
//! \nick
//! \param[in] lhs The real number on the "left hand side" of pow
//! \param[in] rhs The object on the "right hand side" of pow
//! \return SurrealPow object to represent pow(lhs,rhs)
//****************************************************************************80
template<class RHSType>
inline SurrealPow<typename RHSType::realT_, RHSType>
pow(const typename RHSType::realT_ lhs,
    const SurrealBase<RHSType, typename RHSType::realT_, RHSType::N_>& rhs){
  return(SurrealPow<typename RHSType::realT_,RHSType>
         (lhs, rhs.CastToDerived()));
}

//****************************************************************************80
//! \brief pow : This operator declaration declares how to raise a surreal
//!              type to the power of a real number
//! \details SurrealPow = pow(SurrealBase, real)
//! \nick
//! \param[in] lhs The object on the "left hand side" of pow
//! \param[in] rhs The real number on the "right hand side" of pow
//! \return SurrealPow object to represent pow(lhs,rhs)
//****************************************************************************80
template<class LHSType>
inline SurrealPow<LHSType, typename LHSType::realT_>
pow(const SurrealBase<LHSType, typename LHSType::realT_, LHSType::N_>& lhs,
    const typename LHSType::realT_ rhs){
  return(SurrealPow<LHSType, typename LHSType::realT_>
         (lhs.CastToDerived(), rhs));
}

#endif

