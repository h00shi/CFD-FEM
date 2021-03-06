// -*-c++-*-
#ifndef SURREALPOW_H
#define SURREALPOW_H

#include "Surreal/SurrealBase.h"

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
template<class LHSType, class RHSType, int N>
class SurrealPow :  public SurrealBase<SurrealPow<LHSType, RHSType, N>, N>
{
public:
  typedef typename LHSType::realT_ realT; //get fundamental real type from left
  typedef realT realT_; //store real type

//****************************************************************************80
//! \brief Constructor for constructing pow operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to const object on "left hand side" of pow
//! \param[in] rhs_in The reference to const object on "right hand side" of pow
//****************************************************************************80
  SurrealPow(const LHSType& lhs_in, const RHSType& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {}

//****************************************************************************80
//! \brief Value : Returns the value of the power operator on the left and
//!                right hand sides
//! \details Value of pow(Surreal1, Surreal2) =
//!                std::pow(Surreal1.value(), Surreal2.Value())
//! \nick
//! \return pow(lhs.Value(), rhs.Value())
//****************************************************************************80
  inline realT Value() const {
    return(std::pow(lhs_.Value(), rhs_.Value()));
  }
//****************************************************************************80
//! \brief Deriv : Returns the derivative of pow(lhs,rhs)
//! \details Derivative of pow(Surreal1,Surreal2) =
//! \nick
//! \param[in] i The index of the derivative you wish to compute
//! \return Derivative of pow(lhs,rhs)
//****************************************************************************80
 inline realT Deriv(const int& i) const {
   return(std::pow(lhs_.Value(), rhs_.Value())*
          ( rhs_.Deriv(i) * std::log(lhs_.Value()) +
            rhs_.Value()/lhs_.Value() * lhs_.Deriv(i) ) );
 }
private:
  const LHSType & lhs_; //!< Reference to constant object on left hand side
  const RHSType & rhs_; //!< Reference to constant object on right hand side
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
template<class RHSType, int N>
class SurrealPow<typename RHSType::realT_, RHSType, N> :
  public SurrealBase<SurrealPow<typename RHSType::realT_, RHSType, N>, N>
{
public:
  typedef typename RHSType::realT_ realT;  //get fundamental real type
  typedef realT realT_;  //store real type
//****************************************************************************80
//! \brief Constructor for constructing add operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to const real   on "left hand side" of pow
//! \param[in] rhs_in The reference to const object on "right hand side" of pow
//****************************************************************************80
  SurrealPow(const typename RHSType::realT_& lhs_in, const RHSType& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {}

//****************************************************************************80
//! \brief Value : Returns the value of the power operator on the left and
//!                right hand sides
//! \details Value of pow(real, Surreal) =
//!                std::pow(real,Surreal.Value())
//! \nick
//! \return pow(lhs, rhs.Value())
//****************************************************************************80
  inline realT Value() const {
    return(std::pow(lhs_,rhs_.Value()));
  }
//****************************************************************************80
//! \brief Deriv : Returns the derivative of pow(lhs,rhs)
//! \details Derivative of pow(real,Surreal) =
//!                        Surreal.Value()*std::log(real)*Surreal.Deriv()
//! \nick
//! \param[in] i The index of the derivative you wish to compute
//! \return Derivative of pow(lhs,rhs)
//****************************************************************************80
 inline realT Deriv(const int& i) const {
   return rhs_.Value()*std::log(lhs_)*rhs_.Deriv(i);
 }
private:
  const typename RHSType::realT_ & lhs_; //left  hand side
  const RHSType& rhs_;                   //right hand side
};

//---> pow(Surreal, real)
//****************************************************************************80
//! \brief SurrealPow : A class template to represent automatic
//!                     differentiation of a surreal raised to a real
//! \details pow(Surreal, real)
//! \nick
//! \tparam LHSType Type used for the argument on the "left hand side" of pow
//****************************************************************************80
template<class LHSType, int N>
class SurrealPow<LHSType, typename LHSType::realT_, N> :
  public SurrealBase<SurrealPow<LHSType, typename LHSType::realT_,N>,N>
{
public:
  typedef typename LHSType::realT_ realT;  //get fundamental real type
  typedef realT realT_;  //store real type

//****************************************************************************80
//! \brief Constructor for constructing add operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to const object on "left hand side" of pow
//! \param[in] rhs_in The reference to const real   on "right hand side" of pow
//****************************************************************************80
  SurrealPow(const LHSType& lhs_in, const typename LHSType::realT_& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {}

//****************************************************************************80
//! \brief Value : Returns the value of the power operator on the left and
//!                right hand sides
//! \details Value of pow(Surreal, real) =
//!                   pow(Surreal.Value(), real)
//! \nick
//! \return pow(lhs.Value(), rhs.Value())
//****************************************************************************80
  inline realT Value() const {
    return(std::pow(lhs_.Value(), rhs_));
  }
//****************************************************************************80
//! \brief Deriv : Returns the derivative of pow(lhs,rhs)
//! \details Derivative of pow(Surreal, real) =
//!           real * pow(Surreal.Value, real-1) * Surreal.Deriv()
//! \nick
//! \param[in] i The index of the derivative you wish to compute
//! \return Derivative of pow(lhs,rhs)
//****************************************************************************80
 inline realT Deriv(const int& i) const {
      return rhs_*std::pow(lhs_.Value(), rhs_-1)*lhs_.Deriv(i);
 }
private:
  const LHSType & lhs_;                   //left  hand side
  const typename LHSType::realT_  & rhs_; //right hand side
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
template<class LHSType, class RHSType, int N>
inline SurrealPow<LHSType, RHSType, N>
pow(const SurrealBase<LHSType,N>& lhs, const SurrealBase<RHSType, N>& rhs);

//****************************************************************************80
//! \brief pow : This operator declaration declares how to raise a real
//!              number to the power of a SurrealBase type
//! \details SurrealPow = pow(real, SurrealBase)
//! \nick
//! \param[in] lhs The real number on the "left hand side" of pow
//! \param[in] rhs The object on the "right hand side" of pow
//! \return SurrealPow object to represent pow(lhs,rhs)
//****************************************************************************80
template<class RHSType, int N>
inline SurrealPow<typename RHSType::realT_, RHSType, N>
pow(const typename RHSType::realT_& lhs, const SurrealBase<RHSType, N>& rhs);

//****************************************************************************80
//! \brief pow : This operator declaration declares how to raise a surreal
//!              type to the power of a real number
//! \details SurrealPow = pow(SurrealBase, real)
//! \nick
//! \param[in] lhs The object on the "left hand side" of pow
//! \param[in] rhs The real number on the "right hand side" of pow
//! \return SurrealPow object to represent pow(lhs,rhs)
//****************************************************************************80
template<class LHSType, int N>
inline SurrealPow<LHSType, typename LHSType::realT_, N>
pow(const SurrealBase<LHSType, N>& lhs, const typename LHSType::realT_ & rhs);

#include "SurrealPow_Imple.h"
#endif

