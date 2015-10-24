#ifndef SURREALRPOW_H
#define SURREALRPOW_H

#include "Surreal/Reverse/SurrealRBase.h"
#include <type_traits>
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
class SurrealRPow : public SurrealRBase<SurrealRPow<LHSType, RHSType>,
                                          typename RHSType::realT_, RHSType::N_>
{
private:
  typedef typename RHSType::realT_ realT;  //get fundamental real type from rght
  LHSType const& lhs_; //!< Reference to constant object on left hand side
  RHSType const& rhs_; //!< Reference to constant object on right hand side

public:

//****************************************************************************80
//! \brief Constructor for constructing pow operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to const object on "left hand side" of pow
//! \param[in] rhs_in The reference to const object on "right hand side" of pow
//****************************************************************************80
  SurrealRPow(LHSType const& lhs_in, RHSType const& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {
    static_assert(std::is_same<typename LHSType::realT_,
                  typename RHSType::realT_>::value,
                  "Surreal binary operations require the same floating-point "
                  "data type on left and right sides");
    static_assert(LHSType::N_ == RHSType::N_,
                  "Surreal binary operations require the same number of "
                  "derivatives on left and right sides");

    this->value_ = std::pow(lhs_.GetValue(), rhs_.GetValue());
  }

//****************************************************************************80
//! \brief Diff : Calculate adjoints of left and right sides for the
//!               pow operator and calls Diff on both sides
//! \details
//! \nick
//****************************************************************************80
  void Diff(realT adjoint,
            realT (&deriv)[RHSType::N_], int (&deriv_ix)[RHSType::N_],
            int (&ix)[RHSType::N_], unsigned int & count) const{
    lhs_.Diff(adjoint *rhs_.GetValue()*
              std::pow(lhs_.GetValue(), rhs_.GetValue()-1.0),
              deriv, deriv_ix, ix, count);
    rhs_.Diff(adjoint* this->value_ * std::log(lhs_.GetValue()),
              deriv, deriv_ix, ix, count);
  }
};

//---> pow(real, Surreal)
//****************************************************************************80
//! \brief SurrealRPow : A class template to represent automatic
//!                     differentiation of a real raised to
//!                     a surreal type
//! \details pow(real,Surreal)
//! \nick
//! \tparam RHSType Type used for the argument on the "right hand side" of pow
//****************************************************************************80
template<class RHSType>
class SurrealRPow<typename RHSType::realT_, RHSType> :
  public SurrealRBase<SurrealRPow<typename RHSType::realT_, RHSType>,
                       typename RHSType::realT_, RHSType:: N_>
{
private:
  typedef typename RHSType::realT_ realT;  //get fundamental real type
  const realT lhs_; //left  hand side
  RHSType const& rhs_;                   //right hand side
public:
//****************************************************************************80
//! \brief Constructor for constructing add operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to const real   on "left hand side" of pow
//! \param[in] rhs_in The reference to const object on "right hand side" of pow
//****************************************************************************80
  SurrealRPow(const typename RHSType::realT_ lhs_in, RHSType const& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {
    this->value_ = std::pow(lhs_, rhs_.GetValue());
  }
//****************************************************************************80
//! \brief Diff : Calculate adjoints of right side for the
//!               pow operator and calls Diff
//! \details
//! \nick
//****************************************************************************80
  void Diff(realT adjoint,
            realT (&deriv)[RHSType::N_], int (&deriv_ix)[RHSType::N_],
            int (&ix)[RHSType::N_], unsigned int & count) const {
    rhs_.Diff(adjoint* this->value_ * std::log(lhs_),
              deriv, deriv_ix, ix, count);
  }

};

//---> pow(Surreal, real)
//****************************************************************************80
//! \brief SurrealRPow : A class template to represent automatic
//!                     differentiation of a surreal raised to a real
//! \details pow(Surreal, real)
//! \nick
//! \tparam LHSType Type used for the argument on the "left hand side" of pow
//****************************************************************************80
template<class LHSType>
class SurrealRPow<LHSType, typename LHSType::realT_> :
  public SurrealRBase<SurrealRPow<LHSType, typename LHSType::realT_>,
                       typename LHSType::realT_, LHSType::N_>
{
private:
  typedef typename LHSType::realT_ realT;  //get fundamental real type
  LHSType const & lhs_;                   //left  hand side
  const realT rhs_; //right hand side

public:
//****************************************************************************80
//! \brief Constructor for constructing add operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to const object on "left hand side" of pow
//! \param[in] rhs_in The reference to const real   on "right hand side" of pow
//****************************************************************************80
  SurrealRPow(LHSType const & lhs_in, const typename LHSType::realT_ rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {
    this->value_ = std::pow(lhs_.GetValue(), rhs_);
  }

//****************************************************************************80
//! \brief Diff : Calculate adjoint of left side for the
//!               pow operator and calls Diff
//! \details
//! \nick
//****************************************************************************80
  void Diff(realT adjoint,
            realT (&deriv)[LHSType::N_], int (&deriv_ix)[LHSType::N_],
            int (&ix)[LHSType::N_], unsigned int & count) const {
    lhs_.Diff(adjoint *rhs_ *std::pow(lhs_.GetValue(), rhs_-1.0),
              deriv, deriv_ix, ix, count);
  }
};

//------------------------------ POW OPERATORS --------------------------------
//****************************************************************************80
//! \brief pow : This operator declaration declares how to raise a SurrealBase
//!              type to the power of a SurrealBase type
//! \details SurrealRPow = pow(real, SurrealBase)
//! \nick
//! \param[in] lhs The object on the "left hand side" of pow
//! \param[in] rhs The object on the "right hand side" of pow
//! \return SurrealRPow object to represent pow(lhs,rhs)
//****************************************************************************80
template<class LHSType, class RHSType, class realT, int N>
inline SurrealRPow<LHSType, RHSType>
pow(SurrealRBase<LHSType, realT, N> const & lhs,
    SurrealRBase<RHSType, realT, N> const & rhs){
  return(SurrealRPow<LHSType,RHSType>
         (lhs.CastToDerived(), rhs.CastToDerived()));
}
//****************************************************************************80
//! \brief pow : This operator declaration declares how to raise a real
//!              number to the power of a SurrealBase type
//! \details SurrealRPow = pow(real, SurrealBase)
//! \nick
//! \param[in] lhs The real number on the "left hand side" of pow
//! \param[in] rhs The object on the "right hand side" of pow
//! \return SurrealRPow object to represent pow(lhs,rhs)
//****************************************************************************80
template<class RHSType>
inline SurrealRPow<typename RHSType::realT_, RHSType>
pow(const typename RHSType::realT_ lhs,
    SurrealRBase<RHSType, typename RHSType::realT_, RHSType::N_> const & rhs){
  return(SurrealRPow<typename RHSType::realT_,RHSType>
         (lhs, rhs.CastToDerived()));
}

//****************************************************************************80
//! \brief pow : This operator declaration declares how to raise a surreal
//!              type to the power of a real number
//! \details SurrealRPow = pow(SurrealBase, real)
//! \nick
//! \param[in] lhs The object on the "left hand side" of pow
//! \param[in] rhs The real number on the "right hand side" of pow
//! \return SurrealRPow object to represent pow(lhs,rhs)
//****************************************************************************80
template<class LHSType>
inline SurrealRPow<LHSType, typename LHSType::realT_>
pow(SurrealRBase<LHSType, typename LHSType::realT_, LHSType::N_> const & lhs,
    const typename LHSType::realT_ rhs){
  return(SurrealRPow<LHSType, typename LHSType::realT_>
         (lhs.CastToDerived(), rhs));
}

#endif

