#ifndef SURREALRSUBTRACT_H
#define SURREALRSUBTRACT_H

#include "Surreal/Reverse/SurrealRBase.h"
#include <type_traits>

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
class SurrealRSubtract :
  public SurrealRBase<SurrealRSubtract<LHSType, RHSType>,
                       typename RHSType::realT_, RHSType::N_>
{
private:
  typedef typename RHSType::realT_ realT; //get fundamental real type from left
  LHSType const& lhs_;//!< Reference to constant object on left hand side
  RHSType const& rhs_;//!< Reference to constant object on right hand side
public:
//****************************************************************************80
//! \brief Constructor for constructing subtraction operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to constant object on left side of -
//! \param[in] rhs_in The reference to constant object on right side of -
//****************************************************************************80
  SurrealRSubtract(LHSType const& lhs_in, RHSType const& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {
    static_assert(std::is_same<typename LHSType::realT_,
                  typename RHSType::realT_>::value,
                  "Surreal binary operations require the same floating-point "
                  "data type on left and right sides");
    static_assert(LHSType::N_ == RHSType::N_,
                  "Surreal binary operations require the same number of "
                  "derivatives on left and right sides");
    SurrealRBase<SurrealRSubtract<LHSType, RHSType>,
   typename RHSType::realT_, RHSType::N_>::value_ =
       lhs_.GetValue() - rhs_.GetValue();
  }

//****************************************************************************80
//! \brief Diff : Calculate adjoints of left and right sides for the
//!               subtraction operator and calls Diff on both sides
//! \details
//! \nick
//****************************************************************************80
  void Diff(realT adjoint,
            realT (&deriv)[RHSType::N_], int (&deriv_ix)[RHSType::N_],
            int (&ix)[RHSType::N_], unsigned int & count) const {
    lhs_.Diff( adjoint, deriv, deriv_ix, ix, count);
    rhs_.Diff(-adjoint, deriv, deriv_ix, ix, count);
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
class SurrealRSubtract<typename RHSType::realT_, RHSType> :
  public SurrealRBase<SurrealRSubtract<typename RHSType::realT_, RHSType>,
                       typename RHSType::realT_, RHSType::N_>
{
private:
  typedef typename RHSType::realT_ realT;  //get fundamental real type from rght
  const realT lhs_; //left hande side
  RHSType const& rhs_;    //right hand side
public:
//****************************************************************************80
//! \brief Constructor for constructing subtraction operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to constant real   on left side of -
//! \param[in] rhs_in The reference to constant object on right side of -
//****************************************************************************80
  SurrealRSubtract(const realT lhs_in, RHSType const& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {
    SurrealRBase<SurrealRSubtract<typename RHSType::realT_, RHSType>,
   typename RHSType::realT_, RHSType::N_>::value_ = lhs_ - rhs_.GetValue();
  }

//****************************************************************************80
//! \brief Diff : Calculate adjoints of right side for the
//!               subtraction operator and calls Diff
//! \details
//! \nick
//****************************************************************************80
  void Diff(realT adjoint,
            realT (&deriv)[RHSType::N_], int (&deriv_ix)[RHSType::N_],
            int (&ix)[RHSType::N_], unsigned int & count) const{
    rhs_.Diff(-adjoint, deriv, deriv_ix, ix, count);
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
class SurrealRSubtract<LHSType, typename LHSType::realT_> :
  public SurrealRBase<SurrealRSubtract<LHSType, typename LHSType::realT_>,
                       typename LHSType::realT_, LHSType::N_>
{
private:
  typedef typename LHSType::realT_ realT; //get fundamental real type from left
  LHSType const& lhs_;    //left  hand side
  const realT rhs_; //right hand side
public:
//****************************************************************************80
//! \brief Constructor for constructing subtraction operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to constant real   on left side of -
//! \param[in] rhs_in The reference to constant object on right side of -
//****************************************************************************80
  SurrealRSubtract
  (LHSType const& lhs_in, const realT rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {
    SurrealRBase<SurrealRSubtract<LHSType, typename LHSType::realT_>,
   typename LHSType::realT_, LHSType::N_>::value_ = lhs_.GetValue() - rhs_;
  }

//****************************************************************************80
//! \brief Diff : Calculate adjoints of left side for the
//!               subtraction operator and calls Diff
//! \details
//! \nick
//****************************************************************************80
  void Diff(realT adjoint,
            realT (&deriv)[LHSType::N_], int (&deriv_ix)[LHSType::N_],
            int (&ix)[LHSType::N_], unsigned int & count) const{
    lhs_.Diff(adjoint, deriv, deriv_ix, ix, count);
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
inline SurrealRSubtract<LHSType, RHSType>
operator-(SurrealRBase<LHSType, realT, N> const& lhs,
          SurrealRBase<RHSType, realT, N> const& rhs){
  return(SurrealRSubtract<LHSType,RHSType>
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
inline SurrealRSubtract<typename RHSType::realT_, RHSType>
operator-
(const typename RHSType::realT_ lhs,
 SurrealRBase<RHSType, typename RHSType::realT_, RHSType::N_> const& rhs){
  return(SurrealRSubtract<typename RHSType::realT_,RHSType>
         (lhs, rhs.CastToDerived()));
}

//****************************************************************************80
//! \brief Operator - : This operator declaration declares how to subtract a
//!                     real from a SurrealBase type
//! \details SurrealSubtract = real - SurrealBase
//! \nick
//! \param[in] lhs The object on left side of - sign
//! \param[in] rhs The real number on right side of - sign
//! \return SurrealSubtract object to represent subtraction of lhs - rhs
//****************************************************************************80
template<class LHSType>
inline SurrealRSubtract<LHSType, typename LHSType::realT_>
operator-
(SurrealRBase<LHSType, typename LHSType::realT_, LHSType::N_> const & lhs,
 const typename LHSType::realT_ rhs){
  return(SurrealRSubtract<LHSType, typename LHSType::realT_>
         (lhs.CastToDerived(), rhs));
}

#endif
