#ifndef SURREALRADD_H
#define SURREALRADD_H

#include "Surreal/Reverse/SurrealRBase.h"
#include <type_traits>

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
template<class LHSType, class RHSType>
class SurrealRAdd : public SurrealRBase<SurrealRAdd<LHSType,RHSType>,
                                          typename RHSType::realT_, RHSType::N_>
{
private:
  typedef typename RHSType::realT_ realT;  //get fundamental real type from rght
  LHSType const& lhs_; //!< Reference to object on left hand side
  RHSType const& rhs_; //!< Reference to object on right hand side

public:
//****************************************************************************80
//! \brief Constructor for constructing add operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to constant object on left side of +
//! \param[in] rhs_in The reference to constant object on right side of +
//****************************************************************************80
  SurrealRAdd(LHSType const & lhs_in, RHSType const & rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {

    static_assert(std::is_same<typename LHSType::realT_,
                  typename RHSType::realT_>::value,
                  "Surreal binary operations require the same floating-point "
                  "data type on left and right sides");

    static_assert(LHSType::N_ == RHSType::N_,
                  "Surreal binary operations require the same number of "
                  "derivatives on left and right sides");

    SurrealRBase<SurrealRAdd<LHSType,RHSType>,
    typename RHSType::realT_, RHSType::N_>::value_ =
    lhs_.GetValue() + rhs_.GetValue();
  }

//****************************************************************************80
//! \brief Diff : Calculate adjoints of left and right sides for the
//!               addition operator and calls Diff on both sides
//! \details
//! \nick
//****************************************************************************80
  void Diff(realT adjoint,
            realT (&deriv)[RHSType::N_], int (&deriv_ix)[RHSType::N_],
            int (&ix)[RHSType::N_], unsigned int & count) const {
    rhs_.Diff(adjoint, deriv, deriv_ix, ix, count);
    lhs_.Diff(adjoint, deriv, deriv_ix, ix, count);
  }
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
template<class RHSType>
class SurrealRAdd<typename RHSType::realT_, RHSType> :
  public SurrealRBase<SurrealRAdd<typename RHSType::realT_, RHSType>,
                       typename RHSType::realT_, RHSType::N_>
{
private:
  typedef typename RHSType::realT_ realT;  //get fundamental real type from rght
  const realT lhs_; //left  hand side
  RHSType const & rhs_;    //right hand side
public:
//****************************************************************************80
//! \brief Constructor for constructing add operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to constant real   on left side of +
//! \param[in] rhs_in The reference to object on right side of +
//****************************************************************************80
  SurrealRAdd(const realT lhs_in, RHSType const& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {

    SurrealRBase<SurrealRAdd<typename RHSType::realT_, RHSType>,
       typename RHSType::realT_, RHSType::N_>::value_ = lhs_ + rhs_.GetValue();

  }

//****************************************************************************80
//! \brief Diff : Calculate adjoints of right side for the
//!               addition operator and calls right side Diff
//! \details
//! \nick
//****************************************************************************80
  void Diff(realT adjoint,
            realT (&deriv)[RHSType::N_], int (&deriv_ix)[RHSType::N_],
            int (&ix)[RHSType::N_], unsigned int & count) const {
    rhs_.Diff(adjoint, deriv, deriv_ix, ix, count);
  }
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
template<class LHSType, class RHSType, class realT, int N>
inline SurrealRAdd<LHSType, RHSType>
operator+(SurrealRBase<LHSType, realT, N> const &  lhs,
          SurrealRBase<RHSType, realT, N> const & rhs)
{
  return(SurrealRAdd<LHSType, RHSType>
         (lhs.CastToDerived(), rhs.CastToDerived()));
}

//****************************************************************************80
//! \brief Operator + : This operator declaration declares how to add a real
//!                     number with a SurrealBase type
//! \details SurrealAdd = real + SurrealBase
//! \nick
//! \param[in] lhs The real number on left side of + sign
//! \param[in] rhs The object on right side of + sign
//! \return SurrealAdd object to represent addition of lhs + rhs
//****************************************************************************80
template<class RHSType>
inline SurrealRAdd<typename RHSType::realT_, RHSType>
operator+
(const typename RHSType::realT_ lhs,
 SurrealRBase<RHSType, typename RHSType::realT_, RHSType::N_> const & rhs)
{
  return(SurrealRAdd<typename RHSType::realT_, RHSType>
         (lhs, rhs.CastToDerived()));
}


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
template<class LHSType>
inline SurrealRAdd<typename LHSType::realT_, LHSType>
operator+
(SurrealRBase<LHSType, typename LHSType::realT_, LHSType::N_> const & lhs,
 const typename LHSType::realT_ rhs)
{
  return(SurrealRAdd<typename LHSType::realT_, LHSType>
         (rhs, lhs.CastToDerived()));
}
#endif

