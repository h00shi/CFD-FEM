#ifndef SURREALRMULTIPLY_H
#define SURREALRMULTIPLY_H

#include "Surreal/Reverse/SurrealRBase.h"
#include <type_traits>

//----------------------------- Surreal * Surreal ------------------------------
//****************************************************************************80
//! \brief SurrealRMultiply : A class template to represent automatic
//!                     differentiation of the multiplication operator
//!                     on two arbitrary (SurrealR) types
//! \details
//! \nick
//! \tparam LHSType Type used for the arguments on the left hand side of *
//! \tparam RHSType Type used for the arguments on the right hand side of *
//****************************************************************************80
template<class LHSType, class RHSType>
class SurrealRMultiply :
  public SurrealRBase<SurrealRMultiply<LHSType, RHSType>,
                       typename RHSType::realT_, RHSType::N_>
{
private:
  typedef typename RHSType::realT_ realT;  //get ease of use
  LHSType const& lhs_; //!< Reference to object on left hand side
  RHSType const& rhs_; //!< Reference to object on right hand side
public:
//****************************************************************************80
//! \brief Constructor for constructing multiply operation from reference to
//!        left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to object on left side of *
//! \param[in] rhs_in The reference to object on right side of *
//****************************************************************************80
  SurrealRMultiply(LHSType const& lhs_in, RHSType const& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {
    static_assert(std::is_same<typename LHSType::realT_,
                  typename RHSType::realT_>::value,
                  "Surreal binary operations require the same floating-point "
                  "data type on left and right sides");
    static_assert(LHSType::N_ == RHSType::N_,
                  "Surreal binary operations require the same number of "
                  "derivatives on left and right sides");

    this->value_ = lhs_.GetValue() * rhs_.GetValue();
  }

//****************************************************************************80
//! \brief Diff : Calculate adjoints of left and right sides for the
//!               multiplication operator and calls Diff on both sides
//! \details
//! \nick
//****************************************************************************80
  void Diff(realT adjoint,
            realT (&deriv)[RHSType::N_], int (&deriv_ix)[RHSType::N_],
            int (&ix)[RHSType::N_], unsigned int & count) const {
    lhs_.Diff(adjoint*rhs_.GetValue(),
              deriv, deriv_ix, ix, count);
    rhs_.Diff(adjoint*lhs_.GetValue(),
              deriv, deriv_ix, ix, count);
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
class SurrealRMultiply<typename RHSType::realT_, RHSType> :
  public SurrealRBase<SurrealRMultiply<typename RHSType::realT_, RHSType>,
                       typename RHSType::realT_, RHSType::N_>
{
private:
  typedef typename RHSType::realT_ realT;  //get fundamental real type from rght
  const realT lhs_; //left  hand side
  RHSType const& rhs_;    //right hand side

public:
//****************************************************************************80
//! \brief Constructor for constructing multiply operation from reference to
//!        constant left and right sides.
//! \details
//! \nick
//! \param[in] lhs_in The reference to constant real   on left side of *
//! \param[in] rhs_in The reference to constant object on right side of *
//****************************************************************************80
  SurrealRMultiply(const realT lhs_in, RHSType const& rhs_in) :
    lhs_(lhs_in), rhs_(rhs_in) {
    this->value_ = lhs_ * rhs_.GetValue();
  }

//****************************************************************************80
//! \brief Diff : Calculate adjoints of right side for the
//!               multiplication operator and calls Diff
//! \details
//! \nick
//! \return ith Derivative of (Surreal1 * Surreal2)
//****************************************************************************80
  void Diff(realT adjoint,
            realT (&deriv)[RHSType::N_], int (&deriv_ix)[RHSType::N_],
            int (&ix)[RHSType::N_], unsigned int & count) const {
    rhs_.Diff(adjoint*lhs_, deriv, deriv_ix, ix, count);
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
inline SurrealRMultiply<LHSType, RHSType>
operator*(SurrealRBase<LHSType, realT, N> const & lhs,
          SurrealRBase<RHSType, realT, N> const & rhs)
{
  return(SurrealRMultiply<LHSType,RHSType>
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
inline SurrealRMultiply<typename RHSType::realT_, RHSType>
operator*(const typename RHSType::realT_ lhs,
          SurrealRBase<RHSType, typename RHSType::realT_, RHSType::N_> const & rhs)
{
  return(SurrealRMultiply<typename RHSType::realT_, RHSType>
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
inline SurrealRMultiply<typename LHSType::realT_, LHSType>
operator*(SurrealRBase<LHSType, typename LHSType::realT_, LHSType::N_> const & lhs,
          const typename LHSType::realT_ rhs)
{
  return(SurrealRMultiply<typename LHSType::realT_, LHSType>
         (rhs, lhs.CastToDerived()));
}



#endif
