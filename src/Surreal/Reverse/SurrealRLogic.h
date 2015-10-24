#ifndef SURREALRLOGIC_H
#define SURREALRLOGIC_H

#include "Surreal/Reverse/SurrealRBase.h"
//----------------------------- LOGICAL OPERATORS ------------------------------
//----------------------------- OPERATOR < -------------------------------------
//----------------------------- Surreal < Surreal ------------------------------
//****************************************************************************80
//! \brief Operator < : This operator declaration declares how to compare
//!                     two SurrealRBase types with <
//! \nick
//! \param[in] lhs The object on left side of < sign
//! \param[in] rhs The object on right side of < sign
//! \return boolean of the comparison of lhs < rhs
//****************************************************************************80
template<class LHSType, class RHSType, class realT, int N>
inline bool operator<
(SurrealRBase<LHSType, realT, N> const & lhs,
 SurrealRBase<RHSType, realT, N> const & rhs){
  return(lhs.CastToDerived().GetValue() < rhs.CastToDerived().GetValue());
}

//----------------------------- Real < Surreal -------------------------------
//****************************************************************************80
//! \brief Operator < : This operator declaration declares how to test
//!                     if a real number is less than a SurrealRBase type
//! \nick
//! \param[in] lhs The real number on left side of < sign
//! \param[in] rhs The object on right side of < sign
//! \return boolean of the comparison of lhs < rhs
//****************************************************************************80
template<class RHSType>
inline bool operator<
(const typename RHSType::realT_ lhs,
 SurrealRBase<RHSType, typename RHSType::realT_, RHSType::N_> const& rhs){
  return(lhs < rhs.CastToDerived().GetValue());
}

//----------------------------- Surreal < Real -------------------------------
//****************************************************************************80
//! \brief Operator < : This operator declaration declares how to test
//!                     if a SurrealRBase type is less thana real number
//! \nick
//! \param[in] lhs The object on left side of < sign
//! \param[in] rhs The real number on right side of < sign
//! \return boolean of the comparison of lhs < rhs
//****************************************************************************80
template<class LHSType>
inline bool operator<
(SurrealRBase<LHSType, typename LHSType::realT_, LHSType::N_> const& lhs,
 const typename LHSType::realT_ rhs){
  return(lhs.CastToDerived().GetValue() < rhs);
}

//----------------------------- OPERATOR <= ------------------------------------
//----------------------------- Surreal <= Surreal -----------------------------
//****************************************************************************80
//! \brief Operator <= : This operator declaration declares how to compare
//!                     two SurrealRBase types with <=
//! \nick
//! \param[in] lhs The object on left side of <= sign
//! \param[in] rhs The object on right side of <= sign
//! \return boolean of the comparison of lhs <= rhs
//****************************************************************************80
 template<class LHSType, class RHSType, class realT, int N>
 inline bool operator<=
 (SurrealRBase<LHSType, realT, N> const& lhs,
  SurrealRBase<RHSType, realT, N> const& rhs){
   return(lhs.CastToDerived().GetValue() <= rhs.CastToDerived().GetValue());
 }

//----------------------------- Real <= Surreal ------------------------------
//****************************************************************************80
//! \brief Operator <= : This operator declaration declares how to test
//!                     if a real number is <= a SurrealRBase type
//! \nick
//! \param[in] lhs The real number on left side of <= sign
//! \param[in] rhs The object on right side of <= sign
//! \return boolean of the comparison of lhs <= rhs
//****************************************************************************80
template<class RHSType>
inline bool operator<=
(const typename RHSType::realT_ lhs,
 SurrealRBase<RHSType, typename RHSType::realT_, RHSType::N_> const& rhs){
  return(lhs <= rhs.CastToDerived().GetValue());
}

//----------------------------- Surreal <= Real ------------------------------
//****************************************************************************80
//! \brief Operator <= : This operator declaration declares how to test
//!                     if a SurrealRBase type is <= a real number
//! \nick
//! \param[in] lhs The object on left side of <= sign
//! \param[in] rhs The real number on right side of <= sign
//! \return boolean of the comparison of lhs <= rhs
//****************************************************************************80
template<class LHSType>
 inline bool operator<=
 (SurrealRBase<LHSType, typename LHSType::realT_, LHSType::N_> const& lhs,
  const typename LHSType::realT_ rhs){
  return(lhs.CastToDerived().GetValue() <= rhs);
}

//----------------------------- OPERATOR > ------------------------------------
//----------------------------- Surreal > Surreal -----------------------------
//****************************************************************************80
//! \brief Operator > : This operator declaration declares how to compare
//!                     two SurrealRBase types with >
//! \nick
//! \param[in] lhs The object on left side of > sign
//! \param[in] rhs The object on right side of > sign
//! \return boolean of the comparison of lhs > rhs
//****************************************************************************80
 template<class LHSType, class RHSType, class realT, int N>
inline bool operator>
(SurrealRBase<LHSType, realT, N> const& lhs,
 SurrealRBase<RHSType, realT, N> const& rhs){
  return(lhs.CastToDerived().GetValue() > rhs.CastToDerived().GetValue());
}

//----------------------------- Real > Surreal ------------------------------
//****************************************************************************80
//! \brief Operator > : This operator declaration declares how to test
//!                     if a real number is > a SurrealRBase type
//! \nick
//! \param[in] lhs The real number on left side of > sign
//! \param[in] rhs The object on right side of > sign
//! \return boolean of the comparison of lhs > rhs
//****************************************************************************80
template<class RHSType>
inline bool operator>
(const typename RHSType::realT_ lhs,
 SurrealRBase<RHSType, typename RHSType::realT_, RHSType::N_> const& rhs){
  return(lhs > rhs.CastToDerived().GetValue());
}

//----------------------------- Surreal > Real ------------------------------
//****************************************************************************80
//! \brief Operator > : This operator declaration declares how to test
//!                     if a SurrealRBase type is > a real number
//! \nick
//! \param[in] lhs The object on left side of > sign
//! \param[in] rhs The real number on right side of > sign
//! \return boolean of the comparison of lhs > rhs
//****************************************************************************80
 template<class LHSType>
 inline bool operator>
 (SurrealRBase<LHSType, typename LHSType::realT_, LHSType::N_> const& lhs,
  const typename LHSType::realT_ rhs){
   return(lhs.CastToDerived().GetValue() > rhs);
 }

//----------------------------- OPERATOR >= ------------------------------------
//----------------------------- Surreal >= Surreal -----------------------------
//****************************************************************************80
//! \brief Operator >= : This operator declaration declares how to compare
//!                     two SurrealRBase types with >=
//! \nick
//! \param[in] lhs The object on left side of >= sign
//! \param[in] rhs The object on right side of >= sign
//! \return boolean of the comparison of lhs >= rhs
//****************************************************************************80
  template<class LHSType, class RHSType, class realT, int N>
  inline bool operator>=
  (SurrealRBase<LHSType, realT, N> const& lhs,
   SurrealRBase<RHSType, realT, N> const& rhs){
  return(lhs.CastToDerived().GetValue() >= rhs.CastToDerived().GetValue());
}

//----------------------------- Real >= Surreal ------------------------------
//****************************************************************************80
//! \brief Operator >= : This operator declaration declares how to test
//!                     if a real number is >= a SurrealRBase type
//! \nick
//! \param[in] lhs The real number on left side of >= sign
//! \param[in] rhs The object on right side of >= sign
//! \return boolean of the comparison of lhs >= rhs
//****************************************************************************80
template<class RHSType>
inline bool operator>=
(const typename RHSType::realT_ lhs,
 SurrealRBase<RHSType, typename RHSType::realT_, RHSType::N_>const& rhs){
  return(lhs >= rhs.CastToDerived().GetValue());
}

//----------------------------- Surreal >= Real ------------------------------
//****************************************************************************80
//! \brief Operator >= : This operator declaration declares how to test
//!                     if a SurrealRBase type is >= a real number
//! \nick
//! \param[in] lhs The object on left side of >= sign
//! \param[in] rhs The real number on right side of >= sign
//! \return boolean of the comparison of lhs >= rhs
//****************************************************************************80
  template<class LHSType>
  inline bool operator>=
  (SurrealRBase<LHSType, typename LHSType::realT_,LHSType::N_> const& lhs,
   const typename LHSType::realT_  rhs){
    return(lhs.CastToDerived().GetValue() >= rhs);
  }
#endif
