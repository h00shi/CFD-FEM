#ifndef SURREALLOGIC_H
#define SURREALLOGIC_H

#include "Surreal/Forward/SurrealBase.h"
//----------------------------- LOGICAL OPERATORS ------------------------------
//----------------------------- OPERATOR < -------------------------------------
//----------------------------- Surreal < Surreal ------------------------------
//****************************************************************************80
//! \brief Operator < : This operator declaration declares how to compare
//!                     two SurrealBase types with <
//! \nick
//! \param[in] lhs The object on left side of < sign
//! \param[in] rhs The object on right side of < sign
//! \return boolean of the comparison of lhs < rhs
//****************************************************************************80
template<class LHSType, class RHSType, class realT, int N>
inline bool operator<(const SurrealBase<LHSType, realT, N>& lhs,
                      const SurrealBase<RHSType, realT, N>& rhs){
  return(lhs.CastToDerived().Value() < rhs.CastToDerived().Value());
}

//----------------------------- Real < Surreal -------------------------------
//****************************************************************************80
//! \brief Operator < : This operator declaration declares how to test
//!                     if a real number is less than a SurrealBase type
//! \nick
//! \param[in] lhs The real number on left side of < sign
//! \param[in] rhs The object on right side of < sign
//! \return boolean of the comparison of lhs < rhs
//****************************************************************************80
template<class RHSType>
inline bool operator<
(const typename RHSType::realT_ lhs,
 const SurrealBase<RHSType,typename RHSType::realT_, RHSType::N_>& rhs){
  return(lhs < rhs.CastToDerived().Value());
}

//----------------------------- Surreal < Real -------------------------------
//****************************************************************************80
//! \brief Operator < : This operator declaration declares how to test
//!                     if a SurrealBase type is less thana real number
//! \nick
//! \param[in] lhs The object on left side of < sign
//! \param[in] rhs The real number on right side of < sign
//! \return boolean of the comparison of lhs < rhs
//****************************************************************************80
template<class LHSType>
inline bool operator<
(const SurrealBase<LHSType, typename LHSType::realT_, LHSType::N_>& lhs,
 const typename LHSType::realT_ rhs){
  return(lhs.CastToDerived().Value() < rhs);
}

//----------------------------- OPERATOR <= ------------------------------------
//----------------------------- Surreal <= Surreal -----------------------------
//****************************************************************************80
//! \brief Operator <= : This operator declaration declares how to compare
//!                     two SurrealBase types with <=
//! \nick
//! \param[in] lhs The object on left side of <= sign
//! \param[in] rhs The object on right side of <= sign
//! \return boolean of the comparison of lhs <= rhs
//****************************************************************************80
template<class LHSType, class RHSType, class realT, int N>
inline bool operator<=
(const SurrealBase<LHSType, realT, N>& lhs,
 const SurrealBase<RHSType, realT, N>& rhs){
  return(lhs.CastToDerived().Value() <= rhs.CastToDerived().Value());
}

//----------------------------- Real <= Surreal ------------------------------
//****************************************************************************80
//! \brief Operator <= : This operator declaration declares how to test
//!                     if a real number is <= a SurrealBase type
//! \nick
//! \param[in] lhs The real number on left side of <= sign
//! \param[in] rhs The object on right side of <= sign
//! \return boolean of the comparison of lhs <= rhs
//****************************************************************************80
template<class RHSType>
inline bool operator<=
(const typename RHSType::realT_ lhs,
 const SurrealBase<RHSType,typename RHSType::realT_, RHSType::N_>& rhs){
  return(lhs <= rhs.CastToDerived().Value());
}

//----------------------------- Surreal <= Real ------------------------------
//****************************************************************************80
//! \brief Operator <= : This operator declaration declares how to test
//!                     if a SurrealBase type is <= a real number
//! \nick
//! \param[in] lhs The object on left side of <= sign
//! \param[in] rhs The real number on right side of <= sign
//! \return boolean of the comparison of lhs <= rhs
//****************************************************************************80
template<class LHSType>
inline bool operator<=
(const SurrealBase<LHSType, typename LHSType::realT_, LHSType::N_>& lhs,
 const typename LHSType::realT_ rhs){
  return(lhs.CastToDerived().Value() <= rhs);
}

//----------------------------- OPERATOR > ------------------------------------
//----------------------------- Surreal > Surreal -----------------------------
//****************************************************************************80
//! \brief Operator > : This operator declaration declares how to compare
//!                     two SurrealBase types with >
//! \nick
//! \param[in] lhs The object on left side of > sign
//! \param[in] rhs The object on right side of > sign
//! \return boolean of the comparison of lhs > rhs
//****************************************************************************80
template<class LHSType, class RHSType, class realT, int N>
inline bool operator>
(const SurrealBase<LHSType, realT, N>& lhs,
 const SurrealBase<RHSType, realT, N>& rhs){
  return(lhs.CastToDerived().Value() > rhs.CastToDerived().Value());
}

//----------------------------- Real > Surreal ------------------------------
//****************************************************************************80
//! \brief Operator > : This operator declaration declares how to test
//!                     if a real number is > a SurrealBase type
//! \nick
//! \param[in] lhs The real number on left side of > sign
//! \param[in] rhs The object on right side of > sign
//! \return boolean of the comparison of lhs > rhs
//****************************************************************************80
template<class RHSType>
inline bool operator>
(const typename RHSType::realT_ lhs,
 const SurrealBase<RHSType,typename RHSType::realT_, RHSType::N_>& rhs){
  return(lhs > rhs.CastToDerived().Value());
}

//----------------------------- Surreal > Real ------------------------------
//****************************************************************************80
//! \brief Operator > : This operator declaration declares how to test
//!                     if a SurrealBase type is > a real number
//! \nick
//! \param[in] lhs The object on left side of > sign
//! \param[in] rhs The real number on right side of > sign
//! \return boolean of the comparison of lhs > rhs
//****************************************************************************80
template<class LHSType>
inline bool operator>
(const SurrealBase<LHSType, typename LHSType::realT_, LHSType::N_>& lhs,
 const typename LHSType::realT_ rhs){
  return(lhs.CastToDerived().Value() > rhs);
}

//----------------------------- OPERATOR >= ------------------------------------
//----------------------------- Surreal >= Surreal -----------------------------
//****************************************************************************80
//! \brief Operator >= : This operator declaration declares how to compare
//!                     two SurrealBase types with >=
//! \nick
//! \param[in] lhs The object on left side of >= sign
//! \param[in] rhs The object on right side of >= sign
//! \return boolean of the comparison of lhs >= rhs
//****************************************************************************80
template<class LHSType, class RHSType, class realT, int N>
inline bool operator>=
(const SurrealBase<LHSType, realT, N>& lhs,
 const SurrealBase<RHSType, realT, N>& rhs){
  return(lhs.CastToDerived().Value() >= rhs.CastToDerived().Value());
}

//----------------------------- Real >= Surreal ------------------------------
//****************************************************************************80
//! \brief Operator >= : This operator declaration declares how to test
//!                     if a real number is >= a SurrealBase type
//! \nick
//! \param[in] lhs The real number on left side of >= sign
//! \param[in] rhs The object on right side of >= sign
//! \return boolean of the comparison of lhs >= rhs
//****************************************************************************80
template<class RHSType>
inline bool operator>=
(const typename RHSType::realT_ lhs,
 const SurrealBase<RHSType,typename RHSType::realT_, RHSType::N_>& rhs){
  return(lhs >= rhs.CastToDerived().Value());
}

//----------------------------- Surreal >= Real ------------------------------
//****************************************************************************80
//! \brief Operator >= : This operator declaration declares how to test
//!                     if a SurrealBase type is >= a real number
//! \nick
//! \param[in] lhs The object on left side of >= sign
//! \param[in] rhs The real number on right side of >= sign
//! \return boolean of the comparison of lhs >= rhs
//****************************************************************************80
template<class LHSType>
inline bool operator>=
(const SurrealBase<LHSType, typename LHSType::realT_, LHSType::N_>& lhs,
 const typename LHSType::realT_ rhs){
  return(lhs.CastToDerived().Value() >= rhs);
}

#endif
