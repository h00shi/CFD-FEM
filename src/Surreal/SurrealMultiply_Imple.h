//-*-c++-*-
//****************************************************************************80
template<class LHSType, class RHSType, int N>
inline SurrealMultiply<LHSType, RHSType, N>
operator*(const SurrealBase<LHSType, N>& lhs,
          const SurrealBase<RHSType, N>& rhs)
{
  return(SurrealMultiply<LHSType,RHSType, N>
         (lhs.CastToDerived(), rhs.CastToDerived()));
}
//****************************************************************************80
template<class RHSType, int N>
inline SurrealMultiply<typename RHSType::realT_, RHSType, N>
operator*
(const typename RHSType::realT_& lhs, const SurrealBase<RHSType, N>& rhs)
{
  return(SurrealMultiply<typename RHSType::realT_,RHSType,N>
         (lhs, rhs.CastToDerived()));
}
//****************************************************************************80
template<class LHSType, int N>
inline SurrealMultiply<typename LHSType::realT_, LHSType, N>
operator*
(const SurrealBase<LHSType, N>& lhs, const typename LHSType::realT_& rhs)
{
  return(SurrealMultiply<typename LHSType::realT_,LHSType, N>
         (rhs, lhs.CastToDerived()));
}
