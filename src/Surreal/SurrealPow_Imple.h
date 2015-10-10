// -*-c++-*-
//****************************************************************************80
template<class LHSType, class RHSType, int N>
inline SurrealPow<LHSType, RHSType, N>
pow(const SurrealBase<LHSType,N>& lhs, const SurrealBase<RHSType, N>& rhs)
{
  return(SurrealPow<LHSType,RHSType,N>
         (lhs.CastToDerived(), rhs.CastToDerived()));
}
//****************************************************************************80
template<class RHSType, int N>
inline SurrealPow<typename RHSType::realT_, RHSType, N>
pow(const typename RHSType::realT_& lhs, const SurrealBase<RHSType, N>& rhs)
{
  return(SurrealPow<typename RHSType::realT_,RHSType, N>
         (lhs, rhs.CastToDerived()));
}
//****************************************************************************80
template<class LHSType, int N>
inline SurrealPow<LHSType, typename LHSType::realT_, N>
pow(const SurrealBase<LHSType, N>& lhs, const typename LHSType::realT_ & rhs)
{
  return(SurrealPow<LHSType, typename LHSType::realT_, N>
         (lhs.CastToDerived(), rhs));
}
