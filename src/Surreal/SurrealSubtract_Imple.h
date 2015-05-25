//-*-c++-*-
//****************************************************************************80
template<class LHSType, class RHSType, int N>
inline SurrealSubtract<LHSType, RHSType, N>
operator-(const SurrealBase<LHSType, N>& lhs,
          const SurrealBase<RHSType, N>& rhs)
{
  return(SurrealSubtract<LHSType,RHSType, N>
         (lhs.CastToDerived(), rhs.CastToDerived()));
}
//****************************************************************************80
template<class RHSType, int N>
inline SurrealSubtract<typename RHSType::realT_, RHSType, N>
operator-
(const typename RHSType::realT_& lhs, const SurrealBase<RHSType, N>& rhs)
{
  return(SurrealSubtract<typename RHSType::realT_,RHSType, N>
         (lhs, rhs.CastToDerived()));
}
//****************************************************************************80
template<class LHSType, int N>
inline SurrealSubtract<LHSType, typename LHSType::realT_, N>
operator-
(const SurrealBase<LHSType, N>& lhs, const typename LHSType::realT_& rhs)
{
  return(SurrealSubtract<LHSType, typename LHSType::realT_, N>
         (lhs.CastToDerived(), rhs));
}
