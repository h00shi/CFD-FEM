//-*-c++-*-
//****************************************************************************80
template<class LHSType, class RHSType, int N>
inline SurrealAdd<LHSType, RHSType, N>
operator+(const SurrealBase< LHSType, N>& lhs,
          const SurrealBase< RHSType, N>& rhs)
{
  /*---> Assuming that all possible Types of LHSType and RHSType are derived
    from SurrealBase<Derived>. Then create an instance of SurrealAdd with
    an argument of type LHSType by casting from SurrealBase<LHSType> to
    LHSType and an argument of type RHSType by casting from
    SurrealBase<RHSType> to RHSType */
  return(SurrealAdd<LHSType,RHSType, N>(lhs.CastToDerived(),
                                        rhs.CastToDerived()));
}
//****************************************************************************80
template<class RHSType, int N>
inline SurrealAdd<typename RHSType::realT_, RHSType, N>
operator+
(const typename RHSType::realT_& lhs, const SurrealBase<RHSType, N>& rhs)
{
  return(SurrealAdd<typename RHSType::realT_,RHSType, N>
         (lhs, rhs.CastToDerived()));
}
//****************************************************************************80
template<class LHSType, int N>
inline SurrealAdd<typename LHSType::realT_, LHSType, N>
operator+
(const SurrealBase<LHSType, N>& lhs, const typename LHSType::realT_& rhs)
{
  return(SurrealAdd<typename LHSType::realT_,LHSType, N>
         (rhs, lhs.CastToDerived()));
}
