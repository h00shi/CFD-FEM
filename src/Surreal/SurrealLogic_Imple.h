// -*-c++-*-

//****************************************************************************80
template<class LHSType, class RHSType, int N>
bool operator<(const SurrealBase<LHSType, N>& lhs,
               const SurrealBase<RHSType, N >& rhs){
  return(lhs.CastToDerived().Value() < rhs.CastToDerived().Value());
}

//****************************************************************************80
template<class RHSType, int N>
bool operator<(const typename RHSType::realT_ & lhs,
               const SurrealBase<RHSType, N>& rhs)
{
  return(lhs < rhs.CastToDerived().Value());
}

//****************************************************************************80
template<class LHSType, int N>
bool operator<(const SurrealBase<LHSType, N>& lhs,
               const typename LHSType::realT_ & rhs)
{
  return(lhs.CastToDerived().Value() < rhs);
}
//****************************************************************************80
template<class LHSType, class RHSType, int N>
bool operator<=(const SurrealBase<LHSType, N>& lhs,
                const SurrealBase<RHSType, N>& rhs)
{
  return(lhs.CastToDerived().Value() <= rhs.CastToDerived().Value());
}
//****************************************************************************80
template<class RHSType, int N>
bool operator<=(const typename RHSType::realT_ & lhs,
                const SurrealBase<RHSType, N>& rhs)
{
  return(lhs <= rhs.CastToDerived().Value());
}
//****************************************************************************80
template<class LHSType, int N>
bool operator<=(const SurrealBase<LHSType, N>& lhs,
                const typename LHSType::realT_ & rhs)
{
  return(lhs.CastToDerived().Value() <= rhs);
}
//****************************************************************************80
template<class LHSType, class RHSType, int N>
bool operator>(const SurrealBase<LHSType, N>& lhs,
               const SurrealBase<RHSType, N>& rhs)
{
  return(lhs.CastToDerived().Value() > rhs.CastToDerived().Value());
}
//****************************************************************************80
template<class RHSType, int N>
bool operator>(const typename RHSType::realT_ & lhs,
               const SurrealBase<RHSType, N>& rhs)
{
  return(lhs > rhs.CastToDerived().Value());
}
//****************************************************************************80
template<class LHSType, int N>
bool operator>(const SurrealBase<LHSType, N>& lhs,
               const typename LHSType::realT_ & rhs)
{
  return(lhs.CastToDerived().Value() > rhs);
}
//****************************************************************************80
template<class LHSType, class RHSType, int N>
bool operator>=(const SurrealBase<LHSType, N>& lhs,
                const SurrealBase<RHSType, N>& rhs)
{
  return(lhs.CastToDerived().Value() >= rhs.CastToDerived().Value());
}
//****************************************************************************80
template<class RHSType, int N>
bool operator>=(const typename RHSType::realT_ & lhs,
                const SurrealBase<RHSType, N>& rhs)
{
  return(lhs >= rhs.CastToDerived().Value());
}
//****************************************************************************80
template<class LHSType, int N>
bool operator>=(const SurrealBase<LHSType, N>& lhs,
               const typename LHSType::realT_ & rhs)
{
  return(lhs.CastToDerived().Value() >= rhs);
}
// //****************************************************************************80
// template<class LHSType, class RHSType, int N>
// bool operator==(const SurrealBase<LHSType, N>& lhs,
//                 const SurrealBase<RHSType, N>& rhs)
// {
//   return(lhs.CastToDerived().Value() == rhs.CastToDerived().Value());
// }
// //****************************************************************************80
// template<class RHSType, int N>
// bool operator==(const typename RHSType::realT_ & lhs,
//                 const SurrealBase<RHSType, N>& rhs)
// {
//   return(lhs == rhs.CastToDerived().Value());
// }
// //****************************************************************************80
// template<class LHSType, int N>
// bool operator==(const SurrealBase<LHSType,N>& lhs,
//                const typename LHSType::realT_ & rhs)
// {
//   return(lhs.CastToDerived().Value() == rhs);
// }
