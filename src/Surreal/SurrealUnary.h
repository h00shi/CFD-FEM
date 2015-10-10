// -*-c++-*-
#ifndef SURREALUNARY_H
#define SURREALUNARY_H

#include "Surreal/SurrealBase.h"

/* preprocessor directive function to help define operator classes
   and corresponding operator overloads for simple unary expressions */
#define SURREALS_FUNC1( NAME, OPERATOR_NAME, FUNC, DERIV )                \
  template<class Expr, int N>                                             \
  class Surreal ## NAME : public SurrealBase<Surreal ## NAME<Expr, N>, N> \
  {                                                                       \
  public:                                                                 \
    typedef typename Expr::realT_ realT;                                  \
    typedef realT realT_;                                                 \
                                                                          \
    /* constructor definition*/                                           \
    Surreal ## NAME(const Expr& e) : e_(e) {}                             \
                                                                          \
  /* Value definition*/                                                   \
  inline realT Value() const { return FUNC; }                             \
                                                                          \
  /* Derivative definition*/                                              \
  inline realT Deriv(const int& i) const { return DERIV*e_.Deriv(i); }    \
  private:                                                                \
  const Expr& e_;                                                         \
  };                                                                      \
                                                                          \
  /* Operator overload definition*/                                       \
  template<class Expr, int N>                                             \
  inline Surreal ## NAME<Expr, N>                                         \
  OPERATOR_NAME(const SurrealBase<Expr, N>& z);                           \
                                                                          \
  /* Operator overload definition*/                                       \
  template<class Expr, int N>                                             \
  inline Surreal ## NAME<Expr, N>                                         \
  OPERATOR_NAME(const SurrealBase<Expr, N>& z) {                          \
    return Surreal ## NAME<Expr, N>( z.CastToDerived() );                 \
  }

//--->define several standard unary functions using the preprocessor directive
//->unary math functions in std
SURREALS_FUNC1( Neg, operator-, -(e_.Value()), -1.0)
SURREALS_FUNC1( Pos, operator+,  (e_.Value()),  1.0)
SURREALS_FUNC1( Log, log, std::log(e_.Value()), 1.0/(e_.Value()))
SURREALS_FUNC1( Sqrt, sqrt, std::sqrt(e_.Value()),
                1.0/2.0 * std::pow(e_.Value(),-1.0/2.0))
SURREALS_FUNC1( Exp, exp, std::exp(e_.Value()), std::exp(e_.Value()))
SURREALS_FUNC1( Abs, abs, std::abs(e_.Value()), std::copysign(1.0, e_.Value()))

//->trig functions in std
SURREALS_FUNC1( Cos, cos, std::cos(e_.Value()), -std::sin(e_.Value()) )
SURREALS_FUNC1( Sin, sin, std::sin(e_.Value()),  std::cos(e_.Value()) )
SURREALS_FUNC1( Tan, tan, std::tan(e_.Value()),
                1.0/(std::cos(e_.Value())*std::cos(e_.Value())) )
SURREALS_FUNC1( Acos, acos, std::acos(e_.Value()),
                -1.0/std::sqrt(1 - e_.Value()*e_.Value()) )
SURREALS_FUNC1( Asin, asin, std::asin(e_.Value()),
                1.0/std::sqrt(1 - e_.Value()*e_.Value()) )
SURREALS_FUNC1( Atan, atan, std::atan(e_.Value()),
                1.0/(1.0 + e_.Value()*e_.Value()) )

#endif

