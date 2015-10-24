#ifndef SURREALUNARY_H
#define SURREALUNARY_H

#include "Surreal/Forward/SurrealBase.h"
#include <cmath>

/* preprocessor directive function to help define operator classes
   and corresponding operator overloads for simple unary expressions */
#define SURREALS_FUNC1( NAME, OPERATOR_NAME, FUNC, DERIV )                    \
  template<class Expr>                                                        \
  class Surreal ## NAME : public SurrealBase<Surreal ## NAME<Expr>,           \
                                             typename Expr::realT_, Expr::N_> \
  {                                                                           \
  private:                                                                    \
    typedef typename Expr::realT_ realT;                                      \
  public:                                                                     \
    /* constructor definition*/                                               \
    Surreal ## NAME(const Expr& e) : e_(e) {                                  \
      this->value_ = FUNC;                                                    \
    }                                                                         \
  /* Derivative definition*/                                                  \
  inline realT Deriv(const int i) const { return DERIV*e_.Deriv(i); }         \
  private:                                                                    \
  const Expr& e_;                                                             \
  };                                                                          \
                                                                              \
  /* Operator overload definition*/                                           \
  template<class Expr>                                                        \
  inline Surreal ## NAME<Expr>                                                \
  OPERATOR_NAME(const SurrealBase<Expr, typename Expr::realT_, Expr::N_>& z) {\
    return Surreal ## NAME<Expr>( z.CastToDerived() );                        \
  }

//--->define several standard unary functions using the preprocessor directive
//->unary math functions in std
SURREALS_FUNC1( Neg, operator-, -(e_.Value()), -1.0)
SURREALS_FUNC1( Pos, operator+,  (e_.Value()),  1.0)
SURREALS_FUNC1( Log, log, std::log(e_.Value()), 1.0/(e_.Value()))
SURREALS_FUNC1( Log10, log10, std::log10(e_.Value()), 1.0/(e_.Value()*log(10.0)))
SURREALS_FUNC1( Sqrt, sqrt, std::sqrt(e_.Value()),
                1.0/2.0 * std::pow(e_.Value(),-1.0/2.0))
SURREALS_FUNC1( Exp, exp, std::exp(e_.Value()), std::exp(e_.Value()))
SURREALS_FUNC1( Abs, abs, std::abs(e_.Value()), std::copysign(1.0, e_.Value()))
SURREALS_FUNC1( Fabs, fabs, std::fabs(e_.Value()), std::copysign(1.0, e_.Value()))
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

