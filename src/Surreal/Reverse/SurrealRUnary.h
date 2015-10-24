#ifndef SURREALRUNARY_H
#define SURREALRUNARY_H
#include <cmath>
#include "Surreal/Reverse/SurrealRBase.h"

/* preprocessor directive function to help define operator classes
   and corresponding operator overloads for simple unary expressions */
#define SURREALR_UNARY( NAME, OPERATOR_NAME, FUNC, DERIV )              \
  template<class Expr>                                                  \
  class SurrealR ## NAME : public SurrealRBase<SurrealR ## NAME<Expr>,  \
                                               typename Expr::realT_, Expr::N_> \
  {                                                                     \
  private:                                                              \
    typedef typename Expr::realT_ realT;                                \
    Expr const & e_;                                                    \
  public:                                                               \
    /* constructor definition*/                                         \
    SurrealR ## NAME(Expr const& e) : e_(e) {                           \
      this->value_ = FUNC;                                              \
    }                                                                   \
                                                                        \
    void Diff(realT adjoint,                                            \
              realT (&deriv)[Expr::N_], int (&deriv_ix)[Expr::N_],      \
              int (&ix)[Expr::N_], unsigned int & count) const{         \
      e_.Diff(adjoint*DERIV, deriv, deriv_ix, ix, count);               \
    }                                                                   \
  };                                                                    \
                                                                        \
  /* Operator overload definition*/                                     \
  template<class Expr>                                                  \
  inline SurrealR ## NAME<Expr>                                         \
  OPERATOR_NAME(SurrealRBase<Expr, typename Expr::realT_, Expr::N_>const& z) { \
    return SurrealR ## NAME<Expr>( z.CastToDerived() );                 \
  }

//--->define several standard unary functions using the preprocessor directive
//->unary math functions in std
SURREALR_UNARY( Pos, operator+,  (e_.GetValue()),  1.0)
SURREALR_UNARY( Neg, operator-, -(e_.GetValue()), -1.0)

SURREALR_UNARY( Log, log, std::log(e_.GetValue()), 1.0/(e_.GetValue()))
SURREALR_UNARY( Log10, log10, std::log10(e_.GetValue()), 1.0/(e_.GetValue()*log(10.0)))

SURREALR_UNARY( Sqrt, sqrt, std::sqrt(e_.GetValue()),
                1.0/2.0 * std::pow(e_.GetValue(),-1.0/2.0))
SURREALR_UNARY( Exp, exp, std::exp(e_.GetValue()), std::exp(e_.GetValue()))
SURREALR_UNARY( Abs, abs, std::abs(e_.GetValue()), std::copysign(1.0, e_.GetValue()))
SURREALR_UNARY( Fabs, fabs, std::fabs(e_.GetValue()), std::copysign(1.0, e_.GetValue()))
//->trig functions in std
SURREALR_UNARY( Cos, cos, std::cos(e_.GetValue()), -std::sin(e_.GetValue()) )
SURREALR_UNARY( Sin, sin, std::sin(e_.GetValue()),  std::cos(e_.GetValue()) )
SURREALR_UNARY( Tan, tan, std::tan(e_.GetValue()),
                1.0/(std::cos(e_.GetValue())*std::cos(e_.GetValue())) )
SURREALR_UNARY( Acos, acos, std::acos(e_.GetValue()),
                -1.0/std::sqrt(1 - e_.GetValue()*e_.GetValue()) )
SURREALR_UNARY( Asin, asin, std::asin(e_.GetValue()),
                1.0/std::sqrt(1 - e_.GetValue()*e_.GetValue()) )
SURREALR_UNARY( Atan, atan, std::atan(e_.GetValue()),
                1.0/(1.0 + e_.GetValue()*e_.GetValue()) )

#endif

