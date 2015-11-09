/*
 * BarH1Basis.h
 *
 *  Created on: Nov 8, 2015
 *      Author: rabbit
 */

#ifndef H1BARBASIS_H_
#define H1BARBASIS_H_
#include "my_incl.h"
#include "consts.h"
#include "BasisFunctions/Polynomials.h"
//****************************************************************************80
//! \brief A class the defines the interface to evaluate H1-hierarchical basis
//!        functions.
//! \nick
//!
//****************************************************************************80
class H1BarBasis
{
public:
//****************************************************************************80
//!
//! \brief H1BarBasis : The default constructor for the class H1BarBasis
//! \details
//! \nick
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//!
//****************************************************************************80
  H1BarBasis();
//****************************************************************************80
//!
//! \brief ~H1BarBasis : The destructor for the class H1BarBasis
//! \details
//! \nick
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//!
//****************************************************************************80
  ~H1BarBasis();
//****************************************************************************80
//!
//! \brief EvalBasis : A public interface to evaluate a specified polynomial
//! \details
//! \nick
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] k The polyomial you want to evaluate
//! \param[in] xi The standard element coordinate \f$\xi\f$ where
//! you want the basis function value.
//! \return phi Value of basis functio n \f$ \phi_{n} \left( \xi\right) \f$
//****************************************************************************80
  realT EvalBasis(const intT& k, const realT& xi);

//****************************************************************************80
//!
//! \brief EvalBasisD : A public interface to evaluate a the derivate of
//!                      a specified polynomail.
//! \details
//! \nick
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] k The polyomial you want to evaluate
//! \param[in] xi The standard element coordinate \f$\xi\f$ where you want the
//!  basis function derivate value.
//! \return phi Value of basis function
//! \f$ \frac{ d \phi_{n}}{d \xi} \left( \xi\right) \f$
//****************************************************************************80
  realT EvalBasisD(const intT& k, const realT& xi);

protected:
//****************************************************************************80
//!
//! \brief EdgePoly : Compute the value of a lobatto polynomial on an edge
//!                   given the affine coordinates L0, L1.
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[in] p The degree of polynomial you want: MUST BE > 2.
//! \param[in] L0 The first affine coordinate or "left" one
//! \param[in] L1 The second affine coorindat or "right" one
//! \return phi The edge polynomial value at evaluated at L1 - L0;
//****************************************************************************80
  realT EdgePoly( const intT& p, const realT& L0, const realT& L1);

//****************************************************************************80
//!
//! \brief edgepolyD : Compute the derivative of the edge polynomial w.r.t
//!                    affine coordintes L0, L1
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[in] p The degree of polynomila you : MUST BE > 2.
//! \param[in] L0 The first affine coordinate or "left" one
//! \param[in] L1 The second affine coorindat or "right" one
//! \param[out] dphidL0 \f$ \frac{ \partial phi}{\partial L0} \f$
//! \param[out] dphidL1 \f$ \frac{ \partial phi}{\partial L1} \f$
//****************************************************************************80
  void EdgePolyD(const intT& p, const realT& L0, const realT& L1,
     realT& dphidL0, realT& dphidL1);
};

#endif /* H1BARBASIS_H_ */
