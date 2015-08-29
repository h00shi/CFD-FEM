// -*-c++-*-
#ifndef POLYNOMIALS
#define POLYNOMIALS
//****************************************************************************80
//! \file Polynomials.h 
//! \namespace Polynomials Polynomials.h
//! \brief This is the header file defining the namespace Polynomials. 
//! \details The functions compute polynomials.  Things like Jacobi polyomials 
//! \nick
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//****************************************************************************80

#include "my_incl.h"
#include "consts.h"

namespace Polynomials {
//****************************************************************************80
//! \brief Jacobi : Computes the value of a Jacobia polynomial.  
//! Which polynomial and where is specified by the inputs 
//! \details Computes formula: \f$ P_{n}^{\alpha, \beta} ( \xi ) \f$ 
//! using a reccurence formula: 
//! \f$ P_{n}^{\alpha, \beta} ( \xi ) = \frac{ a_{n}^{2} + a_{n}^{3}) P_{n-1}^{\alpha, \beta} (\xi) - a_{n}^{4} P_{n-2}^{\alpha, \beta} (\xi) } {a_{n}^{1} } \f$
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \tparam intT Integer data type
//! \tparam realT Real data type
//! \param[in] n The order of the Jacobi polynomial
//! \param[in] alpha The parameter \f$ \alpha \f$ in the Jacobi polynomial
//! \param[in] beta The parameter \f$ \beta \f$ in the Jacobi polynomial
//! \param[in] xi The \f$ \xi \f$ coordinate \f$ \xi \in [-1,1] \f$ 
//! \return p The result of the function 
//****************************************************************************80
  realT Jacobi(const intT& n, const intT& alpha, const intT& beta, 
	       const realT& xi); 
  
//****************************************************************************80
//!
//! \brief JacobiD : The derivative of the Jacobi polynomial w.r.t \f$ \xi \f$
//! \details Computes \f$ \frac{d P_{n}^{\alpha, \beta}(\xi)}{d \xi} \f$.  
//! The derivative of a Jacobi polynomial is defined by: 
//! \f$ \frac{1}{2}(\alpha + \beta + n + 1)P_{n-1}^{\alpha + 1, \beta + 1}(\xi) \f$
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \tparam intT Integer data type
//! \tparam realT Real data type
//! \param[in] n The order of the Jacobi polynomial
//! \param[in] alpha The parameter \f$ \alpha \f$ in the Jacobi polynomial
//! \param[in] beta The parameter \f$ \beta \f$ in the Jacobi polynomial
//! \param[in] xi The \f$ \xi \f$ coordinate \f$ \xi \in [-1,1] \f$ 
//! \return dp \formula{\frac{d p}{d \xi} }
//****************************************************************************80
  realT JacobiD(const intT& n, const intT& alpha, const intT& beta, 
		const realT& xi);
 
//****************************************************************************80
//!
//! \brief LobattoKern : Kernel function for lobatto functions. 
//!  NOT the actual lobatto function but the KERNEL.  
//! \details Computes the kernel function for a Lobatto polynomial.  
//!  Say you want a Lobatto polynomial of k > 2.  
//!  \f$ l_{k>2} = l_{0}l_{1}*\psi_{k-2} \f$. This function returns the 
//!  \f$ \psi_{k-2}(\xi) \f$ part.
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] km2 \f$ k-2 \f$ Degree of kernel function
//! \param[in] xi \f$ \xi \f$
//! \return psi \f$ \psi_{k-2}(\xi) \f$ i.e. psi_{km2} (xi).  
//****************************************************************************80
  realT LobattoKern (const intT& km2, const realT& xi);
 
//****************************************************************************80
//!
//! \brief LobattoKernD : Derivative of LobattoKern w.r.t. \f$ \xi \f$
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] km2 \f$ k-2 \f$ Degree of kernel function
//! \param[in] xi \f$ \xi \f$
//! \return psi \f$ \psi_{k-2}(\xi) \f$ i.e. psi_{km2} (xi).  
//****************************************************************************80
  realT LobattoKernD(const intT& km2, const realT& xi);
  
} // End namespace Polynomials

#endif
