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
//! \brief jacobi : Computes the value of a jacobia polynomial.  
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
  template<typename intT, typename realT>
  realT jacobi(const intT& n, const intT& alpha, const intT& beta, 
	       const realT& xi) 
  {
    //---> Local Variables
    realT p0;
    realT p1;
    realT an1;
    realT an2;
    realT an3;
    realT an4;
    realT ar;
    realT br;
    realT im1;
    
    //---> Function result variable 
    realT p;
    
    /*---> Cast the variables alpha, beta to doubles to ensure proper precision 
      of computation */ 
    ar = (realT)alpha;
    br = (realT)beta;

    //---> Initialize recurrence formula. 
    p0 = 1.0;
    p1 = half*(ar - br + (ar + br + 2.0)*xi);
    
    if ( n == 0 ) { // order_check
      //---> If we need 0th polynomial then just make p = p0
      p = p0;
    }
    else if( n == 1 ) {
      //---> If we need the 1st polynomila then just make p = p1
      p = p1;
    }
    else {
      //---> Otherwise evaluate the recurrence relation
      
      
      for(intT i = 2; i <= n; i++) { //poly_loop 
	realT N = (realT)i - 1.0;
	an1 = 2.0*(N + 1.0)*(N + ar + br + 1.0)*(2.0*N + ar + br);
	an2 = (2.0*N + ar + br + 1.0)*(ar*ar - br*br);
	an3 = (2.0*N + ar + br)*(2.0*N + ar + br + 1.0)*(2.0*N + 
							 ar + br + 2.0);
	an4 = 2.0*(N + ar)*(N + br)*(2.0*N + ar + br + 2.0);
	
	p = 1.0/an1*((an2 + an3*xi)*p1 - an4*p0);
	p0 = p1;
	p1 = p;

      } // End poly_loop
      
    } // End order_check
    
    //---> Return result of fuction 
    return(p);

  } // End jacobi
  
//****************************************************************************80
//!
//! \brief jacobiD : The derivative of the jacobi polynomial w.r.t \f$ \xi \f$
//! \details Computes \f$ \frac{d P_{n}^{\alpha, \beta}(\xi)}{d \xi} \f$.  
//! The derivative of a jacobi polynomial is defined by: 
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
 template<typename intT, typename realT>
 realT jacobiD(const intT& n, const intT& alpha, const intT& beta, 
  	       const realT& xi) 
  {
    //---> Local Variables
    intT bn1;
    intT bn2;
    realT pnm1;
    
    //---> Function result variable
    realT dp;
    

    if( n == 0 ) { // order_check
      dp = 0.0;
    }
    else {
      bn1 = alpha + 1;
      bn2 = beta + 1;
      pnm1 = jacobi(n - 1, bn1, bn2, xi);
      
      dp = half*( (realT)alpha + (realT)beta + (realT)n + 1.0)*pnm1;
    } // end order_check
    
    //---> Resturn result of fuction 
    return(dp);
    
  } // end jacobiD


//****************************************************************************80
//!
//! \brief lobatto_kern : Kernel function for lobatto functions. 
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
  template<typename intT, typename realT> 
  realT lobatto_kern (const intT& km2, const realT& xi)
  {
    //---> Local Variables
    intT k = km2 + 2; /*---> Coefficients of the kernel function are 
			functions of k */
    realT sigma = 1.0/sqrt(2.0/(2.0*(realT)k - 1.0)); //---> Coefficient
    
    //---> Function return variable
    realT psi;

    //---> Define kernel function
    psi = -2.0*sigma/((realT)k - 1.0)*jacobi<intT,realT>(km2, 1, 1, xi);
    
    //---> Return psi to user, which kernel evaluated at xi
    return(psi);
    
  } // End lobatto_kern

//****************************************************************************80
//!
//! \brief lobatto_kernD : Derivative of lobatto_kern w.r.t. \f$ \xi \f$
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] km2 \f$ k-2 \f$ Degree of kernel function
//! \param[in] xi \f$ \xi \f$
//! \return psi \f$ \psi_{k-2}(\xi) \f$ i.e. psi_{km2} (xi).  
//****************************************************************************80
  template<typename intT, typename realT> 
  realT lobatto_kernD(const intT& km2, const realT& xi)
  {
    //---> Local Variables
    intT k = km2 + 2; /*---> Coefficients of the kernel function are 
			functions of k */
    realT sigma = 1.0/sqrt(2.0/(2.0*(realT)k - 1.0)); //---> Coefficient
   
    //---> Function return variable
    realT dpsi;
       
    //---> Define kernel function derivative
    dpsi = -2.0*sigma/((realT)k - 1.0)*jacobiD<intT, realT>(km2, 1, 1, xi);
    
    //---> Return dpsi to user, which is derivative of kernel at xi;
    return(dpsi);
  } // End lobatto_kernD

} // End namespace Polynomials

#endif
