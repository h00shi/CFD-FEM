
// -*-c++-*-
#include "my_incl.h"
#include "Polynomials.h"

//****************************************************************************80
realT Polynomials::Jacobi(const intT& n, const intT& alpha, const intT& beta, 
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

} // End Jacobi

//****************************************************************************80
realT Polynomials::JacobiD(const intT& n, const intT& alpha, const intT& beta, 
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
    pnm1 = Jacobi(n - 1, bn1, bn2, xi);
      
    dp = half*( (realT)alpha + (realT)beta + (realT)n + 1.0)*pnm1;
  } // end order_check
    
    //---> Resturn result of fuction 
  return(dp);
    
} // end JacobiD
//****************************************************************************80
realT Polynomials::LobattoKern (const intT& km2, const realT& xi)
{
  //---> Local Variables
  intT k = km2 + 2; /*---> Coefficients of the kernel function are 
		      functions of k */
  realT sigma = 1.0/sqrt(2.0/(2.0*(realT)k - 1.0)); //---> Coefficient
    
  //---> Function return variable
  realT psi;

  //---> Define kernel function
  psi = -2.0*sigma/((realT)k - 1.0)*Jacobi(km2, 1, 1, xi);
    
  //---> Return psi to user, which kernel evaluated at xi
  return(psi);
    
} // End LobattoKern

//****************************************************************************80
realT Polynomials::LobattoKernD(const intT& km2, const realT& xi)
{
  //---> Local Variables
  intT k = km2 + 2; /*---> Coefficients of the kernel function are 
		      functions of k */
  realT sigma = 1.0/sqrt(2.0/(2.0*(realT)k - 1.0)); //---> Coefficient
   
  //---> Function return variable
  realT dpsi;
       
  //---> Define kernel function derivative
  dpsi = -2.0*sigma/((realT)k - 1.0)*JacobiD(km2, 1, 1, xi);
    
  //---> Return dpsi to user, which is derivative of kernel at xi;
  return(dpsi);
} // End LobattoKernD
