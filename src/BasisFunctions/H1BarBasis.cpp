/*
 * H1BarBasis.cpp
 *
 *  Created on: Nov 8, 2015
 *      Author: rabbit
 */

#include "BasisFunctions/H1BarBasis.h"
//****************************************************************************80
H1BarBasis::H1BarBasis(){}
//****************************************************************************80
H1BarBasis::~H1BarBasis(){}
//****************************************************************************80
realT H1BarBasis::EvalBasis(const intT& k, const realT& xi)
{
  //---> Local Variables
  realT L0;
  realT L1;

  //---> Function result variable
  realT l;

  //---> Vertex functions
  L0 = (1.0 - xi)*half;
  L1 = (1.0 + xi)*half;

  //---> Based on the specified function  number requested
  switch(k) { // func_select
  case  0: // User is asking for mode 0
    l = L0;
    break;
  case 1: // User is asking for mode 1
    l = L1;
    break;
  default: // User is asking for mode > 1
    /*---> In 1D the order is k, since there a 2 vertex func + 1
      edge bubble */
    l = EdgePoly(k, L0, L1);

    break;
  }// end func_select

  //---> Return result of function
  return(l);

} // End BarH1

//****************************************************************************80
realT H1BarBasis::EvalBasisD(const intT& k, const realT& xi)
{
  //---> Local Variables
  realT L0;
  realT L1;
  realT dL0dxi;
  realT dL1dxi;
  realT dldL0;
  realT dldL1;

  //---> Return variable
  realT dl;

  //---> Vertex functions
  L0 = (1.0 - xi)*half;
  L1 = (1.0 + xi)*half;

  dL0dxi = -half;
  dL1dxi = half;

  //---> Based on the specified function  number requested
  switch(k) { // func_select
  case  0: // User is asking for mode 0
    dl = dL0dxi;
    break;
  case 1: // User is asking for mode 1
    dl = dL1dxi;
    break;
  default: // User is asking for mode > 1

    EdgePolyD(k, L0,  L1, dldL0, dldL1);

    dl = dldL0*dL0dxi + dldL1*dL1dxi;
    //dL0dxi*L1*psi + L0*dL1dxi*psi + half*L1*dpsi;
  }

  //---> Return the derivative to user
  return(dl);
}// end BarH1D

//****************************************************************************80
realT H1BarBasis::EdgePoly( const intT& p, const realT& L0, const realT& L1)
{
  //---> Return to user the formula for an edge polynomial
  return( L0*L1*Polynomials::LobattoKern(p - 2, L1 - L0) );

} // End EdgePoly

//****************************************************************************80
void H1BarBasis::EdgePolyD(const intT& p, const realT& L0, const realT& L1,
         realT& dphidL0, realT& dphidL1)
{
  //---> Compute value of kernel function
  realT psi = Polynomials::LobattoKern(p - 2, L1 - L0);

  //---> Compute value of derivative of kernel function w.r.t (L1 - L0);
  realT dpsi = Polynomials::LobattoKernD(p - 2, L1 - L0);

  //---> Form partial derivative w.r.t. L0
  dphidL0 = L1*psi - L1*L0*dpsi;
  //---> Form partial derivative w.r.t. L1
  dphidL1 = L0*psi + L1*L0*dpsi;
  return;
} // End EdgePolyD


