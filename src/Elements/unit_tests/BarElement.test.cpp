#include "BarElement.h"
#include "Array1D.h"
#include "gtest/gtest.h"
#include "my_incl.h"
#include "Polynomials.h"
#include "Element.h"

void quadpoly (const Array1D<realT>&, Array1D<realT>& );
void p3poly (const Array1D<realT>&, Array1D<realT>& );
void p4poly (const Array1D<realT>&, Array1D<realT>& );

//****************************************************************************80
//! \file barelement_test.cpp
//! \brief Collection of google tests for namespace Polynomials
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! 
//****************************************************************************80
const realT tol = 9.9e-13;
TEST(BarElement, eval_basis) {
  realT phi = -9.9e99;
  BarElement sbar;

  //-------------------------------- k = 0 -------------------------------------
  phi = sbar.EvalBasis(0, 0.0);
  EXPECT_NEAR(.5, phi, tol);
  
  phi = sbar.EvalBasis(0, -1.0);
  EXPECT_NEAR(1.0, phi, tol);
  
  phi = sbar.EvalBasis(0, 1.0);
  EXPECT_NEAR(0.0, phi, tol);

  phi = sbar.EvalBasis(0, .5526);
  EXPECT_NEAR(.2237, phi, tol);
  
  //-------------------------------- k = 1 -------------------------------------
  phi = sbar.EvalBasis(1, 0.0);
  EXPECT_NEAR(.5, phi, tol);
  
  phi = sbar.EvalBasis(1, -1.0);
  EXPECT_NEAR(0.0, phi, tol);
  
  phi = sbar.EvalBasis(1, 1.0);
  EXPECT_NEAR(1.0, phi, tol);

  phi = sbar.EvalBasis(1, .5526);
  EXPECT_NEAR(.7763, phi, tol);

  //-------------------------------- k = 2 -------------------------------------
  phi = sbar.EvalBasis(2, 0.0);
  EXPECT_NEAR(-.6123724356957945, phi, tol);
  
  phi = sbar.EvalBasis(2, -1.0);
  EXPECT_NEAR(0.0, phi, tol);
  
  phi = sbar.EvalBasis(2, 1.0);
  EXPECT_NEAR(0.0, phi, tol);

  phi = sbar.EvalBasis(2, .5526);
  EXPECT_NEAR(-.4253742490940614, phi, tol);

  //-------------------------------- k = 3 -------------------------------------
  phi = sbar.EvalBasis(3, 0.0);
  EXPECT_NEAR(0.0, phi, tol);
  
  phi = sbar.EvalBasis(3, -1.0);
  EXPECT_NEAR(0.0, phi, tol);
  
  phi = sbar.EvalBasis(3, 1.0);
  EXPECT_NEAR(0.0, phi, tol);

  phi = sbar.EvalBasis(3, .5526);
  EXPECT_NEAR(-.3034634918835378, phi, tol);

}

TEST(BarElement, EvalBasisD) {
  realT dphi = -9.9e99;
  BarElement sbar;
  //-------------------------------- k = 0 -------------------------------------
  dphi = sbar.EvalBasisD(0, 0.0);
  EXPECT_NEAR(-.5, dphi, tol);
  
  dphi = sbar.EvalBasisD(0, -1.0);
  EXPECT_NEAR(-.5, dphi, tol);
  
  dphi = sbar.EvalBasisD(0, 1.0);
  EXPECT_NEAR(-.5, dphi, tol);

  dphi = sbar.EvalBasisD(0, .5526);
  EXPECT_NEAR(-.5, dphi, tol);
  
  //-------------------------------- k = 1 -------------------------------------
  dphi = sbar.EvalBasisD(1, 0.0);
  EXPECT_NEAR(.5, dphi, tol);
  
  dphi = sbar.EvalBasisD(1, -1.0);
  EXPECT_NEAR(.5, dphi, tol);
  
  dphi = sbar.EvalBasisD(1, 1.0);
  EXPECT_NEAR(.5, dphi, tol);

  dphi = sbar.EvalBasisD(1, .5526);
  EXPECT_NEAR(.5, dphi, tol);

  //------------------------------- k = 2 -------------------------------------
  dphi = sbar.EvalBasisD(2, -0.774596669241483  );
  EXPECT_NEAR(-0.9486832980505131, dphi, tol);
  
  dphi = sbar.EvalBasisD(2,0.0 );
  EXPECT_NEAR(0.0, dphi, tol);
  
  dphi = sbar.EvalBasisD(2, 0.774596669241483 );
  EXPECT_NEAR(-0.9486832980505131, dphi, tol);

 //-------------------------------- k = 3 -------------------------------------
  dphi = sbar.EvalBasisD(3, 0.0);
  EXPECT_NEAR(-.7905694150420948, dphi, tol);
  
  dphi = sbar.EvalBasisD(3, -1.0);
  EXPECT_NEAR(1.581138830084190, dphi, tol);
  
}

TEST(BarElement, IntegrateFunction) {
  BarElement bar;
  Array2D<realT> x(1,2);
  Array2D<realT> xmap(1,2);
  Array1D<realT> ans(1);
  realT a;
  
  xmap(0,0) = -3.0;
  xmap(0,1) = 3.0;
  x(0,0) = xmap(0,0);
  x(0,1) = xmap(0,1);
  
  bar.initialize(1, 1, 2);
  bar.IntegrateFunction(1, x, xmap, quadpoly, ans);
  a = ans(0);
  EXPECT_NEAR(18.0, a, 9.9e-15);

  
}

TEST(BarElement, IntegrateFunction2) {
  BarElement bar;
  Array2D<realT> x(1,2);
  Array2D<realT> xmap(1,2);
  Array1D<realT> ans(1);
  realT a;
  
  xmap(0,0) = -3.0;
  xmap(0,1) = 3.0;
  x(0,0) = xmap(0,0);
  x(0,1) = xmap(0,1);
  
  bar.initialize(1, 1, 3);
  bar.IntegrateFunction(1, x, xmap, p3poly, ans);
  a = ans(0);
  EXPECT_NEAR(24.0, a,tol);

  
}
TEST(BarElement, IntegrateFunction3) {
  BarElement bar;
  Array2D<realT> x(1,2);
  Array2D<realT> xmap(1,2);
  Array1D<realT> ans(1);
  realT a;
  
  xmap(0,0) = -3.0;
  xmap(0,1) = 3.0;
  x(0,0) = xmap(0,0);
  x(0,1) = xmap(0,1);
  
  bar.initialize(1, 1, 4);
  bar.IntegrateFunction(1, x, xmap, p4poly, ans);
  a = ans(0);
  EXPECT_NEAR(115.2, a, tol);

  
}
TEST(BarElement, IntegrateFunction4) {
  BarElement bar;
  Array2D<realT> x(1,2);
  Array2D<realT> xmap(1,2);
  Array1D<realT> ans(1);
  realT a;
  
  xmap(0,0) = -3.0;
  xmap(0,1) = 3.0;
  x(0,0) = xmap(0,0);
  x(0,1) = xmap(0,1);
  
  bar.initialize(1, 1, 5);
  bar.IntegrateFunction(1, x, xmap, p4poly, ans);
  a = ans(0);
  EXPECT_NEAR(115.2, a,tol);

  
}

TEST(BarElement, bc_geom) {
  BarElement bar;
  Array2D<realT> xmap(1,2);
  Array1D<realT> norm(1), tang(1);

  xmap(0,0) = -.25;
  xmap(0,1) = .678;
  
  bar.initialize(1, 1, 2);
  bar.ComputeFaceVectors(0, 0, xmap, norm, tang);
  EXPECT_NEAR(0.0, tang(0), 9.9e-15);
  EXPECT_NEAR(-1.0, norm(0), 9.9e-15);

  bar.ComputeFaceVectors(1, 0, xmap, norm, tang);
  EXPECT_NEAR(0.0, tang(0), 9.9e-15);
  EXPECT_NEAR(1.0, norm(0), 9.9e-15);
  
}

void quadpoly (const Array1D<realT>& x, Array1D<realT>& y)
{
  y(0) = x(0)*x(0);
}

void p3poly(const Array1D<realT>& x, Array1D<realT>& y)
{
  y(0) = x(0)*x(0)*x(0) + x(0)*x(0) + 1;
}

void p4poly(const Array1D<realT>& x, Array1D<realT>& y)
{
  y(0) = x(0)*x(0)*x(0)*x(0) + x(0)*x(0)*x(0) + x(0)*x(0);
}
