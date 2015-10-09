#include "Elements/TriElement.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "gtest/gtest.h"
#include "my_incl.h"

void quadpoly (const Array1D<double>&, Array1D<double>& );

//****************************************************************************80
//! \file trielement_test.cpp
//! \brief Collection of google tests for namespace Polynomials
//! \details
//! \ryan
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! 
//****************************************************************************80
const double tol = 9.9e-15;
TEST(TriElement, EvalBasis) {
  double phi = -9.9e99;
  TriElement stri;

  //-------------------------------- k = 0 -------------------------------------
  phi = stri.EvalBasis(0,-1.0,-1.0);
  EXPECT_EQ(1.0,phi);
  
  phi = stri.EvalBasis(0,-1.0,0.0);
  EXPECT_EQ(0.5,phi);
  
  
  //-------------------------------- k = 1 -------------------------------------
  phi = stri.EvalBasis(1,1.0,1.0);
  EXPECT_EQ(1.0,phi);
  
  phi = stri.EvalBasis(1,0.5,1.0);
  EXPECT_EQ(0.75,phi);
  
  //-------------------------------- k = 2 -------------------------------------
  phi = stri.EvalBasis(2,1.0,1.0);
  EXPECT_EQ(1.0,phi);
  
  phi = stri.EvalBasis(2,1.0,0.5);
  EXPECT_EQ(0.75,phi);

  //-------------------------------- k = 3 -------------------------------------
  phi = stri.EvalBasis(3, -.3, .25);
  ASSERT_NEAR(-.02143303524935281, phi, tol);
    
  //-------------------------------- k = 4 -------------------------------------
  phi = stri.EvalBasis(4,-.3, .25);
  ASSERT_NEAR(-.5358258812338202, phi, tol);
  
  //-------------------------------- k = 5 -------------------------------------
  phi = stri.EvalBasis(5, -.3, .25);
  ASSERT_NEAR(-.03827327723098716, phi, tol);

  //-------------------------------- k = 9 -------------------------------------
  phi = stri.EvalBasis(9, -.3, .25);
  ASSERT_NEAR(.0328125, phi, tol);

  //------------------------------- k = 18 -------------------------------------
  phi = stri.EvalBasis(18, -.3, .25);
  ASSERT_NEAR(.01002438433271590, phi, tol);

}

TEST(TriElement, EvalBasisD) {
  TriElement stri;
  double dphidxi, dphideta;
  //-------------------------------- k = 0 -------------------------------------
  stri.EvalBasisD(0, -1.0, -1.0, dphidxi, dphideta);
  EXPECT_EQ(-0.5,dphidxi);
  EXPECT_EQ(-0.5,dphideta);
  
  //-------------------------------- k = 1 -------------------------------------
  stri.EvalBasisD(1, -1.0, -1.0, dphidxi, dphideta);
  EXPECT_EQ(0.5, dphidxi);
  EXPECT_EQ(0.0, dphideta);

  //-------------------------------- k = 2 -------------------------------------
  stri.EvalBasisD(2, -1.0, -1.0, dphidxi, dphideta);
  EXPECT_EQ(0.0, dphidxi);
  EXPECT_EQ(0.5, dphideta);

  //-------------------------------- k = 3 -------------------------------------
  stri.EvalBasisD(3, -.3, .25, dphidxi, dphideta);
  ASSERT_NEAR(.3980420832022664, dphidxi, tol);
  ASSERT_NEAR(.4286607049870562, dphideta, tol);

  //-------------------------------- k = 4 -------------------------------------
  stri.EvalBasisD(4, -.3, .25, dphidxi, dphideta);
  ASSERT_NEAR(-.7654655446197431, dphidxi, tol);
  ASSERT_NEAR(-.4286607049870562, dphideta, tol);
 
  //-------------------------------- k = 5 -------------------------------------
  stri.EvalBasisD(5, -.3, .25, dphidxi, dphideta);
  ASSERT_NEAR(.7654655446197431, dphidxi, tol);
  ASSERT_NEAR(.7348469228349534, dphideta, tol);

  //-------------------------------- k = 9 -------------------------------------
  stri.EvalBasisD(9, -.3, .25, dphidxi, dphideta);
  ASSERT_NEAR(-.609375, dphidxi, tol);
  ASSERT_NEAR(-.63, dphideta, tol);

  //-------------------------------- k = 18 ------------------------------------
  stri.EvalBasisD(18, -.3, .25, dphidxi, dphideta);
  ASSERT_NEAR(-.14857569635989634, dphidxi, tol);
  ASSERT_NEAR(-.11728529669277601, dphideta, tol);

  //-------------------------------- k = 19 ------------------------------------
  stri.EvalBasisD(19, -.3, .25, dphidxi, dphideta);
  ASSERT_NEAR(.15634765625, dphidxi, tol);
  ASSERT_NEAR(.17057031250, dphideta, tol);

}

TEST(TriElement, IntegrateFunction) {
  TriElement stri;
  Array2D<double> q(2,3);
  Array2D<double> xmap(2,3);
  Array1D<double> ans(2);
  double a;
  
  xmap(0,0) = 2.0;
  xmap(1,0) = 1.0;
  xmap(0,1) = 2.0;
  xmap(1,1) = 2.0;
  xmap(0,2) = 1.0;
  xmap(1,2) = 2.0;
  q(0,0) = xmap(0,0);
  q(1,0) = xmap(1,0);
  q(0,1) = xmap(0,1);
  q(1,1) = xmap(1,1);
  q(0,2) = xmap(0,2);
  q(1,2) = xmap(1,2);
  
  stri.initialize(1,1, 1);
  std::cout << "HERE" << std::endl;
  stri.IntegrateFunction(2, q, xmap, quadpoly, ans);
  a = ans(0);
  ASSERT_NEAR(5.0/3.0, a,9.9e-15);

}
TEST(TriElement2, IntegrateFunction2) {
  TriElement stri;
  Array2D<double> q(2,3);
  Array2D<double> xmap(2,3);
  Array1D<double> ans(2);
  double a;
  
  xmap(0,0) = 2.0;
  xmap(1,0) = 1.0;
  xmap(0,1) = 2.0;
  xmap(1,1) = 2.0;
  xmap(0,2) = 1.0;
  xmap(1,2) = 2.0;
  q(0,0) = xmap(0,0);
  q(1,0) = xmap(1,0);
  q(0,1) = xmap(0,1);
  q(1,1) = xmap(1,1);
  q(0,2) = xmap(0,2);
  q(1,2) = xmap(1,2);
  
  stri.initialize(1,1, 2);
  stri.IntegrateFunction(2, q, xmap, quadpoly, ans);
  a = ans(0);
  ASSERT_NEAR(5.0/3.0, a,9.9e-15);

}
TEST(TriElement3, IntegrateFunction3) {
  TriElement stri;
  Array2D<double> q(2,3);
  Array2D<double> xmap(2,3);
  Array1D<double> ans(2);
  double a;
  
  xmap(0,0) = 2.0;
  xmap(1,0) = 1.0;
  xmap(0,1) = 2.0;
  xmap(1,1) = 2.0;
  xmap(0,2) = 1.0;
  xmap(1,2) = 2.0;
  q(0,0) = xmap(0,0);
  q(1,0) = xmap(1,0);
  q(0,1) = xmap(0,1);
  q(1,1) = xmap(1,1);
  q(0,2) = xmap(0,2);
  q(1,2) = xmap(1,2);
  
  stri.initialize(1,1, 3);
  stri.IntegrateFunction(2, q, xmap, quadpoly, ans);
  a = ans(0);
  ASSERT_NEAR(5.0/3.0, a,9.9e-15);

}


TEST(TriElement, bc_geom) {
  TriElement tri;
  Array2D<double> xmap(2,3);
  Array1D<double> norm(2), tang(2);
  double Nx, Ny;
  
  xmap(0,0) = .5;
  xmap(0,1) = 1.15;
  xmap(0,2) = -.15;
  
  xmap(1,0) = .3;
  xmap(1,1) = .68;
  xmap(1,2) = 1.21;
  
  tri.initialize(1, 1, 2);
  
  //----------------------------- Side 0 ---------------------------------------
  Nx = half*(xmap(1,1) - xmap(1,0));
  Ny = -half*(xmap(0,1) - xmap(0,0));
  //---> point 0
  tri.ComputeFaceVectors(0, 0, xmap, norm, tang);
  
  ASSERT_NEAR(half*(xmap(0,1) - xmap(0,0)), tang(0), 9.9e-15);
  ASSERT_NEAR(half*(xmap(1,1) - xmap(1,0)), tang(1), 9.9e-15);
  
  ASSERT_NEAR(Nx, norm(0), 9.9e-15);
  ASSERT_NEAR(Ny, norm(1), 9.9e-15);
  
  //---> point 1
  tri.ComputeFaceVectors(0, 1, xmap, norm, tang);
  
  ASSERT_NEAR(half*(xmap(0,1) - xmap(0,0)), tang(0), 9.9e-15);
  ASSERT_NEAR(half*(xmap(1,1) - xmap(1,0)), tang(1), 9.9e-15);
  
  ASSERT_NEAR(Nx, norm(0), 9.9e-15);
  ASSERT_NEAR(Ny, norm(1), 9.9e-15);
  
  //------------------------------- Side 1 ------------------------------------
  //---> point 0
  Nx = half*(xmap(1,2) - xmap(1,1));
  Ny = -half*(xmap(0,2) - xmap(0,1));
  tri.ComputeFaceVectors(1, 0, xmap, norm, tang);
  
  ASSERT_NEAR(half*(xmap(0,2) - xmap(0,1)), tang(0), 9.9e-15);
  ASSERT_NEAR(half*(xmap(1,2) - xmap(1,1)), tang(1), 9.9e-15);

  ASSERT_NEAR(Nx, norm(0), 9.9e-15);
  ASSERT_NEAR(Ny, norm(1), 9.9e-15);
  
  //---> point 1;
  tri.ComputeFaceVectors(1, 1, xmap, norm, tang);
  
  ASSERT_NEAR(half*(xmap(0,2) - xmap(0,1)), tang(0), 9.9e-15);
  ASSERT_NEAR(half*(xmap(1,2) - xmap(1,1)), tang(1), 9.9e-15);
  
  ASSERT_NEAR(Nx, norm(0), 9.9e-15);
  ASSERT_NEAR(Ny, norm(1), 9.9e-15);
  
  //------------------------------- Side 2 ------------------------------------
  //---> point 0
  Nx = half*(xmap(1,0) - xmap(1,2));
  Ny = -half*(xmap(0,0) - xmap(0,2));
  tri.ComputeFaceVectors(2, 0, xmap, norm, tang);
  
  ASSERT_NEAR(half*(xmap(0,0) - xmap(0,2)), tang(0), 9.9e-15);
  ASSERT_NEAR(half*(xmap(1,0) - xmap(1,2)), tang(1), 9.9e-15);

  ASSERT_NEAR(Nx, norm(0), 9.9e-15);
  ASSERT_NEAR(Ny, norm(1), 9.9e-15);
  
  //---> point 1;
  tri.ComputeFaceVectors(2, 1, xmap, norm, tang);
  
  ASSERT_NEAR(half*(xmap(0,0) - xmap(0,2)), tang(0), 9.9e-15);
  ASSERT_NEAR(half*(xmap(1,0) - xmap(1,2)), tang(1), 9.9e-15);
  
  ASSERT_NEAR(Nx, norm(0), 9.9e-15);
  ASSERT_NEAR(Ny, norm(1), 9.9e-15);
  
}

TEST(Triangle, BasisSummation) {
 TriElement stri;
 int p = 2;
 stri.initialize(p, p, 1);

 double xi = -.3;
 double eta = 0.0;
 std::cout << stri.get_ndof() << std::endl;
 Array1D<double> basis(stri.get_ndof());

 
 for(int i = 0; i < stri.get_ndof(); i++){
   basis(i) = stri.EvalBasis(i,xi,eta);
 }

 std::cout << basis << std::endl;
 std::cout << basis(0) + basis(1) + basis(2) << std::endl;
 std::cout << basis(3) + basis(4) + basis(5) << std::endl;

}




void quadpoly (const Array1D<double>& x, Array1D<double>& y)
{
  y(0) = x(0) + x(1) ;
}
