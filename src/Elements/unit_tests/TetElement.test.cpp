#include "Elements/TetElement.h"
#include "DataStructures/Array1D.h"
#include "gtest/gtest.h"
#include "my_incl.h"

void quadpoly (const Array1D<double>& x, Array1D<double>& y)
{
  y(0) = x(0) + x(1) + x(2) ;
}
//****************************************************************************80
//! \file tetelement_test.cpp
//! \brief Collection of google tests for namespace Polynomials
//! \details
//! \ryan
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! 
//****************************************************************************80
const double tol = 9.9e-15;
TEST(TetElement, EvalBasis) {
  double phi = -9.9e99;
  TetElement stet;

  //-------------------------------- k = 0 -------------------------------------
  phi = stet.EvalBasis(0, -1.0, -1.0, -1.0);
  ASSERT_NEAR(1.0, phi, tol);
  
  phi = stet.EvalBasis(0, 1.0, -1.0, -1.0);
  ASSERT_NEAR(0.0, phi, tol);
  
  phi = stet.EvalBasis(0, -1.0, 1.0, -1.0);
  ASSERT_NEAR(0.0, phi, tol);
  
  phi = stet.EvalBasis(0, -1.0, -1.0, 1.0);
  ASSERT_NEAR(0.0, phi, tol);

  //-------------------------------- k = 1 -------------------------------------
  phi = stet.EvalBasis(1, -1.0, -1.0, -1.0);
  ASSERT_NEAR(0.0, phi, tol);
  
  phi = stet.EvalBasis(1, 1.0, -1.0, -1.0);
  ASSERT_NEAR(1.0, phi, tol);
  
  phi = stet.EvalBasis(1, -1.0, 1.0, -1.0);
  ASSERT_NEAR(0.0, phi, tol);
  
  phi = stet.EvalBasis(1, -1.0, -1.0, 1.0);
  ASSERT_NEAR(0.0, phi, tol);

  //-------------------------------- k = 2 -------------------------------------
  phi = stet.EvalBasis(2, -1.0, -1.0, -1.0);
  ASSERT_NEAR(0.0, phi, tol);
  
  phi = stet.EvalBasis(2, 1.0, -1.0, -1.0);
  ASSERT_NEAR(0.0, phi, tol);
  
  phi = stet.EvalBasis(2, -1.0, 1.0, -1.0);
  ASSERT_NEAR(1.0, phi, tol);
  
  phi = stet.EvalBasis(2, -1.0, -1.0, 1.0);
  ASSERT_NEAR(0.0, phi, tol);

  //-------------------------------- k = 3 -------------------------------------
  phi = stet.EvalBasis(3, -1.0, -1.0, -1.0);
  ASSERT_NEAR(0.0, phi, tol);
  
  phi = stet.EvalBasis(3, 1.0, -1.0, -1.0);
  ASSERT_NEAR(0.0, phi, tol);
  
  phi = stet.EvalBasis(3, -1.0, 1.0, -1.0);
  ASSERT_NEAR(0.0, phi, tol);
  
  phi = stet.EvalBasis(3, -1.0, -1.0, 1.0);
  ASSERT_NEAR(1.0, phi, tol);

  //------------------------------- k = 10 -------------------------------------
  const double xiq = 0.170820393249937;
  const double etaq = -0.723606797749979;
  const double zetaq  = -0.723606797749979;
  phi = stet.EvalBasis(10, xiq, etaq, zetaq);
  ASSERT_NEAR(-.1144122805635369, phi, tol);

  //------------------------------- k = 11 -------------------------------------
  phi = stet.EvalBasis(11, xiq, etaq, zetaq);
  ASSERT_NEAR(.1144122805635369, phi, tol);

  //------------------------------- k = 12 -------------------------------------
  phi = stet.EvalBasis(12, xiq, etaq, zetaq);
  ASSERT_NEAR(0.0, phi, tol);

  //------------------------------- k = 13 -------------------------------------
  phi = stet.EvalBasis(13, xiq, etaq, zetaq);
  ASSERT_NEAR(0.0, phi, tol);

  //------------------------------- k = 14 -------------------------------------
  phi = stet.EvalBasis(14, xiq, etaq, zetaq);
  ASSERT_NEAR(.1144122805635369, phi, tol);

  //------------------------------- k = 15 -------------------------------------
  phi = stet.EvalBasis(15, xiq, etaq, zetaq);
  ASSERT_NEAR(0.0, phi, tol);

  //------------------------------- k = 16 -------------------------------------
  phi = stet.EvalBasis(16, xiq, etaq, zetaq);
  ASSERT_NEAR(.06708203932499370, phi, tol);

  //------------------------------- k = 17 -------------------------------------
  phi = stet.EvalBasis(17, xiq, etaq, zetaq);
  ASSERT_NEAR(.06708203932499370, phi, tol); 

  //------------------------------- k = 18 -------------------------------------
  phi = stet.EvalBasis(18, xiq, etaq, zetaq);
  ASSERT_NEAR(.01583592135001261, phi, tol); 

  //------------------------------- k = 19 -------------------------------------
  phi = stet.EvalBasis(19, xiq, etaq, zetaq);
  ASSERT_NEAR(.06708203932499364, phi, tol); 

  //------------------------------- k = 34 -------------------------------------
  phi = stet.EvalBasis(34, xiq, etaq, zetaq);
  ASSERT_NEAR(-.02270801874201362, phi, tol); 

  phi = stet.EvalBasis(52, -1.0, -1.0, -1.0
);
  

}

TEST(TetElement, EvalBasisD) {
  TetElement stet;
  double dphidxi,dphideta,dphidzeta;

  //-------------------------------- k = 0 -------------------------------------
  stet.EvalBasisD(0, -1.0, -1.0, -1.0, dphidxi, dphideta, dphidzeta);
  EXPECT_EQ(-0.5, dphidxi);
  EXPECT_EQ(-0.5, dphideta);
  EXPECT_EQ(-0.5, dphidzeta);

 //-------------------------------- k = 1 -------------------------------------
  stet.EvalBasisD(1, -1.0, -1.0, -1.0, dphidxi, dphideta, dphidzeta);
  ASSERT_NEAR(0.5, dphidxi, tol);
  ASSERT_NEAR(0.0, dphideta, tol);
  ASSERT_NEAR(0.0, dphidzeta, tol);
 
  //-------------------------------- k = 2 -------------------------------------
  stet.EvalBasisD(2, -1.0, -1.0, -1.0, dphidxi, dphideta, dphidzeta);
  ASSERT_NEAR(0.0, dphidxi, tol);
  ASSERT_NEAR(0.5, dphideta, tol);
  ASSERT_NEAR(0.0, dphidzeta, tol);

  //-------------------------------- k = 3 -------------------------------------
  stet.EvalBasisD(3, -1.0, -1.0, -1.0, dphidxi, dphideta, dphidzeta);
  ASSERT_NEAR(0.0, dphidxi, tol);
  ASSERT_NEAR(0.0, dphideta, tol);
  ASSERT_NEAR(0.5, dphidzeta, tol);

  //------------------------------- k = 10 -------------------------------------
  const double xiq = 0.170820393249937;
  const double etaq = -0.723606797749979;
  const double zetaq  = -0.723606797749979;
  stet.EvalBasisD(10, xiq, etaq, zetaq, dphidxi, dphideta, dphidzeta);
  ASSERT_NEAR(0.0603941292159917, dphidxi, tol);
  ASSERT_NEAR(0.2860307014088421, dphideta, tol);
  ASSERT_NEAR(0.2860307014088421, dphidzeta, tol);
  //------------------------------- k = 11 -------------------------------------
  stet.EvalBasisD(11, xiq, etaq, zetaq, dphidxi, dphideta, dphidzeta);
  ASSERT_NEAR(0.2256365721928506, dphidxi, tol);
  ASSERT_NEAR(0.2860307014088421, dphideta, tol);
  ASSERT_NEAR(0.0, dphidzeta, tol);

  //------------------------------- k = 12 -------------------------------------
  stet.EvalBasisD(12, xiq, etaq, zetaq, dphidxi, dphideta, dphidzeta);
  ASSERT_NEAR(0.03019706460799585, dphidxi, tol);
  ASSERT_NEAR(0.06039412921599151, dphideta, tol);
  ASSERT_NEAR(0.03019706460799585, dphidzeta, tol);

  //------------------------------- k = 13 -------------------------------------
  stet.EvalBasisD(13, xiq, etaq, zetaq, dphidxi, dphideta, dphidzeta);
  ASSERT_NEAR(-0.03019706460799577, dphidxi, tol);
  ASSERT_NEAR(-0.03019706460799577, dphideta, tol);
  ASSERT_NEAR(-0.06039412921599154, dphidzeta, tol);

  //------------------------------- k = 14 -------------------------------------
  stet.EvalBasisD(14, xiq, etaq, zetaq, dphidxi, dphideta, dphidzeta);
  ASSERT_NEAR(0.2256365721928506, dphidxi, tol);
  ASSERT_NEAR(0.0, dphideta, tol);
  ASSERT_NEAR(0.2860307014088421, dphidzeta, tol);

  //------------------------------- k = 15 -------------------------------------
  stet.EvalBasisD(15, xiq, etaq, zetaq, dphidxi, dphideta, dphidzeta);
  ASSERT_NEAR(0.0, dphidxi, tol);
  ASSERT_NEAR(0.03019706460799577, dphideta, tol);
  ASSERT_NEAR(-0.03019706460799577, dphidzeta, tol);

  //------------------------------- k = 16 -------------------------------------
  stet.EvalBasisD(16, xiq, etaq, zetaq, dphidxi, dphideta, dphidzeta);
  ASSERT_NEAR(-0.1854101966249684, dphidxi, tol);
  ASSERT_NEAR(-0.2427050983124842, dphideta, tol);
  ASSERT_NEAR(0.0, dphidzeta, tol);

  //------------------------------- k = 17 -------------------------------------
  stet.EvalBasisD(17, xiq, etaq, zetaq, dphidxi, dphideta, dphidzeta);
  ASSERT_NEAR(0.05729490168751575, dphidxi, tol);
  ASSERT_NEAR(0.2427050983124842, dphideta, tol);
  ASSERT_NEAR(0.2427050983124842, dphidzeta, tol);

  //------------------------------- k = 18 -------------------------------------
  stet.EvalBasisD(18, xiq, etaq, zetaq, dphidxi, dphideta, dphidzeta);
  ASSERT_NEAR(-0.05729490168751575, dphidxi, tol);
  ASSERT_NEAR(0.0, dphideta, tol);
  ASSERT_NEAR(0.0, dphidzeta, tol);

  //------------------------------- k = 19 -------------------------------------
  stet.EvalBasisD(19, xiq, etaq, zetaq, dphidxi, dphideta, dphidzeta);
  ASSERT_NEAR(-0.1854101966249684, dphidxi, tol);
  ASSERT_NEAR(0.0, dphideta, tol);
  ASSERT_NEAR(-0.2427050983124842, dphidzeta, tol);

  //------------------------------- k = 34 -------------------------------------
  stet.EvalBasisD(34, xiq, etaq, zetaq, dphidxi, dphideta, dphidzeta);
  ASSERT_NEAR(0.06276342016858640, dphidxi, tol);
  ASSERT_NEAR(0.0, dphideta, tol);
  ASSERT_NEAR(0.0, dphidzeta, tol);

}

TEST(TetElement, IntegrateFunction) {
  TetElement stet;
  Array2D<double> q(3,4);
  Array2D<double> xmap(3,4);
  Array1D<double> ans(3);
  double a;
  
  xmap(0,0) = 0.0;
  xmap(1,0) = 0.0;
  xmap(2,0) = 0.0;
  xmap(0,1) = 1.0;
  xmap(1,1) = 0.0;
  xmap(2,1) = 0.0;
  xmap(0,2) = 0.0;
  xmap(1,2) = 1.0;
  xmap(2,2) = 0.0;
  xmap(0,3) = 0.0;
  xmap(1,3) = 0.0;
  xmap(2,3) = 1.0;
  q(0,0) = xmap(0,0);
  q(1,0) = xmap(1,0);
  q(2,0) = xmap(2,0);
  q(0,1) = xmap(0,1);
  q(1,1) = xmap(1,1);
  q(2,1) = xmap(2,1);
  q(0,2) = xmap(0,2);
  q(1,2) = xmap(1,2);
  q(2,2) = xmap(2,2);
  q(0,3) = xmap(0,3);
  q(1,3) = xmap(1,3);
  q(2,3) = xmap(2,3);
 
  stet.initialize(1, 1, 1);
  stet.IntegrateFunction(3, q, xmap, quadpoly, ans);
  a = ans(0);
  ASSERT_NEAR(0.125, a,9.9e-15);

}

TEST(TetElement, IntegrateFunction2) {
  TetElement stet;
  Array2D<double> q(3,4);
  Array2D<double> xmap(3,4);
  Array1D<double> ans(3);
  double a;
  
  xmap(0,0) = 0.0;
  xmap(1,0) = 0.0;
  xmap(2,0) = 0.0;
  xmap(0,1) = 1.0;
  xmap(1,1) = 0.0;
  xmap(2,1) = 0.0;
  xmap(0,2) = 0.0;
  xmap(1,2) = 1.0;
  xmap(2,2) = 0.0;
  xmap(0,3) = 0.0;
  xmap(1,3) = 0.0;
  xmap(2,3) = 1.0;
  q(0,0) = xmap(0,0);
  q(1,0) = xmap(1,0);
  q(2,0) = xmap(2,0);
  q(0,1) = xmap(0,1);
  q(1,1) = xmap(1,1);
  q(2,1) = xmap(2,1);
  q(0,2) = xmap(0,2);
  q(1,2) = xmap(1,2);
  q(2,2) = xmap(2,2);
  q(0,3) = xmap(0,3);
  q(1,3) = xmap(1,3);
  q(2,3) = xmap(2,3);
  
  stet.initialize(1, 1, 2);
  stet.IntegrateFunction(3, q, xmap, quadpoly, ans);
  a = ans(0);
  ASSERT_NEAR(0.125, a,9.9e-15);

}
TEST(TetElement, IntegrateFunction3) {
  TetElement stet;
  Array2D<double> q(3,4);
  Array2D<double> xmap(3,4);
  Array1D<double> ans(3);
  double a;
  
  xmap(0,0) = 0.0;
  xmap(1,0) = 0.0;
  xmap(2,0) = 0.0;
  xmap(0,1) = 1.0;
  xmap(1,1) = 0.0;
  xmap(2,1) = 0.0;
  xmap(0,2) = 0.0;
  xmap(1,2) = 1.0;
  xmap(2,2) = 0.0;
  xmap(0,3) = 0.0;
  xmap(1,3) = 0.0;
  xmap(2,3) = 1.0;
  q(0,0) = xmap(0,0);
  q(1,0) = xmap(1,0);
  q(2,0) = xmap(2,0);
  q(0,1) = xmap(0,1);
  q(1,1) = xmap(1,1);
  q(2,1) = xmap(2,1);
  q(0,2) = xmap(0,2);
  q(1,2) = xmap(1,2);
  q(2,2) = xmap(2,2);
  q(0,3) = xmap(0,3);
  q(1,3) = xmap(1,3);
  q(2,3) = xmap(2,3);
  
  stet.initialize(1, 1, 3);
  stet.IntegrateFunction(3, q, xmap, quadpoly, ans);
  a = ans(0);
  ASSERT_NEAR(0.125, a,9.9e-15);

}
TEST(TetElement, IntegrateFunction4) {
  TetElement stet;
  Array2D<double> q(3,4);
  Array2D<double> xmap(3,4);
  Array1D<double> ans(3);
  double a;
  
  xmap(0,0) = 0.0;
  xmap(1,0) = 0.0;
  xmap(2,0) = 0.0;
  xmap(0,1) = 1.0;
  xmap(1,1) = 0.0;
  xmap(2,1) = 0.0;
  xmap(0,2) = 0.0;
  xmap(1,2) = 1.0;
  xmap(2,2) = 0.0;
  xmap(0,3) = 0.0;
  xmap(1,3) = 0.0;
  xmap(2,3) = 1.0;
  q(0,0) = xmap(0,0);
  q(1,0) = xmap(1,0);
  q(2,0) = xmap(2,0);
  q(0,1) = xmap(0,1);
  q(1,1) = xmap(1,1);
  q(2,1) = xmap(2,1);
  q(0,2) = xmap(0,2);
  q(1,2) = xmap(1,2);
  q(2,2) = xmap(2,2);
  q(0,3) = xmap(0,3);
  q(1,3) = xmap(1,3);
  q(2,3) = xmap(2,3);
  
  stet.initialize(1, 1, 4);
  stet.IntegrateFunction(3, q, xmap, quadpoly, ans);
  a = ans(0);
  ASSERT_NEAR(0.125, a,9.9e-15);

}

TEST(TetElement, IntegrateFunction5) {
  TetElement stet;
  Array2D<double> q(3,4);
  Array2D<double> xmap(3,4);
  Array1D<double> ans(3);
  double a;
  
  xmap(0,0) = 0.0;
  xmap(1,0) = 0.0;
  xmap(2,0) = 0.0;
  xmap(0,1) = 1.0;
  xmap(1,1) = 0.0;
  xmap(2,1) = 0.0;
  xmap(0,2) = 0.0;
  xmap(1,2) = 1.0;
  xmap(2,2) = 0.0;
  xmap(0,3) = 0.0;
  xmap(1,3) = 0.0;
  xmap(2,3) = 1.0;
  q(0,0) = xmap(0,0);
  q(1,0) = xmap(1,0);
  q(2,0) = xmap(2,0);
  q(0,1) = xmap(0,1);
  q(1,1) = xmap(1,1);
  q(2,1) = xmap(2,1);
  q(0,2) = xmap(0,2);
  q(1,2) = xmap(1,2);
  q(2,2) = xmap(2,2);
  q(0,3) = xmap(0,3);
  q(1,3) = xmap(1,3);
  q(2,3) = xmap(2,3);
  
  stet.initialize(1, 1, 5);
  stet.IntegrateFunction(3, q, xmap, quadpoly, ans);
  a = ans(0);
  ASSERT_NEAR(0.125, a,9.9e-15);

}
TEST(TetElement, bc_geom) {
  TetElement tet;
  Array2D<double> xmap(3,4);
  Array1D<double> norm(3), tang(3);
  double Nx, Ny, Nz, Tx, Ty, Tz;
  
  xmap(0,0) = .5;
  xmap(0,1) = 1.2;
  xmap(0,2) = .75;
  xmap(0,3) = .98;
  
  xmap(1,0) = -.1;
  xmap(1,1) = .3;
  xmap(1,2) = 1.1;
  xmap(1,3) = .65;
  
  xmap(2,0) = -.12;
  xmap(2,1) = .05;
  xmap(2,2) = .14;
  xmap(2,3) = .98;

  tet.initialize(1, 1, 2);
  
  //----------------------------- Side 0 ---------------------------------------
  Tx = half*.7;
  Ty = half*.4;
  Tz = half*.17;
  Nx = half*0.156250000000000;
  Ny = half*(-.3442);
  Nz = half*(.1665);

  //---> point 0
  tet.ComputeFaceVectors(0, 0, xmap, norm, tang);
  ASSERT_NEAR(Tx, tang(0),9.9e-15);
  ASSERT_NEAR(Ty, tang(1),9.9e-15);
  ASSERT_NEAR(Tz, tang(2),9.9e-15);

  ASSERT_NEAR(Nx, norm(0),9.9e-15);
  ASSERT_NEAR(Ny, norm(1),9.9e-15);
  ASSERT_NEAR(Nz, norm(2),9.9e-15);
   
  //---> point 1
  tet.ComputeFaceVectors(0, 1, xmap, norm, tang);
  ASSERT_NEAR(Tx, tang(0),9.9e-15);
  ASSERT_NEAR(Ty, tang(1),9.9e-15);
  ASSERT_NEAR(Tz, tang(2),9.9e-15);

  ASSERT_NEAR(Nx, norm(0),9.9e-15);
  ASSERT_NEAR(Ny, norm(1),9.9e-15);
  ASSERT_NEAR(Nz, norm(2),9.9e-15);

  //---> point 2
  tet.ComputeFaceVectors(0, 2, xmap, norm, tang);
  ASSERT_NEAR(Tx, tang(0),9.9e-15);
  ASSERT_NEAR(Ty, tang(1),9.9e-15);
  ASSERT_NEAR(Tz, tang(2),9.9e-15);

  ASSERT_NEAR(Nx, norm(0),9.9e-15);
  ASSERT_NEAR(Ny, norm(1),9.9e-15);
  ASSERT_NEAR(Nz, norm(2),9.9e-15);
  
  //------------------------------- Side 1 ------------------------------------
  Tx = half*(-.45);
  Ty = half*.8;
  Tz = half*.09;
  Nx = half*.35625;
  Ny = half*(.19935);
  Nz = half*(.00925);
  //---> point 0
  tet.ComputeFaceVectors(1, 0, xmap, norm, tang);
  ASSERT_NEAR(Tx, tang(0),9.9e-15);
  ASSERT_NEAR(Ty, tang(1),9.9e-15);
  ASSERT_NEAR(Tz, tang(2),9.9e-15);

  ASSERT_NEAR(Nx, norm(0),9.9e-15);
  ASSERT_NEAR(Ny, norm(1),9.9e-15);
  ASSERT_NEAR(Nz, norm(2),9.9e-15);
   
  //---> point 1
  tet.ComputeFaceVectors(1, 1, xmap, norm, tang);
  ASSERT_NEAR(Tx, tang(0),9.9e-15);
  ASSERT_NEAR(Ty, tang(1),9.9e-15);
  ASSERT_NEAR(Tz, tang(2),9.9e-15);

  ASSERT_NEAR(Nx, norm(0),9.9e-15);
  ASSERT_NEAR(Ny, norm(1),9.9e-15);
  ASSERT_NEAR(Nz, norm(2),9.9e-15);

  //---> point 2
  tet.ComputeFaceVectors(1, 2, xmap, norm, tang);
  ASSERT_NEAR(Tx, tang(0),9.9e-15);
  ASSERT_NEAR(Ty, tang(1),9.9e-15);
  ASSERT_NEAR(Tz, tang(2),9.9e-15);

  ASSERT_NEAR(Nx, norm(0),9.9e-15);
  ASSERT_NEAR(Ny, norm(1),9.9e-15);
  ASSERT_NEAR(Nz, norm(2),9.9e-15);
  //------------------------------- Side 2 ------------------------------------
  Tx = half*(-.25);
  Ty = half*(-1.2);
  Tz = half*(-.26);
  Nx = half*(-.5625);
  Ny = half*(.0751);
  Nz = half*(.19425);
  //---> point 0
  tet.ComputeFaceVectors(2, 0, xmap, norm, tang);
  ASSERT_NEAR(Tx, tang(0),9.9e-15);
  ASSERT_NEAR(Ty, tang(1),9.9e-15);
  ASSERT_NEAR(Tz, tang(2),9.9e-15);

  ASSERT_NEAR(Nx, norm(0),9.9e-15);
  ASSERT_NEAR(Ny, norm(1),9.9e-15);
  ASSERT_NEAR(Nz, norm(2),9.9e-15);
   
  //---> point 1
  tet.ComputeFaceVectors(2, 1, xmap, norm, tang);
  ASSERT_NEAR(Tx, tang(0),9.9e-15);
  ASSERT_NEAR(Ty, tang(1),9.9e-15);
  ASSERT_NEAR(Tz, tang(2),9.9e-15);

  ASSERT_NEAR(Nx, norm(0),9.9e-15);
  ASSERT_NEAR(Ny, norm(1),9.9e-15);
  ASSERT_NEAR(Nz, norm(2),9.9e-15);

  //---> point 2
  tet.ComputeFaceVectors(2, 2, xmap, norm, tang);
  ASSERT_NEAR(Tx, tang(0),9.9e-15);
  ASSERT_NEAR(Ty, tang(1),9.9e-15);
  ASSERT_NEAR(Tz, tang(2),9.9e-15);

  ASSERT_NEAR(Nx, norm(0),9.9e-15);
  ASSERT_NEAR(Ny, norm(1),9.9e-15);
  ASSERT_NEAR(Nz, norm(2),9.9e-15);

  //------------------------------- Side 3 ------------------------------------
  Tx = half*(.2500);
  Ty = half*(1.200);
  Tz = half*(.2600);
  Nx = half*(.05);
  Ny = half*(.06975);
  Nz = half*(-.3700);
  //---> point 0
  tet.ComputeFaceVectors(3, 0, xmap, norm, tang);
  ASSERT_NEAR(Tx, tang(0),9.9e-15);
  ASSERT_NEAR(Ty, tang(1),9.9e-15);
  ASSERT_NEAR(Tz, tang(2),9.9e-15);

  ASSERT_NEAR(Nx, norm(0),9.9e-15);
  ASSERT_NEAR(Ny, norm(1),9.9e-15);
  ASSERT_NEAR(Nz, norm(2),9.9e-15);
   
  //---> point 1
  tet.ComputeFaceVectors(3, 1, xmap, norm, tang);
  ASSERT_NEAR(Tx, tang(0),9.9e-15);
  ASSERT_NEAR(Ty, tang(1),9.9e-15);
  ASSERT_NEAR(Tz, tang(2),9.9e-15);

  ASSERT_NEAR(Nx, norm(0),9.9e-15);
  ASSERT_NEAR(Ny, norm(1),9.9e-15);
  ASSERT_NEAR(Nz, norm(2),9.9e-15);

  //---> point 2
  tet.ComputeFaceVectors(3, 2, xmap, norm, tang);
  ASSERT_NEAR(Tx, tang(0),9.9e-15);
  ASSERT_NEAR(Ty, tang(1),9.9e-15);
  ASSERT_NEAR(Tz, tang(2),9.9e-15);

  ASSERT_NEAR(Nx, norm(0),9.9e-15);
  ASSERT_NEAR(Ny, norm(1),9.9e-15);
  ASSERT_NEAR(Nz, norm(2),9.9e-15);
}




