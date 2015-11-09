/*
 * H1BarBasis.test.cpp
 *
 *  Created on: Nov 8, 2015
 *      Author: rabbit
 */
#include "DataStructures/Array1D.h"
#include "gtest/gtest.h"
#include "my_incl.h"
#include "BasisFunctions/Polynomials.h"
#include "BasisFunctions/H1BarBasis.h"
#include "Elements/Element.h"

const realT tol = 9.9e-13;
TEST(H1BarBasis, eval_basis) {
  realT phi = -9.9e99;
  H1BarBasis sbar;

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

TEST(H1BarBasis, EvalBasisD) {
  realT dphi = -9.9e99;
  H1BarBasis sbar;
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
  EXPECT_NEAR(0.9486832980505131, dphi, tol);

 //-------------------------------- k = 3 -------------------------------------
  dphi = sbar.EvalBasisD(3, 0.0);
  EXPECT_NEAR(-.7905694150420948, dphi, tol);

  dphi = sbar.EvalBasisD(3, -1.0);
  EXPECT_NEAR(1.581138830084190, dphi, tol);

}
