#include "Polynomials.h"
#include "gtest/gtest.h"
#include "my_incl.h"

//****************************************************************************80
//! \file polynomials_test.cpp
//! \brief Collection of google tests for namespace Polynomials
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! 
//****************************************************************************80
TEST(Polynomials, Jacobi) {
  double xi;
  double phi = 0.0;
  double dphi = 0.0;

  xi = .6;
  phi  = Polynomials::Jacobi(0, 1, 1, xi);
  dphi = Polynomials::JacobiD(0,1, 1, xi);
  EXPECT_EQ(1.0, phi);
  EXPECT_EQ(0.0, dphi);

  xi = .5;

  phi  = Polynomials::Jacobi(1, 1, 1, xi);
  dphi = Polynomials::JacobiD(1,1, 1, xi);
  EXPECT_EQ(1.0, phi);
  EXPECT_EQ(2.0, dphi);
  
  phi  = Polynomials::Jacobi(2, 1, 1, xi);
  dphi = Polynomials::JacobiD(2,1, 1, xi);
  EXPECT_EQ(.1875, phi);
  EXPECT_EQ(3.75, dphi);
 
  phi  = Polynomials::Jacobi(3, 1, 1, xi);
  dphi = Polynomials::JacobiD(3,1, 1, xi);
  EXPECT_EQ(-.625, phi);
  EXPECT_EQ(2.25, dphi);
  
  phi  = Polynomials::Jacobi(4, 1, 1, xi);
  dphi = Polynomials::JacobiD(4,1, 1, xi);
  EXPECT_EQ(-.7421875, phi);
  EXPECT_EQ(-2.1875, dphi);
  
  phi  = Polynomials::Jacobi(5, 1, 1, xi);
  dphi = Polynomials::JacobiD(5,1, 1, xi);
  EXPECT_EQ(-.1640625, phi);
  EXPECT_EQ(-5.390625, dphi);
  phi  = Polynomials::LobattoKern(5,  xi);
}
