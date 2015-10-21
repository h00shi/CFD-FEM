/*
 *  * Poisson.test.cpp
 *
 *  Created on: Oct 16, 2015
 *      Author: rabbit
 */
#include "gtest/gtest.h"
#include "PDE/Poisson.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"

TEST(Poisson,OneD)
{
 Poisson poisson_equation;
 Array1D<realT> q(Poisson::nfld);
 Array2D<realT> dq(Poisson::nfld,1);
 Array1D<realT> flux(Poisson::nfld);
 Array1D<realT> accum(Poisson::nfld);
 Array2D<realT> D(1,1);
 Array1D<realT> c(1);
 Array1D<realT> norm(1);

 q(0) = 10.0;
 dq(0,0) = -185.1;
 D(0,0) = 1.0e-3;
 norm(0) = 2.0;
 c(0) = 0.0;
 realT Vr = 156.0;

 poisson_equation.CalcAccumulation(q,Vr,accum);
 poisson_equation.CalcNormFlux<realT>(1, norm, c, D, q, dq, flux);

 EXPECT_DOUBLE_EQ(10.0*156.0, accum(0));

 EXPECT_DOUBLE_EQ(-D(0,0)*dq(0,0)*norm(0), flux(0));

}



