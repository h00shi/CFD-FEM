#include "DataStructures/SquareMatrix.h"
#include "DataStructures/Array1D.h"
#include "gtest/gtest.h"
#include "Surreal/Forward/Surreal.h"
#define MATRIX_N 4
//****************************************************************************80
//! \file squarematrix_test.cpp
//! \brief Collection of google tests for block matrix
//! \details
//! \ryan
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! 
//****************************************************************************80
TEST(SquareMatrix, eval_LU) {

  SquareMatrix<double> unit(MATRIX_N);
  unit(0,0) = 4.;
  unit(0,1) = 6.;
  unit(0,2) = 7.;
  unit(0,3) = 8.;

  unit(1,0) = 1.;
  unit(1,1) = 0.;
  unit(1,2) = 2.;
  unit(1,3) = 4.;

  unit(2,0) = 6.;
  unit(2,1) = 5.;
  unit(2,2) = 7.;
  unit(2,3) = 8.;

  unit(3,0) = 9.;
  unit(3,1) = 8.;
  unit(3,2) = 7.;
  unit(3,3) = 6.;

  Array1D<int> pivot(MATRIX_N);
  unit.LUdecomp(pivot);
 

  EXPECT_EQ(3,pivot(0));
  EXPECT_EQ(0,pivot(1));
  EXPECT_EQ(2,pivot(2));
  EXPECT_EQ(1,pivot(3));
 
  Array2D<double> ExpLUResult(MATRIX_N,MATRIX_N);
  ExpLUResult(0,0) = 9.;
  ExpLUResult(0,1) = 8.;
  ExpLUResult(0,2) = 7.;
  ExpLUResult(0,3) = 6.;

  ExpLUResult(1,0) = 0.444444444444444;
  ExpLUResult(1,1) = 2.44444444444444;
  ExpLUResult(1,2) = 3.88888888888889;
  ExpLUResult(1,3) = 5.33333333333333;

  ExpLUResult(2,0) = 0.666666666666667;
  ExpLUResult(2,1) = -0.136363636363636;
  ExpLUResult(2,2) = 2.86363636363636;
  ExpLUResult(2,3) = 4.72727272727273;

  ExpLUResult(3,0) = 0.111111111111111;
  ExpLUResult(3,1) = -0.363636363636364;
  ExpLUResult(3,2) = 0.920634920634921;
  ExpLUResult(3,3) = 0.920634920634922;
  
  double small = 1.e-13;

  EXPECT_NEAR(ExpLUResult(0,0),unit(pivot(0),0),small);
  EXPECT_NEAR(ExpLUResult(0,1),unit(pivot(0),1),small);
  EXPECT_NEAR(ExpLUResult(0,2),unit(pivot(0),2),small);
  EXPECT_NEAR(ExpLUResult(0,3),unit(pivot(0),3),small);
   
  EXPECT_NEAR(ExpLUResult(1,0),unit(pivot(1),0),small);
  EXPECT_NEAR(ExpLUResult(1,1),unit(pivot(1),1),small);
  EXPECT_NEAR(ExpLUResult(1,2),unit(pivot(1),2),small);
  EXPECT_NEAR(ExpLUResult(1,3),unit(pivot(1),3),small);
  
  EXPECT_NEAR(ExpLUResult(2,0),unit(pivot(2),0),small);
  EXPECT_NEAR(ExpLUResult(2,1),unit(pivot(2),1),small);
  EXPECT_NEAR(ExpLUResult(2,2),unit(pivot(2),2),small);
  EXPECT_NEAR(ExpLUResult(2,3),unit(pivot(2),3),small);

  EXPECT_NEAR(ExpLUResult(3,0),unit(pivot(3),0),small);
  EXPECT_NEAR(ExpLUResult(3,1),unit(pivot(3),1),small);
  EXPECT_NEAR(ExpLUResult(3,2),unit(pivot(3),2),small);
  EXPECT_NEAR(ExpLUResult(3,3),unit(pivot(3),3),small);
}

TEST(SquareMatrix, eval_inv1) {

  SquareMatrix<double> unit(MATRIX_N);
  unit(0,0) = 4.;
  unit(0,1) = 6.;
  unit(0,2) = 7.;
  unit(0,3) = 8.;

  unit(1,0) = 1.;
  unit(1,1) = 0.;
  unit(1,2) = 2.;
  unit(1,3) = 4.;

  unit(2,0) = 6.;
  unit(2,1) = 5.;
  unit(2,2) = 7.;
  unit(2,3) = 8.;

  unit(3,0) = 9.;
  unit(3,1) = 8.;
  unit(3,2) = 7.;
  unit(3,3) = 6.;

  unit.invert();

  Array2D<double> ExpINVResult(MATRIX_N,MATRIX_N);
  ExpINVResult(0,0) = -0.275862068965517;
  ExpINVResult(0,1) = 0.241379310344828;
  ExpINVResult(0,2) = 0.0;
  ExpINVResult(0,3) = 0.206896551724138;

  ExpINVResult(1,0) = 0.448275862068965;
  ExpINVResult(1,1) = 0.482758620689654;
  ExpINVResult(1,2) = -0.999999999999999;
  ExpINVResult(1,3) = 0.413793103448275;

  ExpINVResult(2,0) = -0.379310344827586;
  ExpINVResult(2,1) = -1.79310344827586;
  ExpINVResult(2,2) = 2.00000000000000;
  ExpINVResult(2,3) = -0.965517241379309;

  ExpINVResult(3,0) = 0.258620689655172;
  ExpINVResult(3,1) = 1.08620689655172;
  ExpINVResult(3,2) = -0.999999999999999;
  ExpINVResult(3,3) = 0.431034482758620;
  
  double small = 1.e-13;

  EXPECT_NEAR(ExpINVResult(0,0),unit(0,0),small);
  EXPECT_NEAR(ExpINVResult(0,1),unit(0,1),small);
  EXPECT_NEAR(ExpINVResult(0,2),unit(0,2),small);
  EXPECT_NEAR(ExpINVResult(0,3),unit(0,3),small);
   
  EXPECT_NEAR(ExpINVResult(1,0),unit(1,0),small);
  EXPECT_NEAR(ExpINVResult(1,1),unit(1,1),small);
  EXPECT_NEAR(ExpINVResult(1,2),unit(1,2),small);
  EXPECT_NEAR(ExpINVResult(1,3),unit(1,3),small);
  
  EXPECT_NEAR(ExpINVResult(2,0),unit(2,0),small);
  EXPECT_NEAR(ExpINVResult(2,1),unit(2,1),small);
  EXPECT_NEAR(ExpINVResult(2,2),unit(2,2),small);
  EXPECT_NEAR(ExpINVResult(2,3),unit(2,3),small);

  EXPECT_NEAR(ExpINVResult(3,0),unit(3,0),small);
  EXPECT_NEAR(ExpINVResult(3,1),unit(3,1),small);
  EXPECT_NEAR(ExpINVResult(3,2),unit(3,2),small);
  EXPECT_NEAR(ExpINVResult(3,3),unit(3,3),small);
}

TEST(SquareMatrix, eval_LU2) {

  SquareMatrix<double> unit(MATRIX_N);
  unit(0,0) = 4.;
  unit(0,1) = 6.;
  unit(0,2) = 7.;
  unit(0,3) = 8.;

  unit(1,0) = 1.;
  unit(1,1) = 0.;
  unit(1,2) = 2.;
  unit(1,3) = 4.;

  unit(2,0) = 6.;
  unit(2,1) = 5.;
  unit(2,2) = 7.;
  unit(2,3) = 8.;

  unit(3,0) = 9.;
  unit(3,1) = 8.;
  unit(3,2) = 7.;
  unit(3,3) = 6.;

  Array1D<int> pivot(MATRIX_N);
  unit.LUdecomp(pivot);
 

  EXPECT_EQ(3,pivot(0));
  EXPECT_EQ(0,pivot(1));
  EXPECT_EQ(2,pivot(2));
  EXPECT_EQ(1,pivot(3));
 
  Array2D<double> ExpLUResult(MATRIX_N,MATRIX_N);
  ExpLUResult(0,0) = 9.;
  ExpLUResult(0,1) = 8.;
  ExpLUResult(0,2) = 7.;
  ExpLUResult(0,3) = 6.;

  ExpLUResult(1,0) = 0.444444444444444;
  ExpLUResult(1,1) = 2.44444444444444;
  ExpLUResult(1,2) = 3.88888888888889;
  ExpLUResult(1,3) = 5.33333333333333;

  ExpLUResult(2,0) = 0.666666666666667;
  ExpLUResult(2,1) = -0.136363636363636;
  ExpLUResult(2,2) = 2.86363636363636;
  ExpLUResult(2,3) = 4.72727272727273;

  ExpLUResult(3,0) = 0.111111111111111;
  ExpLUResult(3,1) = -0.363636363636364;
  ExpLUResult(3,2) = 0.920634920634921;
  ExpLUResult(3,3) = 0.920634920634922;
  
  double small = 1.e-13;

  EXPECT_NEAR(ExpLUResult(0,0),unit(pivot(0),0),small);
  EXPECT_NEAR(ExpLUResult(0,1),unit(pivot(0),1),small);
  EXPECT_NEAR(ExpLUResult(0,2),unit(pivot(0),2),small);
  EXPECT_NEAR(ExpLUResult(0,3),unit(pivot(0),3),small);
   
  EXPECT_NEAR(ExpLUResult(1,0),unit(pivot(1),0),small);
  EXPECT_NEAR(ExpLUResult(1,1),unit(pivot(1),1),small);
  EXPECT_NEAR(ExpLUResult(1,2),unit(pivot(1),2),small);
  EXPECT_NEAR(ExpLUResult(1,3),unit(pivot(1),3),small);
  
  EXPECT_NEAR(ExpLUResult(2,0),unit(pivot(2),0),small);
  EXPECT_NEAR(ExpLUResult(2,1),unit(pivot(2),1),small);
  EXPECT_NEAR(ExpLUResult(2,2),unit(pivot(2),2),small);
  EXPECT_NEAR(ExpLUResult(2,3),unit(pivot(2),3),small);

  EXPECT_NEAR(ExpLUResult(3,0),unit(pivot(3),0),small);
  EXPECT_NEAR(ExpLUResult(3,1),unit(pivot(3),1),small);
  EXPECT_NEAR(ExpLUResult(3,2),unit(pivot(3),2),small);
  EXPECT_NEAR(ExpLUResult(3,3),unit(pivot(3),3),small);
}

// TEST(SquareMatrix, eval_LU3) {
//   Array1D<int> pivot(5);
//   SquareMatrix<double> unit(5), ExpLUResult(5);
  
//   unit(0,0) = 3.88762;
//   unit(0,1) = 1.69154;
//   unit(0,2) = -2.07394;
//   unit(0,3) = 1.6199;
//   unit(0,4) = 4.50208;
  
//   unit(1,0) = -2.83669;
//   unit(1,1) = 10.0125;
//   unit(1,2) = -3.9293;
//   unit(1,3) = 2.95569;
//   unit(1,4) = 3.82807;
  
//   unit(2,0) = 1.34247;
//   unit(2,1) = -2.47754;
//   unit(2,2) = 9.67165;
//   unit(2,3) = -2.10167;
//   unit(2,4) = -0.829575;
  
//   unit(3,0) = -0.991873;
//   unit(3,1) = 1.82176;
//   unit(3,2) = -2.10167;
//   unit(3,3) = 7.51674;
//   unit(3,4) = 0.64796; 
 
//   unit(4,0) = -4.76361;
//   unit(4,1) = 6.3756;
//   unit(4,2) = -6.89025;
//   unit(4,3) = 5.3251;
//   unit(4,4) = 17.9026; 

//   unit.LUdecomp(pivot);
 
//   ExpLUResult(0,0) = -4.76361000000000 ;
//   ExpLUResult(0,1) = 6.37560000000000;
//   ExpLUResult(0,2) = -6.890250000000000;
//   ExpLUResult(0,3) = 5.325100000000000;
//   ExpLUResult(0,4) = 17.902600000000000;

//   ExpLUResult(1,0) =-0.816107951742481;
//   ExpLUResult(1,1) =6.89471785712936;
//   ExpLUResult(1,2) =-7.69712781449363;
//   ExpLUResult(1,3) =5.96575645382389;
//   ExpLUResult(1,4) =19.1125342168649;
  
//   ExpLUResult(2,0) =0.595491654438546;
//   ExpLUResult(2,1) =0.901542824052500;
//   ExpLUResult(2,2) =7.11307671896682;
//   ExpLUResult(2,3) =-5.59374753004051;
//   ExpLUResult(2,4) =-24.0635469654240;
  
//   ExpLUResult(3,0) = 0.208218766859588;
//   ExpLUResult(3,1) = 0.0716839238749926;
//   ExpLUResult(3,2) = -0.0161997855341865;
//   ExpLUResult(3,3) = 5.88968790278431;
//   ExpLUResult(3,4) = -4.83958304347208;

//   ExpLUResult(4,0) = -0.281817781052605;
//   ExpLUResult(4,1) = -0.0987397264439268;
//   ExpLUResult(4,2) = 0.979863287669927;
//   ExpLUResult(4,3) = 0.928606568229258;
//   ExpLUResult(4,4) = 34.1759172514949;

//   double small = 1.e-13;

//   EXPECT_NEAR(ExpLUResult(0,0),unit(pivot(0),0),small);
//   EXPECT_NEAR(ExpLUResult(0,1),unit(pivot(0),1),small);
//   EXPECT_NEAR(ExpLUResult(0,2),unit(pivot(0),2),small);
//   EXPECT_NEAR(ExpLUResult(0,3),unit(pivot(0),3),small);
//   EXPECT_NEAR(ExpLUResult(0,4),unit(pivot(0),4),small);
 
//   EXPECT_NEAR(ExpLUResult(1,0),unit(pivot(1),0),small);
//   EXPECT_NEAR(ExpLUResult(1,1),unit(pivot(1),1),small);
//   EXPECT_NEAR(ExpLUResult(1,2),unit(pivot(1),2),small);
//   EXPECT_NEAR(ExpLUResult(1,3),unit(pivot(1),3),small);
//   EXPECT_NEAR(ExpLUResult(1,4),unit(pivot(1),4),small);
   
//   EXPECT_NEAR(ExpLUResult(2,0),unit(pivot(2),0),small);
//   EXPECT_NEAR(ExpLUResult(2,1),unit(pivot(2),1),small);
//   EXPECT_NEAR(ExpLUResult(2,2),unit(pivot(2),2),small);
//   EXPECT_NEAR(ExpLUResult(2,3),unit(pivot(2),3),small);
//   EXPECT_NEAR(ExpLUResult(2,4),unit(pivot(2),4),small);

//   EXPECT_NEAR(ExpLUResult(3,0),unit(pivot(3),0),small);
//   EXPECT_NEAR(ExpLUResult(3,1),unit(pivot(3),1),small);
//   EXPECT_NEAR(ExpLUResult(3,2),unit(pivot(3),2),small);
//   EXPECT_NEAR(ExpLUResult(3,3),unit(pivot(3),3),small);
//   EXPECT_NEAR(ExpLUResult(3,4),unit(pivot(3),4),small);
  
//   EXPECT_NEAR(ExpLUResult(4,0),unit(pivot(4),0),small);
//   EXPECT_NEAR(ExpLUResult(4,1),unit(pivot(4),1),small);
//   EXPECT_NEAR(ExpLUResult(4,2),unit(pivot(4),2),small);
//   EXPECT_NEAR(ExpLUResult(4,3),unit(pivot(4),3),small);
//   EXPECT_NEAR(ExpLUResult(4,4),unit(pivot(4),4),small);
  
// }

TEST(SquareMatrix, inv3){

 Array1D<int> pivot(5);
 SquareMatrix<double> unit(5);
 Array1D<double> b(5), x(5);
 
  b(0) = 5.24556;
  b(1) = 3.971;
  b(2) = 4.12354;
  b(3) = .6892154;
  b(4) = 1.6195;
  
  unit(0,0) = 3.88762;
  unit(0,1) = 1.69154;
  unit(0,2) = -2.07394;
  unit(0,3) = 1.6199;
  unit(0,4) = 4.50208;
  
  unit(1,0) = -2.83669;
  unit(1,1) = 10.0125;
  unit(1,2) = -3.9293;
  unit(1,3) = 2.95569;
  unit(1,4) = 3.82807;
  
  unit(2,0) = 1.34247;
  unit(2,1) = -2.47754;
  unit(2,2) = 9.67165;
  unit(2,3) = -2.10167;
  unit(2,4) = -0.829575;
  
  unit(3,0) = -0.991873;
  unit(3,1) = 1.82176;
  unit(3,2) = -2.10167;
  unit(3,3) = 7.51674;
  unit(3,4) = 0.64796; 
 
  unit(4,0) = -4.76361;
  unit(4,1) = 6.3756;
  unit(4,2) = -6.89025;
  unit(4,3) = 5.3251;
  unit(4,4) = 17.9026; 

  unit.squareSolve(b,x);
  
  double small = 1.e-13;
  
  EXPECT_NEAR(0.965360116727089,x(0),small);
  EXPECT_NEAR(0.739468484520957,x(1),small);
  EXPECT_NEAR(0.539374998411889,x(2),small);
  EXPECT_NEAR(0.169887158672791,x(3),small);
  EXPECT_NEAR(0.241043552885589,x(4),small);

}


TEST(SquareMatrix, surreal_eval_LU) {

  SquareMatrix<Surreal<double,MATRIX_N> > unit(MATRIX_N);
  Surreal<double,MATRIX_N> q[MATRIX_N];
  q[0] = 1.0;
  q[1] = 2.0;
  q[2] = 3.0;
  q[3] = 4.0;
  q[0].Deriv(0) = 1.0;
  q[1].Deriv(1) = 1.0;
  q[2].Deriv(2) = 1.0;
  q[3].Deriv(3) = 1.0;

  unit(0,0) = q[0]*q[1];
  unit(0,1) = q[1]*q[2];
  unit(0,2) = q[3]*q[0];
  unit(0,3) = q[2]/q[1];

  unit(1,0) = q[1]*q[1];
  unit(1,1) = q[3]*q[2];
  unit(1,2) = q[2]*q[0];
  unit(1,3) = q[0]/q[1];

  unit(2,0) = q[1]*q[3];
  unit(2,1) = q[3]*q[1];
  unit(2,2) = q[2]*q[2];
  unit(2,3) = q[0]/q[0];

  unit(3,0) = q[1]*q[0];
  unit(3,1) = q[2]*q[1];
  unit(3,2) = q[3]*q[2];
  unit(3,3) = q[0]/q[3];

  Array1D<int> pivot(MATRIX_N);
  unit.LUdecomp(pivot);
 

  EXPECT_EQ(2,pivot(0));
  EXPECT_EQ(1,pivot(1));
  EXPECT_EQ(3,pivot(2));
  EXPECT_EQ(0,pivot(3));
 
  Array2D<double> ExpLUResult(MATRIX_N,MATRIX_N);
  ExpLUResult(0,0) = 8.;
  ExpLUResult(0,1) = 8.;
  ExpLUResult(0,2) = 9.;
  ExpLUResult(0,3) = 1.;

  ExpLUResult(1,0) = 0.5;
  ExpLUResult(1,1) = 8.0;
  ExpLUResult(1,2) = -1.5;
  ExpLUResult(1,3) = 0.0;

  ExpLUResult(2,0) = 0.25;
  ExpLUResult(2,1) = 0.5;
  ExpLUResult(2,2) = 10.5;
  ExpLUResult(2,3) = 0.0;

  ExpLUResult(3,0) = 0.25;
  ExpLUResult(3,1) = 0.5;
  ExpLUResult(3,2) = 0.238095238095238;
  ExpLUResult(3,3) = 1.25;
  
  double small = 1.e-13;

  EXPECT_NEAR(ExpLUResult(0,0),unit(pivot(0),0).Value(),small);
  EXPECT_NEAR(ExpLUResult(0,1),unit(pivot(0),1).Value(),small);
  EXPECT_NEAR(ExpLUResult(0,2),unit(pivot(0),2).Value(),small);
  EXPECT_NEAR(ExpLUResult(0,3),unit(pivot(0),3).Value(),small);
   
  EXPECT_NEAR(ExpLUResult(1,0),unit(pivot(1),0).Value(),small);
  EXPECT_NEAR(ExpLUResult(1,1),unit(pivot(1),1).Value(),small);
  EXPECT_NEAR(ExpLUResult(1,2),unit(pivot(1),2).Value(),small);
  EXPECT_NEAR(ExpLUResult(1,3),unit(pivot(1),3).Value(),small);
  
  EXPECT_NEAR(ExpLUResult(2,0),unit(pivot(2),0).Value(),small);
  EXPECT_NEAR(ExpLUResult(2,1),unit(pivot(2),1).Value(),small);
  EXPECT_NEAR(ExpLUResult(2,2),unit(pivot(2),2).Value(),small);
  EXPECT_NEAR(ExpLUResult(2,3),unit(pivot(2),3).Value(),small);

  EXPECT_NEAR(ExpLUResult(3,0),unit(pivot(3),0).Value(),small);
  EXPECT_NEAR(ExpLUResult(3,1),unit(pivot(3),1).Value(),small);
  EXPECT_NEAR(ExpLUResult(3,2),unit(pivot(3),2).Value(),small);
  EXPECT_NEAR(ExpLUResult(3,3),unit(pivot(3),3).Value(),small);
}

TEST(SquareMatrix, surreal_eval_inv) {

  SquareMatrix<Surreal<double,MATRIX_N> > unit(MATRIX_N);
  Surreal<double,MATRIX_N> q[MATRIX_N];
  q[0] = 1.0;
  q[1] = 2.0;
  q[2] = 3.0;
  q[3] = 4.0;
  q[0].Deriv(0) = 1.0;
  q[1].Deriv(1) = 1.0;
  q[2].Deriv(2) = 1.0;
  q[3].Deriv(3) = 1.0;

  unit(0,0) = q[0]*q[1];
  unit(0,1) = q[1]*q[2];
  unit(0,2) = q[3]*q[0];
  unit(0,3) = q[2]/q[1];

  unit(1,0) = q[1]*q[1];
  unit(1,1) = q[3]*q[2];
  unit(1,2) = q[2]*q[0];
  unit(1,3) = q[0]/q[1];

  unit(2,0) = q[1]*q[3];
  unit(2,1) = q[3]*q[1];
  unit(2,2) = q[2]*q[2];
  unit(2,3) = q[0]/q[0];

  unit(3,0) = q[1]*q[0];
  unit(3,1) = q[2]*q[1];
  unit(3,2) = q[3]*q[2];
  unit(3,3) = q[0]/q[3];
  
  unit.invert();

  Array2D<double> ExpINVResult(MATRIX_N,MATRIX_N);
  ExpINVResult(0,0) = -0.100000000000000;
  ExpINVResult(0,1) = -2.440476190476192e-002;
  ExpINVResult(0,2) = 0.187500000000000;
  ExpINVResult(0,3) = -0.101190476190476;

  ExpINVResult(1,0) =  0.000000000000000e+000;
  ExpINVResult(1,1) = 0.116071428571429;
  ExpINVResult(1,2) = -6.250000000000000e-002;
  ExpINVResult(1,3) = 1.785714285714286e-002;

  ExpINVResult(2,0) = 0.000000000000000e+000;
  ExpINVResult(2,1) = -4.761904761904762e-002;
  ExpINVResult(2,2) = 0.000000000000000e+000;
  ExpINVResult(2,3) = 9.523809523809523e-002;

  ExpINVResult(3,0) = 0.800000000000000;
  ExpINVResult(3,1) = -0.304761904761905;
  ExpINVResult(3,2) = 0.000000000000000e+000;
  ExpINVResult(3,3) = -0.190476190476190;
  
  double small = 1.e-13;

  EXPECT_NEAR(ExpINVResult(0,0),unit(0,0).Value(),small);
  EXPECT_NEAR(ExpINVResult(0,1),unit(0,1).Value(),small);
  EXPECT_NEAR(ExpINVResult(0,2),unit(0,2).Value(),small);
  EXPECT_NEAR(ExpINVResult(0,3),unit(0,3).Value(),small);
   
  EXPECT_NEAR(ExpINVResult(1,0),unit(1,0).Value(),small);
  EXPECT_NEAR(ExpINVResult(1,1),unit(1,1).Value(),small);
  EXPECT_NEAR(ExpINVResult(1,2),unit(1,2).Value(),small);
  EXPECT_NEAR(ExpINVResult(1,3),unit(1,3).Value(),small);
  
  EXPECT_NEAR(ExpINVResult(2,0),unit(2,0).Value(),small);
  EXPECT_NEAR(ExpINVResult(2,1),unit(2,1).Value(),small);
  EXPECT_NEAR(ExpINVResult(2,2),unit(2,2).Value(),small);
  EXPECT_NEAR(ExpINVResult(2,3),unit(2,3).Value(),small);

  EXPECT_NEAR(ExpINVResult(3,0),unit(3,0).Value(),small);
  EXPECT_NEAR(ExpINVResult(3,1),unit(3,1).Value(),small);
  EXPECT_NEAR(ExpINVResult(3,2),unit(3,2).Value(),small);
  EXPECT_NEAR(ExpINVResult(3,3),unit(3,3).Value(),small);



  Array2D<double> ExpINVResultderiv0(MATRIX_N,MATRIX_N);
  ExpINVResultderiv0(0,0) = -1.023809523809524E-002;
  ExpINVResultderiv0(0,1) = -4.378259637188209E-002;
  ExpINVResultderiv0(0,2) = 7.544642857142858E-002;
  ExpINVResultderiv0(0,3) = -2.791950113378688E-003;

  ExpINVResultderiv0(1,0) = -4.642857142857143E-002; 
  ExpINVResultderiv0(1,1) = 3.650085034013605E-002;
  ExpINVResultderiv0(1,2) = -6.696428571428573E-003;
  ExpINVResultderiv0(1,3) = -1.764455782312925E-002;

  ExpINVResultderiv0(2,0) = 1.904761904761905E-002;
  ExpINVResultderiv0(2,1) = -2.154195011337868E-003;
  ExpINVResultderiv0(2,2) = -3.571428571428572E-002;
  ExpINVResultderiv0(2,3) = 3.287981859410431E-002;

  ExpINVResultderiv0(3,0) = 0.281904761904762;
  ExpINVResultderiv0(3,1) = 7.764172335600909E-002;
  ExpINVResultderiv0(3,2) = -0.228571428571429;
  ExpINVResultderiv0(3,3) = -0.132426303854875;
  

  EXPECT_NEAR(ExpINVResultderiv0(0,0),unit(0,0).Deriv(0),small);
  EXPECT_NEAR(ExpINVResultderiv0(0,1),unit(0,1).Deriv(0),small);
  EXPECT_NEAR(ExpINVResultderiv0(0,2),unit(0,2).Deriv(0),small);
  EXPECT_NEAR(ExpINVResultderiv0(0,3),unit(0,3).Deriv(0),small);
   
  EXPECT_NEAR(ExpINVResultderiv0(1,0),unit(1,0).Deriv(0),small);
  EXPECT_NEAR(ExpINVResultderiv0(1,1),unit(1,1).Deriv(0),small);
  EXPECT_NEAR(ExpINVResultderiv0(1,2),unit(1,2).Deriv(0),small);
  EXPECT_NEAR(ExpINVResultderiv0(1,3),unit(1,3).Deriv(0),small);
  
  EXPECT_NEAR(ExpINVResultderiv0(2,0),unit(2,0).Deriv(0),small);
  EXPECT_NEAR(ExpINVResultderiv0(2,1),unit(2,1).Deriv(0),small);
  EXPECT_NEAR(ExpINVResultderiv0(2,2),unit(2,2).Deriv(0),small);
  EXPECT_NEAR(ExpINVResultderiv0(2,3),unit(2,3).Deriv(0),small);

  EXPECT_NEAR(ExpINVResultderiv0(3,0),unit(3,0).Deriv(0),small);
  EXPECT_NEAR(ExpINVResultderiv0(3,1),unit(3,1).Deriv(0),small);
  EXPECT_NEAR(ExpINVResultderiv0(3,2),unit(3,2).Deriv(0),small);
  EXPECT_NEAR(ExpINVResultderiv0(3,3),unit(3,3).Deriv(0),small);


  Array2D<double> ExpINVResultderiv1(MATRIX_N,MATRIX_N);
  ExpINVResultderiv0(0,0) = -1.976190476190476E-002;
  ExpINVResultderiv0(0,1) = 1.873157596371882E-002;
  ExpINVResultderiv0(0,2) = -7.544642857142858E-002;
  ExpINVResultderiv0(0,3) = 5.848922902494331E-002;

  ExpINVResultderiv0(1,0) = 4.642857142857143E-002; 
  ExpINVResultderiv0(1,1) = 1.962159863945578E-002;
  ExpINVResultderiv0(1,2) = -5.580357142857142E-002;
  ExpINVResultderiv0(1,3) = 2.147108843537415E-002;

  ExpINVResultderiv0(2,0) = -1.904761904761905E-002;
  ExpINVResultderiv0(2,1) = -3.185941043083901E-002;
  ExpINVResultderiv0(2,2) = 3.571428571428572E-002;
  ExpINVResultderiv0(2,3) = -1.247165532879818E-002;

  ExpINVResultderiv0(3,0) = 0.358095238095238;
  ExpINVResultderiv0(3,1) = -0.386757369614512;
  ExpINVResultderiv0(3,2) = 0.228571428571429;
  ExpINVResultderiv0(3,3) = -0.194104308390023;
  

  EXPECT_NEAR(ExpINVResultderiv0(0,0),unit(0,0).Deriv(1),small);
  EXPECT_NEAR(ExpINVResultderiv0(0,1),unit(0,1).Deriv(1),small);
  EXPECT_NEAR(ExpINVResultderiv0(0,2),unit(0,2).Deriv(1),small);
  EXPECT_NEAR(ExpINVResultderiv0(0,3),unit(0,3).Deriv(1),small);
   
  EXPECT_NEAR(ExpINVResultderiv0(1,0),unit(1,0).Deriv(1),small);
  EXPECT_NEAR(ExpINVResultderiv0(1,1),unit(1,1).Deriv(1),small);
  EXPECT_NEAR(ExpINVResultderiv0(1,2),unit(1,2).Deriv(1),small);
  EXPECT_NEAR(ExpINVResultderiv0(1,3),unit(1,3).Deriv(1),small);
  
  EXPECT_NEAR(ExpINVResultderiv0(2,0),unit(2,0).Deriv(1),small);
  EXPECT_NEAR(ExpINVResultderiv0(2,1),unit(2,1).Deriv(1),small);
  EXPECT_NEAR(ExpINVResultderiv0(2,2),unit(2,2).Deriv(1),small);
  EXPECT_NEAR(ExpINVResultderiv0(2,3),unit(2,3).Deriv(1),small);

  EXPECT_NEAR(ExpINVResultderiv0(3,0),unit(3,0).Deriv(1),small);
  EXPECT_NEAR(ExpINVResultderiv0(3,1),unit(3,1).Deriv(1),small);
  EXPECT_NEAR(ExpINVResultderiv0(3,2),unit(3,2).Deriv(1),small);
  EXPECT_NEAR(ExpINVResultderiv0(3,3),unit(3,3).Deriv(1),small);



  Array2D<double> ExpINVResultderiv2(MATRIX_N,MATRIX_N);
  ExpINVResultderiv0(0,0) = 4.000000000000000E-002;
  ExpINVResultderiv0(0,1) = 7.593253968253967E-002;
  ExpINVResultderiv0(0,2) = -3.125000000000000E-002;
  ExpINVResultderiv0(0,3) = -6.686507936507935E-002;

  ExpINVResultderiv0(1,0) = 0.000000000000000E+000; 
  ExpINVResultderiv0(1,1) = -6.696428571428571E-002;
  ExpINVResultderiv0(1,2) = 3.125000000000000E-002;
  ExpINVResultderiv0(1,3) = 8.928571428571428E-003;

  ExpINVResultderiv0(2,0) = 0.000000000000000E+000;
  ExpINVResultderiv0(2,1) = 1.587301587301588E-002;
  ExpINVResultderiv0(2,2) = 0.000000000000000E+000;
  ExpINVResultderiv0(2,3) = -3.174603174603175E-002;

  ExpINVResultderiv0(3,0) = -0.320000000000000;
  ExpINVResultderiv0(3,1) = 7.111111111111111E-002;
  ExpINVResultderiv0(3,2) = 0.000000000000000E+000;
  ExpINVResultderiv0(3,3) = 0.177777777777778;
  

  EXPECT_NEAR(ExpINVResultderiv0(0,0),unit(0,0).Deriv(2),small);
  EXPECT_NEAR(ExpINVResultderiv0(0,1),unit(0,1).Deriv(2),small);
  EXPECT_NEAR(ExpINVResultderiv0(0,2),unit(0,2).Deriv(2),small);
  EXPECT_NEAR(ExpINVResultderiv0(0,3),unit(0,3).Deriv(2),small);
   
  EXPECT_NEAR(ExpINVResultderiv0(1,0),unit(1,0).Deriv(2),small);
  EXPECT_NEAR(ExpINVResultderiv0(1,1),unit(1,1).Deriv(2),small);
  EXPECT_NEAR(ExpINVResultderiv0(1,2),unit(1,2).Deriv(2),small);
  EXPECT_NEAR(ExpINVResultderiv0(1,3),unit(1,3).Deriv(2),small);
  
  EXPECT_NEAR(ExpINVResultderiv0(2,0),unit(2,0).Deriv(2),small);
  EXPECT_NEAR(ExpINVResultderiv0(2,1),unit(2,1).Deriv(2),small);
  EXPECT_NEAR(ExpINVResultderiv0(2,2),unit(2,2).Deriv(2),small);
  EXPECT_NEAR(ExpINVResultderiv0(2,3),unit(2,3).Deriv(2),small);

  EXPECT_NEAR(ExpINVResultderiv0(3,0),unit(3,0).Deriv(2),small);
  EXPECT_NEAR(ExpINVResultderiv0(3,1),unit(3,1).Deriv(2),small);
  EXPECT_NEAR(ExpINVResultderiv0(3,2),unit(3,2).Deriv(2),small);
  EXPECT_NEAR(ExpINVResultderiv0(3,3),unit(3,3).Deriv(2),small);


  Array2D<double> ExpINVResultderiv3(MATRIX_N,MATRIX_N);
  ExpINVResultderiv0(0,0) = 3.244047619047619E-002;
  ExpINVResultderiv0(0,1) = -4.316716269841270E-002;
  ExpINVResultderiv0(0,2) = -5.145089285714286E-002;
  ExpINVResultderiv0(0,3) = 7.219742063492063E-002;

  ExpINVResultderiv0(1,0) = -1.160714285714286E-002; 
  ExpINVResultderiv0(1,1) = -2.674851190476191E-002;
  ExpINVResultderiv0(1,2) = 3.738839285714286E-002;
  ExpINVResultderiv0(1,3) = -2.194940476190476E-002;

  ExpINVResultderiv0(2,0) = 4.761904761904762E-003;
  ExpINVResultderiv0(2,1) = 2.837301587301587E-002;
  ExpINVResultderiv0(2,2) = -8.928571428571430E-003;
  ExpINVResultderiv0(2,3) = -2.579365079365080E-002;

  ExpINVResultderiv0(3,0) = -9.523809523809523E-003;
  ExpINVResultderiv0(3,1) = 0.120634920634921;
  ExpINVResultderiv0(3,2) = -5.714285714285715E-002;
  ExpINVResultderiv0(3,3) = -3.174603174603162E-003;
  

  EXPECT_NEAR(ExpINVResultderiv0(0,0),unit(0,0).Deriv(3),small);
  EXPECT_NEAR(ExpINVResultderiv0(0,1),unit(0,1).Deriv(3),small);
  EXPECT_NEAR(ExpINVResultderiv0(0,2),unit(0,2).Deriv(3),small);
  EXPECT_NEAR(ExpINVResultderiv0(0,3),unit(0,3).Deriv(3),small);
   
  EXPECT_NEAR(ExpINVResultderiv0(1,0),unit(1,0).Deriv(3),small);
  EXPECT_NEAR(ExpINVResultderiv0(1,1),unit(1,1).Deriv(3),small);
  EXPECT_NEAR(ExpINVResultderiv0(1,2),unit(1,2).Deriv(3),small);
  EXPECT_NEAR(ExpINVResultderiv0(1,3),unit(1,3).Deriv(3),small);
  
  EXPECT_NEAR(ExpINVResultderiv0(2,0),unit(2,0).Deriv(3),small);
  EXPECT_NEAR(ExpINVResultderiv0(2,1),unit(2,1).Deriv(3),small);
  EXPECT_NEAR(ExpINVResultderiv0(2,2),unit(2,2).Deriv(3),small);
  EXPECT_NEAR(ExpINVResultderiv0(2,3),unit(2,3).Deriv(3),small);

  EXPECT_NEAR(ExpINVResultderiv0(3,0),unit(3,0).Deriv(3),small);
  EXPECT_NEAR(ExpINVResultderiv0(3,1),unit(3,1).Deriv(3),small);
  EXPECT_NEAR(ExpINVResultderiv0(3,2),unit(3,2).Deriv(3),small);
  EXPECT_NEAR(ExpINVResultderiv0(3,3),unit(3,3).Deriv(3),small);
}


TEST(SquareMatrix, surreal_eval_solve) {

  SquareMatrix<Surreal<double,MATRIX_N> > unit(MATRIX_N);
  Surreal<double,MATRIX_N> q[MATRIX_N];
  q[0] = 1.0;
  q[1] = 2.0;
  q[2] = 3.0;
  q[3] = 4.0;
  q[0].Deriv(0) = 1.0;
  q[1].Deriv(1) = 1.0;
  q[2].Deriv(2) = 1.0;
  q[3].Deriv(3) = 1.0;
  
  unit(0,0) = q[0]*q[1];
  unit(0,1) = q[1]*q[2];
  unit(0,2) = q[3]*q[0];
  unit(0,3) = q[2]/q[1];

  unit(1,0) = q[1]*q[1];
  unit(1,1) = q[3]*q[2];
  unit(1,2) = q[2]*q[0];
  unit(1,3) = q[0]/q[1];

  unit(2,0) = q[1]*q[3];
  unit(2,1) = q[3]*q[1];
  unit(2,2) = q[2]*q[2];
  unit(2,3) = q[0]/q[0];

  unit(3,0) = q[1]*q[0];
  unit(3,1) = q[2]*q[1];
  unit(3,2) = q[3]*q[2];
  unit(3,3) = q[0]/q[3];

  Array1D<Surreal<double,MATRIX_N> > x(MATRIX_N);
  Array1D<Surreal<double,MATRIX_N> > b(MATRIX_N);

  b(0) = q[0]*q[3];
  b(1) = q[1]*q[2];
  b(2) = q[3]*q[0];
  b(3) = q[2]*q[1];

  unit.squareSolve(b,x);

  Array1D<double> ExpsolveResult(MATRIX_N);
  ExpsolveResult(0) = -0.403571428571429;
  ExpsolveResult(1) = 0.553571428571429;
  ExpsolveResult(2) = 0.285714285714286;
  ExpsolveResult(3) = 0.228571428571429;

  double small = 1.e-13;

  EXPECT_NEAR(ExpsolveResult(0),x(0).Value(),small);
  EXPECT_NEAR(ExpsolveResult(1),x(1).Value(),small);
  EXPECT_NEAR(ExpsolveResult(2),x(2).Value(),small);
  EXPECT_NEAR(ExpsolveResult(3),x(3).Value(),small);

  EXPECT_NEAR(0.331386054421769,x(0).Deriv(0),small);
  EXPECT_NEAR(-0.349362244897959,x(1).Deriv(0),small);
  EXPECT_NEAR(0.117687074829932,x(2).Deriv(0),small);
  EXPECT_NEAR(3.08462585034014,x(3).Deriv(0),small);
   
}
