#include "gtest/gtest.h"
#include "Surreal.h"

//----------------------------- SURREAL LOGIC ----------------------------------
//--->These tests test the SurrealLogic class of the Surreal package
//--> TESTING < OPERATOR
//->testing surreal < surreal
TEST(SurrealLogicLessThan, SurrealVsSurreal) {
  Surreal<double, 3> a, b;
  double eps = 1.0e-6;

  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  //test when values are equal
  b.Value()  = a.Value();
  EXPECT_FALSE(a < b);

  //test when rhs is greater
  b.Value()  = a.Value() + eps;
  EXPECT_TRUE(a < b);

  //test when rhs is smaller
  b.Value()  = a.Value() - eps;
  EXPECT_FALSE(a < b);
}

//->testing double < surreal
TEST(SurrealLogicLessThan, DoubleVsSurreal){
  Surreal<double, 3> a;
  double b;
  double eps = 1.0e-6;

  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  //test when values are equal
  b = a.Value();
  EXPECT_FALSE(a < b);

  //test when rhs is greater
  b = a.Value() + eps;
  EXPECT_TRUE(a < b);

  //test when rhs is smaller
  b = a.Value() - eps;
  EXPECT_FALSE(a < b);
}

//->testing surreal < double
TEST(SurrealLogicLessThan, SurrealVsDouble){
  Surreal<double, 3> b;
  double a;
  double eps = 1.0e-6;

  b.Value()  = 2.0;
  b.Deriv(0) = 2.0;
  b.Deriv(1) = 3.0;
  b.Deriv(2) = 4.0;

  //test when values are equal
  a = b.Value();
  EXPECT_FALSE(a < b);

  //test when rhs is greater
  a = b.Value() - eps;
  EXPECT_TRUE(a < b);

  //test when rhs is smaller
  a = b.Value() + eps;
  EXPECT_FALSE(a < b);
}

//--> TESTING <= OPERATOR
//->testing surreal <= surreal
TEST(SurrealLogicEqualOrLessThan, SurrealVsSurreal) {
  Surreal<double, 3> a, b;
  double eps = 1.0e-6;

  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  //test when values are equal
  b.Value()  = a.Value();
  EXPECT_TRUE(a <= b);

  //test when rhs is greater
  b.Value()  = a.Value() + eps;
  EXPECT_TRUE(a <= b);

  //test when rhs is smaller
  b.Value()  = a.Value() - eps;
  EXPECT_FALSE(a <= b);
}

//->testing double <= surreal
TEST(SurrealLogicEqualOrLessThan, DoubleVsSurreal) {

  Surreal<double, 3> a;
  double b;
  double eps = 1.0e-6;

  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  //test when values are equal
  b = a.Value();
  EXPECT_TRUE(a <= b);

  //test when rhs is greater
  b = a.Value() + eps;
  EXPECT_TRUE(a <= b);

  //test when rhs is smaller
  b = a.Value() - eps;
  EXPECT_FALSE(a <= b);
}

//->testing surreal <= double
TEST(SurrealLogicEqualOrLessThan, SurrealVsDouble){

  Surreal<double, 3> b;
  double a;
  double eps = 1.0e-6;

  b.Value()  = 2.0;
  b.Deriv(0) = 2.0;
  b.Deriv(1) = 3.0;
  b.Deriv(2) = 4.0;

  //test when values are equal
  a = b.Value();
  EXPECT_TRUE(a <= b);

  //test when rhs is greater
  a = b.Value() - eps;
  EXPECT_TRUE(a <= b);

  //test when rhs is smaller
  a = b.Value() + eps;
  EXPECT_FALSE(a <= b);
}

//--> TESTING > OPERATOR
//->testing surreal > surreal
TEST(SurrealLogicGreaterThan, SurrealVsSurreal) {
  Surreal<double, 3> a, b;
  double eps = 1.0e-6;

  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  //test when values are equal
  b.Value()  = a.Value();
  EXPECT_FALSE(a > b);

  //test when rhs is greater
  b.Value()  = a.Value() + eps;
  EXPECT_FALSE(a > b);

  //test when rhs is smaller
  b.Value()  = a.Value() - eps;
  EXPECT_TRUE(a > b);
}

//->testing double > surreal
TEST(SurrealLogicGreaterThan, DoubleVsSurreal) {

  Surreal<double, 3> a;
  double b;
  double eps = 1.0e-6;

  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  //test when values are equal
  b = a.Value();
  EXPECT_FALSE(a > b);

  //test when rhs is greater
  b = a.Value() + eps;
  EXPECT_FALSE(a > b);

  //test when rhs is smaller
  b = a.Value() - eps;
  EXPECT_TRUE(a > b);
}

//->testing surreal > double
TEST(SurrealLogicGreaterThan, SurrealVsDouble){

  Surreal<double, 3> b;
  double a;
  double eps = 1.0e-6;

  b.Value()  = 2.0;
  b.Deriv(0) = 2.0;
  b.Deriv(1) = 3.0;
  b.Deriv(2) = 4.0;

  //test when values are equal
  a = b.Value();
  EXPECT_FALSE(a > b);

  //test when rhs is greater
  a = b.Value() - eps;
  EXPECT_FALSE(a > b);

  //test when rhs is smaller
  a = b.Value() + eps;
  EXPECT_TRUE(a > b);
}

//--> TESTING >= OPERATOR
//->testing surreal >= surreal
TEST(SurrealLogicEqualOrGreaterThan, SurrealVsSurreal) {
  Surreal<double, 3> a, b;
  double eps = 1.0e-6;

  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  //test when values are equal
  b.Value()  = a.Value();
  EXPECT_TRUE(a >= b);

  //test when rhs is greater
  b.Value()  = a.Value() + eps;
  EXPECT_FALSE(a >= b);

  //test when rhs is smaller
  b.Value()  = a.Value() - eps;
  EXPECT_TRUE(a >= b);
}

//->testing double >= surreal
TEST(SurrealLogicEqualOrGreaterThan, DoubleVsSurreal) {
  Surreal<double, 3> a;
  double b;
  double eps = 1.0e-6;


  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  //test when values are equal
  b = a.Value();
  EXPECT_TRUE(a >= b);

  //test when rhs is greater
  b = a.Value() + eps;
  EXPECT_FALSE(a >= b);

  //test when rhs is smaller
  b = a.Value() - eps;
  EXPECT_TRUE(a >= b);
}

//->testing surreal >= double
TEST(SurrealLogicEqualOrGreaterThan, SurrealVsDouble){
  Surreal<double, 3> b;
  double a;
  double eps = 1.0e-6;

  b.Value()  = 2.0;
  b.Deriv(0) = 2.0;
  b.Deriv(1) = 3.0;
  b.Deriv(2) = 4.0;

  //test when values are equal
  a = b.Value();
  EXPECT_TRUE(a >= b);

  //test when rhs is greater
  a = b.Value() - eps;
  EXPECT_FALSE(a >= b);

  //test when rhs is smaller
  a = b.Value() + eps;
  EXPECT_TRUE(a >= b);
}

// //--> TESTING == OPERATOR
// //->testing surreal == surreal
// TEST(SurrealLogicEqualTo, SurrealVsSurreal) {
//   Surreal<double, 3> a, b;
//   double eps = 1.0e-6;

//   a.Value()  = 2.0;
//   a.Deriv(0) = 2.0;
//   a.Deriv(1) = 3.0;
//   a.Deriv(2) = 4.0;

//   //test when values are equal
//   b.Value()  = a.Value();
//   EXPECT_TRUE(a == b);

//   //test when rhs is greater
//   b.Value()  = a.Value() + eps;
//   EXPECT_FALSE(a == b);

//   //test when rhs is smaller
//   b.Value()  = a.Value() - eps;
//   EXPECT_FALSE(a == b);
// }

// //->testing double == surreal
// TEST(SurrealLogicEqualTo, DoubleVsSurreal) {
//   Surreal<double, 3> a;
//   double b;
//   double eps = 1.0e-6;


//   a.Value()  = 2.0;
//   a.Deriv(0) = 2.0;
//   a.Deriv(1) = 3.0;
//   a.Deriv(2) = 4.0;

//   //test when values are equal
//   b = a.Value();
//   EXPECT_TRUE(a == b);

//   //test when rhs is greater
//   b = a.Value() + eps;
//   EXPECT_FALSE(a == b);

//   //test when rhs is smaller
//   b = a.Value() - eps;
//   EXPECT_FALSE(a == b);
// }

// //->testing surreal == double
// TEST(SurrealLogicEqualTo, SurrealVsDouble){
//   Surreal<double, 3> b;
//   double a;
//   double eps = 1.0e-6;

//   b.Value()  = 2.0;
//   b.Deriv(0) = 2.0;
//   b.Deriv(1) = 3.0;
//   b.Deriv(2) = 4.0;

//   //test when values are equal
//   a = b.Value();
//   EXPECT_TRUE(a == b);

//   //test when rhs is greater
//   a = b.Value() - eps;
//   EXPECT_FALSE(a == b);

//   //test when rhs is smaller
//   a = b.Value() + eps;
//   EXPECT_FALSE(a == b);
// }
