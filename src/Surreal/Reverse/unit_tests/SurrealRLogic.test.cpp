#include "gtest/gtest.h"
#include "Surreal/Reverse/SurrealR.h"

//----------------------------- SURREAL LOGIC ----------------------------------
//--->These tests test the SurrealLogic class of the Surreal package
//--> TESTING < OPERATOR
//->testing surreal < surreal
TEST(SurrealLogicLessThan, SurrealVsSurreal) {
  SurrealR<double, 3> a, b;
  double eps = 1.0e-6;

  a.SetValue(2.0);

  //test when values are equal
  b.SetValue(a.GetValue());
  EXPECT_FALSE(a < b);

  //test when rhs is greater
  b.SetValue(a.GetValue() + eps);
  EXPECT_TRUE(a < b);

  //test when rhs is smaller
  b.SetValue(a.GetValue() - eps);
  EXPECT_FALSE(a < b);
}

//->testing double < surreal
TEST(SurrealLogicLessThan, DoubleVsSurreal){
  SurrealR<double, 3> a;
  double b;
  double eps = 1.0e-6;

  a.SetValue(2.0);

  //test when values are equal
  b = a.GetValue();
  EXPECT_FALSE(a < b);

  //test when rhs is greater
  b = a.GetValue() + eps;
  EXPECT_TRUE(a < b);

  //test when rhs is smaller
  b = a.GetValue() - eps;
  EXPECT_FALSE(a < b);
}

//->testing surreal < double
TEST(SurrealLogicLessThan, SurrealVsDouble){
  SurrealR<double, 3> b;
  double a;
  double eps = 1.0e-6;

  b.SetValue(2.0);

  //test when values are equal
  a = b.GetValue();
  EXPECT_FALSE(a < b);

  //test when rhs is greater
  a = b.GetValue() - eps;
  EXPECT_TRUE(a < b);

  //test when rhs is smaller
  a = b.GetValue() + eps;
  EXPECT_FALSE(a < b);
}

//--> TESTING <= OPERATOR
//->testing surreal <= surreal
TEST(SurrealLogicEqualOrLessThan, SurrealVsSurreal) {
  SurrealR<double, 3> a, b;
  double eps = 1.0e-6;

  a.SetValue(2.0);

  //test when values are equal
  b.SetValue(a.GetValue());
  EXPECT_TRUE(a <= b);

  //test when rhs is greater
  b.SetValue(a.GetValue() + eps);
  EXPECT_TRUE(a <= b);

  //test when rhs is smaller
  b.SetValue(a.GetValue() - eps);
  EXPECT_FALSE(a <= b);
}

//->testing double <= surreal
TEST(SurrealLogicEqualOrLessThan, DoubleVsSurreal) {

  SurrealR<double, 3> a;
  double b;
  double eps = 1.0e-6;

  a.SetValue(2.0);

  //test when values are equal
  b = a.GetValue();
  EXPECT_TRUE(a <= b);

  //test when rhs is greater
  b = a.GetValue() + eps;
  EXPECT_TRUE(a <= b);

  //test when rhs is smaller
  b = a.GetValue() - eps;
  EXPECT_FALSE(a <= b);
}

//->testing surreal <= double
TEST(SurrealLogicEqualOrLessThan, SurrealVsDouble){

  SurrealR<double, 3> b;
  double a;
  double eps = 1.0e-6;

  b.SetValue(2.0);

  //test when values are equal
  a = b.GetValue();
  EXPECT_TRUE(a <= b);

  //test when rhs is greater
  a = b.GetValue() - eps;
  EXPECT_TRUE(a <= b);

  //test when rhs is smaller
  a = b.GetValue() + eps;
  EXPECT_FALSE(a <= b);
}

//--> TESTING > OPERATOR
//->testing surreal > surreal
TEST(SurrealLogicGreaterThan, SurrealVsSurreal) {
  SurrealR<double, 3> a, b;
  double eps = 1.0e-6;

  a.SetValue(2.0);

  //test when values are equal
  b.SetValue(a.GetValue());
  EXPECT_FALSE(a > b);

  //test when rhs is greater
  b.SetValue(a.GetValue() + eps);
  EXPECT_FALSE(a > b);

  //test when rhs is smaller
  b.SetValue(a.GetValue() - eps);
  EXPECT_TRUE(a > b);
}

//->testing double > surreal
TEST(SurrealLogicGreaterThan, DoubleVsSurreal) {

  SurrealR<double, 3> a;
  double b;
  double eps = 1.0e-6;

  a.SetValue(2.0);

  //test when values are equal
  b = a.GetValue();
  EXPECT_FALSE(a > b);

  //test when rhs is greater
  b = a.GetValue() + eps;
  EXPECT_FALSE(a > b);

  //test when rhs is smaller
  b = a.GetValue() - eps;
  EXPECT_TRUE(a > b);
}

//->testing surreal > double
TEST(SurrealLogicGreaterThan, SurrealVsDouble){

  SurrealR<double, 3> b;
  double a;
  double eps = 1.0e-6;

  b.SetValue(2.0);

  //test when values are equal
  a = b.GetValue();
  EXPECT_FALSE(a > b);

  //test when rhs is greater
  a = b.GetValue() - eps;
  EXPECT_FALSE(a > b);

  //test when rhs is smaller
  a = b.GetValue() + eps;
  EXPECT_TRUE(a > b);
}

//--> TESTING >= OPERATOR
//->testing surreal >= surreal
TEST(SurrealLogicEqualOrGreaterThan, SurrealVsSurreal) {
  SurrealR<double, 3> a, b;
  double eps = 1.0e-6;

  a.SetValue(2.0);

  //test when values are equal
  b.SetValue(a.GetValue());
  EXPECT_TRUE(a >= b);

  //test when rhs is greater
  b.SetValue(a.GetValue() + eps);
  EXPECT_FALSE(a >= b);

  //test when rhs is smaller
  b.SetValue(a.GetValue() - eps);
  EXPECT_TRUE(a >= b);
}

//->testing double >= surreal
TEST(SurrealLogicEqualOrGreaterThan, DoubleVsSurreal) {
  SurrealR<double, 3> a;
  double b;
  double eps = 1.0e-6;

  a.SetValue(2.0);

  //test when values are equal
  b = a.GetValue();
  EXPECT_TRUE(a >= b);

  //test when rhs is greater
  b = a.GetValue() + eps;
  EXPECT_FALSE(a >= b);

  //test when rhs is smaller
  b = a.GetValue() - eps;
  EXPECT_TRUE(a >= b);
}

//->testing surreal >= double
TEST(SurrealLogicEqualOrGreaterThan, SurrealVsDouble){
  SurrealR<double, 3> b;
  double a;
  double eps = 1.0e-6;
  b.SetValue(2.0);

  //test when values are equal
  a = b.GetValue();
  EXPECT_TRUE(a >= b);

  //test when rhs is greater
  a = b.GetValue() - eps;
  EXPECT_FALSE(a >= b);

  //test when rhs is smaller
  a = b.GetValue() + eps;
  EXPECT_TRUE(a >= b);
}
