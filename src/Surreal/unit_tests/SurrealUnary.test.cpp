#include "gtest/gtest.h"
#include "Surreal.h"

//------------------------------ SURREAL OPERATOR+------------------------------
//--->This tests the SurrealPos class of the Surreal package
TEST(SurrealPos, PositiveOfSurreal) {
  Surreal<double, 3> a, c;

  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  c = +(a);
  //test value
  EXPECT_DOUBLE_EQ(2.0, c.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(2.0, c.Deriv(0));
  EXPECT_DOUBLE_EQ(3.0, c.Deriv(1));
  EXPECT_DOUBLE_EQ(4.0, c.Deriv(2));
}

//------------------------------ SURREAL OPERATOR- -----------------------------
//--->This tests the SurrealNeg class of the Surreal package
TEST(SurrealNeg, NegationOfSurreal) {
  Surreal<double, 3> a, c;

  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  c = -(a);
  //test value
  EXPECT_DOUBLE_EQ(-2.0, c.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(-2.0, c.Deriv(0));
  EXPECT_DOUBLE_EQ(-3.0, c.Deriv(1));
  EXPECT_DOUBLE_EQ(-4.0, c.Deriv(2));
}

//------------------------------ SURREAL LOG -----------------------------------
//--->This tests the SurrealLog class of the Surreal package
TEST(SurrealLog, LogOfSurreal) {
  Surreal<double, 3> a, c;

  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  c = log(a);
  //test value
  EXPECT_DOUBLE_EQ(std::log(2.0), c.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ( 1.0/2.0 * 2.0, c.Deriv(0));
  EXPECT_DOUBLE_EQ( 1.0/2.0 * 3.0, c.Deriv(1));
  EXPECT_DOUBLE_EQ( 1.0/2.0 * 4.0, c.Deriv(2));
}


//------------------------------ SURREAL SQRT-----------------------------------
//--->This tests the SurrealSqrt class of the Surreal package
TEST(SurrealSqrt, SqrtOfSurreal) {
  Surreal<double, 3> a, c;

  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  c = sqrt(a);
  //test value
  EXPECT_DOUBLE_EQ(std::sqrt(2.0), c.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ( 1.0/2.0 * std::pow(2.0,-1.0/2.0) * 2.0, c.Deriv(0));
  EXPECT_DOUBLE_EQ( 1.0/2.0 * std::pow(2.0,-1.0/2.0) * 3.0, c.Deriv(1));
  EXPECT_DOUBLE_EQ( 1.0/2.0 * std::pow(2.0,-1.0/2.0) * 4.0, c.Deriv(2));
}

//------------------------------ SURREAL EXP-----------------------------------
//--->This tests the SurrealExp class of the Surreal package
TEST(SurrealExp, ExpofSurreal) {
  Surreal<double, 3> a, c;

  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  c = exp(a);
  //test value
  EXPECT_DOUBLE_EQ(std::exp(2.0), c.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(std::exp(2.0) * 2.0, c.Deriv(0));
  EXPECT_DOUBLE_EQ(std::exp(2.0) * 3.0, c.Deriv(1));
  EXPECT_DOUBLE_EQ(std::exp(2.0) * 4.0, c.Deriv(2));
}

//------------------------------ SURREAL ABS -----------------------------------
//--->This tests the SurrealExp class of the Surreal package
TEST(SurrealAbs, AbsofSurreal) {
  Surreal<double, 3> a, c;

  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  c = abs(a);
  //test value
  EXPECT_DOUBLE_EQ(std::abs(2.0), c.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(2.0, c.Deriv(0));
  EXPECT_DOUBLE_EQ(3.0, c.Deriv(1));
  EXPECT_DOUBLE_EQ(4.0, c.Deriv(2));

  a.Value()  = -2.0;
  c = abs(a);
  //test value
  EXPECT_DOUBLE_EQ(std::abs(2.0), c.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(-2.0, c.Deriv(0));
  EXPECT_DOUBLE_EQ(-3.0, c.Deriv(1));
  EXPECT_DOUBLE_EQ(-4.0, c.Deriv(2));
}

//------------------------------ SURREAL TRIG ----------------------------------
//--->This tests the SurrealTrig class of the Surreal package
//--> TESTING COS FUNCTION
TEST(SurrealTrig, Cos) {
  Surreal<double, 3> a, c;

  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  c = cos(a);

  //test value
  EXPECT_DOUBLE_EQ(std::cos(2.0), c.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(-std::sin(2.0)*2.0, c.Deriv(0));
  EXPECT_DOUBLE_EQ(-std::sin(2.0)*3.0, c.Deriv(1));
  EXPECT_DOUBLE_EQ(-std::sin(2.0)*4.0, c.Deriv(2));
}
//--> TESTING SIN FUNCTION
TEST(SurrealTrig, Sin) {
  Surreal<double, 3> a, c;

  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  c = sin(a);

  //test value
  EXPECT_DOUBLE_EQ(std::sin(2.0), c.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(std::cos(2.0)*2.0, c.Deriv(0));
  EXPECT_DOUBLE_EQ(std::cos(2.0)*3.0, c.Deriv(1));
  EXPECT_DOUBLE_EQ(std::cos(2.0)*4.0, c.Deriv(2));
}
//--> TESTING TAN FUNCTION
TEST(SurrealTrig, Tan) {
  Surreal<double, 3> a, c;

  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  c = tan(a);

  //test value
  EXPECT_DOUBLE_EQ(std::tan(2.0), c.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(1.0 /(std::cos(2.0)*std::cos(2.0)) *2.0, c.Deriv(0));
  EXPECT_DOUBLE_EQ(1.0 /(std::cos(2.0)*std::cos(2.0)) *3.0, c.Deriv(1));
  EXPECT_DOUBLE_EQ(1.0 /(std::cos(2.0)*std::cos(2.0)) *4.0, c.Deriv(2));
}
//--> TESTING ACOS FUNCTION
TEST(SurrealTrig, Acos) {
  Surreal<double, 3> a, c;

  a.Value()  = 0.2;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  c = acos(a);

  //test value
  EXPECT_DOUBLE_EQ(std::acos(0.2), c.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(-1.0 /std::sqrt(1.0 - std::pow(0.2,2)) *2.0, c.Deriv(0));
  EXPECT_DOUBLE_EQ(-1.0 /std::sqrt(1.0 - std::pow(0.2,2)) *3.0, c.Deriv(1));
  EXPECT_DOUBLE_EQ(-1.0 /std::sqrt(1.0 - std::pow(0.2,2)) *4.0, c.Deriv(2));
}
//--> TESTING ASIN FUNCTION
TEST(SurrealTrig, Asin) {
  Surreal<double, 3> a, c;

  a.Value()  = 0.2;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  c = asin(a);

  //test value
  EXPECT_DOUBLE_EQ(std::asin(0.2), c.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(1.0 /std::sqrt(1.0 - std::pow(0.2,2)) *2.0, c.Deriv(0));
  EXPECT_DOUBLE_EQ(1.0 /std::sqrt(1.0 - std::pow(0.2,2)) *3.0, c.Deriv(1));
  EXPECT_DOUBLE_EQ(1.0 /std::sqrt(1.0 - std::pow(0.2,2)) *4.0, c.Deriv(2));
}
//--> TESTING ATAN FUNCTION
TEST(SurrealTrig, Atan) {
  Surreal<double, 3> a, c;

  a.Value()  = 0.2;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  c = atan(a);

  //test value
  EXPECT_DOUBLE_EQ(std::atan(0.2), c.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(1.0 /(1.0 + std::pow(0.2,2)) *2.0, c.Deriv(0));
  EXPECT_DOUBLE_EQ(1.0 /(1.0 + std::pow(0.2,2)) *3.0, c.Deriv(1));
  EXPECT_DOUBLE_EQ(1.0 /(1.0 + std::pow(0.2,2)) *4.0, c.Deriv(2));
}
