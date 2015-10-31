#include "Surreal/Reverse/SurrealR.h"
#include "gtest/gtest.h"
#include <cmath>
//------------------------------ SURREAL OPERATOR+------------------------------
//--->This tests the SurrealPos class of the Surreal package
//-->TEST 1.0: +(Surreal)
TEST(SurrealPos, PositiveOfSurreal) {
  SurrealR<double, 3> a, c;

  a.SetValue(2.0);
  a.SetDeriv(0,2.0);
  a.SetDeriv(1,3.0);
  a.SetDeriv(2,4.0);

  c = +(a);
  //test value
  EXPECT_DOUBLE_EQ(2.0, c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(2.0, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(3.0, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(4.0, c.GetDeriv(2));
}

//------------------------------ SURREAL OPERATOR- -----------------------------
//--->This tests the SurrealNeg class of the Surreal package
TEST(SurrealNeg, NegationOfSurreal) {
  SurrealR<double, 3> a, c;

  a.SetValue(2.0);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);

  c = -(a);
  //test value
  EXPECT_DOUBLE_EQ(-2.0, c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(-2.0, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(-3.0, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(-4.0, c.GetDeriv(2));
}

//------------------------------ SURREAL LOG -----------------------------------
//--->This tests the SurrealLog class of the Surreal package
TEST(SurrealLog, LogOfSurreal) {
  SurrealR<double, 3> a, c;

  a.SetValue(2.0);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);

  c = log(a);
  //test value
  EXPECT_DOUBLE_EQ(std::log(2.0), c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ( 1.0/2.0 * 2.0, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ( 1.0/2.0 * 3.0, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ( 1.0/2.0 * 4.0, c.GetDeriv(2));
}

//------------------------------ SURREAL LOG10 -----------------------------------
//--->This tests the SurrealLog10 class of the Surreal package
TEST(SurrealLog10, Log10OfSurreal) {
  SurrealR<double, 3> a, c;

  a.SetValue(100.0);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);

  c = log10(a);
  //test value
  EXPECT_DOUBLE_EQ(std::log10(100.0), c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ( 1.0/100.0 /log(10.0)* 2.0, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ( 1.0/100.0 /log(10.0)* 3.0, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ( 1.0/100.0 /log(10.0)* 4.0, c.GetDeriv(2));

}

//------------------------------ SURREAL SQRT-----------------------------------
//--->This tests the SurrealSqrt class of the Surreal package
TEST(SurrealSqrt, SqrtOfSurreal) {
  SurrealR<double, 3> a, c;

  a.SetValue(2.0);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);

  c = sqrt(a);
  //test value
  EXPECT_DOUBLE_EQ(std::sqrt(2.0), c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ( 1.0/2.0 * std::pow(2.0,-1.0/2.0) * 2.0, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ( 1.0/2.0 * std::pow(2.0,-1.0/2.0) * 3.0, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ( 1.0/2.0 * std::pow(2.0,-1.0/2.0) * 4.0, c.GetDeriv(2));
}

//------------------------------ SURREAL EXP-----------------------------------
//--->This tests the SurrealExp class of the Surreal package
TEST(SurrealExp, ExpofSurreal) {
  SurrealR<double, 3> a, c;

  a.SetValue(2.0);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);

  c = exp(a);
  //test value
  EXPECT_DOUBLE_EQ(std::exp(2.0), c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(std::exp(2.0) * 2.0, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(std::exp(2.0) * 3.0, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(std::exp(2.0) * 4.0, c.GetDeriv(2));
}

//------------------------------ SURREAL ABS -----------------------------------
//--->This tests the SurrealAbs class of the Surreal package
TEST(SurrealAbs, AbsofSurreal) {
  SurrealR<double, 3> a, c;

  a.SetValue(2.0);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);

  c = abs(a);
  //test value
  EXPECT_DOUBLE_EQ(std::abs(2.0), c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(2.0, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(3.0, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(4.0, c.GetDeriv(2));

  a.SetValue(-2.0);
  c = abs(a);
  //test value
  EXPECT_DOUBLE_EQ(std::abs(2.0), c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(-2.0, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(-3.0, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(-4.0, c.GetDeriv(2));
}

//------------------------------ SURREAL FABS ----------------------------------
//--->This tests the SurrealFabs class of the Surreal package
TEST(SurrealFabs, FabsofSurreal) {
  SurrealR<double, 3> a, c;

  a.SetValue(2.);
  a.SetDeriv(0, 2.);
  a.SetDeriv(1, 3.);
  a.SetDeriv(2, 4.);

  c = fabs(a);
  //test value
  EXPECT_DOUBLE_EQ(std::fabs(2.0), c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(2.0, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(3.0, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(4.0, c.GetDeriv(2));

  a.SetValue(-2.);
  c = fabs(a);
  //test value
  EXPECT_DOUBLE_EQ(std::fabs(2.0), c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(-2.0, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(-3.0, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(-4.0, c.GetDeriv(2));
}


//------------------------------ SURREAL TRIG ----------------------------------
//--->This tests the SurrealTrig class of the Surreal package
//--> TESTING COS FUNCTION
TEST(SurrealTrig, Cos) {
  SurrealR<double, 3> a, c;

  a.SetValue(2.0);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);

  c = cos(a);

  //test value
  EXPECT_DOUBLE_EQ(std::cos(2.0), c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(-std::sin(2.0)*2.0, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(-std::sin(2.0)*3.0, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(-std::sin(2.0)*4.0, c.GetDeriv(2));
}

//--> TESTING SIN FUNCTION
TEST(SurrealTrig, Sin) {
  SurrealR<double, 3> a, c;

  a.SetValue(2.0);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);

  c = sin(a);

  //test value
  EXPECT_DOUBLE_EQ(std::sin(2.0), c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(std::cos(2.0)*2.0, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(std::cos(2.0)*3.0, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(std::cos(2.0)*4.0, c.GetDeriv(2));
}

//--> TESTING TAN FUNCTION
TEST(SurrealTrig, Tan) {
  SurrealR<double, 3> a, c;

  a.SetValue(2.0);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);

  c = tan(a);

  //test value
  EXPECT_DOUBLE_EQ(std::tan(2.0), c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(1.0 /(std::cos(2.0)*std::cos(2.0)) *2.0, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(1.0 /(std::cos(2.0)*std::cos(2.0)) *3.0, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(1.0 /(std::cos(2.0)*std::cos(2.0)) *4.0, c.GetDeriv(2));
}

//--> TESTING ACOS FUNCTION
TEST(SurrealTrig, Acos) {
  SurrealR<double, 3> a, c;

  a.SetValue(0.2);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);

  c = acos(a);

  //test value
  EXPECT_DOUBLE_EQ(std::acos(0.2), c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(-1.0 /std::sqrt(1.0 - std::pow(0.2,2)) *2.0, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(-1.0 /std::sqrt(1.0 - std::pow(0.2,2)) *3.0, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(-1.0 /std::sqrt(1.0 - std::pow(0.2,2)) *4.0, c.GetDeriv(2));
}

//--> TESTING ASIN FUNCTION
TEST(SurrealTrig, Asin) {
  SurrealR<double, 3> a, c;

  a.SetValue(0.2);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);

  c = asin(a);

  //test value
  EXPECT_DOUBLE_EQ(std::asin(0.2), c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(1.0 /std::sqrt(1.0 - std::pow(0.2,2)) *2.0, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(1.0 /std::sqrt(1.0 - std::pow(0.2,2)) *3.0, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(1.0 /std::sqrt(1.0 - std::pow(0.2,2)) *4.0, c.GetDeriv(2));
}

//--> TESTING ATAN FUNCTION
TEST(SurrealTrig, Atan) {
  SurrealR<double, 3> a, c;

  a.SetValue(0.2);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);

  c = atan(a);

  //test value
  EXPECT_DOUBLE_EQ(std::atan(0.2), c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(1.0 /(1.0 + std::pow(0.2,2)) *2.0, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(1.0 /(1.0 + std::pow(0.2,2)) *3.0, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(1.0 /(1.0 + std::pow(0.2,2)) *4.0, c.GetDeriv(2));
}