#include "gtest/gtest.h"
#include "Surreal/Forward/Surreal.h"

//------------------------------ SURREAL MULTIPLY ------------------------------
//--->This tests the SurrealMultiply class of the Surreal package
//-->TEST 1: Surreal * Surreal
TEST(SurrealMultiply, SurrealTimesSurreal) {
  Surreal<double, 3> a,b,c;
  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  b.Value()  = 3.0;
  b.Deriv(0) = 3.0;
  b.Deriv(1) = 4.0;
  b.Deriv(2) = 5.0;

  c = a * b;
  //test value
  EXPECT_DOUBLE_EQ(6.0, c.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(12.0, c.Deriv(0));
  EXPECT_DOUBLE_EQ(17.0, c.Deriv(1));
  EXPECT_DOUBLE_EQ(22.0, c.Deriv(2));
}

//-->TEST 2: double * Surreal
TEST(SurrealMultiply, DoubleTimesSurreal) {
  Surreal<double, 3> b,c;
  double x = 1.1;
  b.Value()  = 3.0;
  b.Deriv(0) = 3.0;
  b.Deriv(1) = 4.0;
  b.Deriv(2) = 5.0;

  c = x * b;
  //test value
  EXPECT_DOUBLE_EQ(3.3, c.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(3.3, c.Deriv(0));
  EXPECT_DOUBLE_EQ(4.4, c.Deriv(1));
  EXPECT_DOUBLE_EQ(5.5, c.Deriv(2));
}

//-->TEST 3: Surreal * double
TEST(SurrealMultiply, SurrealTimesDouble) {
  Surreal<double, 3> a,c;
  double x = 1.1;
  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  c = a * x;
  //test value
  EXPECT_DOUBLE_EQ(2.2, c.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(2.2, c.Deriv(0));
  EXPECT_DOUBLE_EQ(3.3, c.Deriv(1));
  EXPECT_DOUBLE_EQ(4.4, c.Deriv(2));
}
