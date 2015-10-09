#include "gtest/gtest.h"
#include "Surreal/Surreal.h"

//------------------------------ SURREAL ADD -----------------------------------
//--->This tests the SurrealAdd class of the Surreal package
//-->TEST 1: Surreal + Surreal
TEST(SurrealAdd, SurrealPlusSurreal) {
  Surreal<double, 3> a,b,c;

  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  b.Value()  = 3.0;
  b.Deriv(0) = 3.0;
  b.Deriv(1) = 4.0;
  b.Deriv(2) = 5.0;

  c = a + b;
  //test value
  EXPECT_DOUBLE_EQ(5.0, c.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(5.0, c.Deriv(0));
  EXPECT_DOUBLE_EQ(7.0, c.Deriv(1));
  EXPECT_DOUBLE_EQ(9.0, c.Deriv(2));
}

//-->TEST 2: double  + Surreal
TEST(SurrealAdd, DoublePlusSurreal) {
  Surreal<double, 3> b,c;
  double x = 1.1;

  b.Value()  = 3.0;
  b.Deriv(0) = 3.0;
  b.Deriv(1) = 4.0;
  b.Deriv(2) = 5.0;

  c = x + b;
  //test value
  EXPECT_DOUBLE_EQ(4.1, c.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(3.0, c.Deriv(0));
  EXPECT_DOUBLE_EQ(4.0, c.Deriv(1));
  EXPECT_DOUBLE_EQ(5.0, c.Deriv(2));
}

//-->TEST 3: Surreal + double
TEST(SurrealAdd, SurrealPlusDouble) {
  Surreal<double, 3> a,c;
  double x = 1.1;

  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  c = a + x;
  //test value
  EXPECT_DOUBLE_EQ(3.1, c.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(2.0, c.Deriv(0));
  EXPECT_DOUBLE_EQ(3.0, c.Deriv(1));
  EXPECT_DOUBLE_EQ(4.0, c.Deriv(2));
}
