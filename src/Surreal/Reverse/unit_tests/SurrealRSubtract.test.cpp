#include "gtest/gtest.h"
#include "Surreal/Reverse/SurrealR.h"

//------------------------------ SURREAL SUBTRACT ------------------------------
//--->This tests the SurrealSubtract class of the Surreal package
//-->TEST 1.0: Surreal - Surreal
TEST(SurrealSubtract, SurrealMinusSurreal) {
  SurrealR<double, 3> a,b,c;
  a.SetValue(2.0);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);

  b.SetValue(3.0);
  b.SetDeriv(0, 3.0);
  b.SetDeriv(1, 5.0);
  b.SetDeriv(2, 7.0);

  c = a - b;
  //test value
  EXPECT_DOUBLE_EQ(-1.0, c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(-1.0, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(-2.0, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(-3.0, c.GetDeriv(2));
}

//-->TEST 2.0: double - Surreal
TEST(SurrealSubtract, DoubleMinusSurreal) {
  SurrealR<double, 3> b,c;
  double x = 1.1;

  b.SetValue(3.0);
  b.SetDeriv(0, 3.0);
  b.SetDeriv(1, 5.0);
  b.SetDeriv(2, 7.0);

  c = x - b;
  //test value
  EXPECT_DOUBLE_EQ(-1.9, c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(-3.0, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(-5.0, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(-7.0, c.GetDeriv(2));
}

//-->TEST 3.0: Surreal - double
TEST(SurrealSubtract, SurrealMinusDouble) {
  SurrealR<double, 3> a,c;
  double x = 1.1;
  a.SetValue(2.0);
  a.SetDeriv(0,2.0);
  a.SetDeriv(1,3.0);
  a.SetDeriv(2,4.0);

  c = a - x;
  //test value
  EXPECT_DOUBLE_EQ(0.9, c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(2.0, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(3.0, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(4.0, c.GetDeriv(2));
}
