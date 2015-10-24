#include <cmath>
#include "gtest/gtest.h"
#include "Surreal/Reverse/SurrealR.h"

//------------------------------ SURREAL POW -----------------------------------
//--->This tests the SurrealPow class of the Surreal package
//-->TEST 1.0: pow(Surreal, Surreal)
TEST(SurrealPow, SurrealRaisedToSurreal) {
  SurrealR<double, 3> a,b,c;
  a.SetValue(2.0);
  a.SetDeriv(0,2.0);
  a.SetDeriv(1,3.0);
  a.SetDeriv(2,4.0);

  b.SetValue(3.0);
  b.SetDeriv(0,3.0);
  b.SetDeriv(1,4.0);
  b.SetDeriv(2,5.0);

  c = pow(a,b);
  double y = std::pow(2.0,3.0);
  double dy1 = y*(b.GetDeriv(0)*std::log(a.GetValue()) +
                  b.GetValue()/a.GetValue()*a.GetDeriv(0));
  double dy2 = y*(b.GetDeriv(1)*std::log(a.GetValue()) +
                  b.GetValue()/a.GetValue()*a.GetDeriv(1));
  double dy3 = y*(b.GetDeriv(2)*std::log(a.GetValue()) +
                  b.GetValue()/a.GetValue()*a.GetDeriv(2));

  //test value
  EXPECT_DOUBLE_EQ(std::pow(2.0,3.0), c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(dy1, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(dy2, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(dy3, c.GetDeriv(2));
}

//-->TEST 2.0: pow(double, Surreal)
TEST(SurrealPow, DoubleRaisedToSurreal) {
  SurrealR<double, 3> b,c;
  double x = 1.1;

  b.SetValue(3.0);
  b.SetDeriv(0,3.0);
  b.SetDeriv(1,4.0);
  b.SetDeriv(2,5.0);

  c = pow(x,b);

  //test value
  EXPECT_DOUBLE_EQ(std::pow(1.1,3.0), c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(std::pow(1.1,3.0)*std::log(x) * 3.0, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(std::pow(1.1,3.0)*std::log(x) * 4.0, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(std::pow(1.1,3.0)*std::log(x) * 5.0, c.GetDeriv(2));
}
//-->TEST 3.0: pow(Surreal lvalue, double)
TEST(SurrealPow, SurrealRaisedToDouble) {
  SurrealR<double, 3> a,c;
  a.SetValue(2.0);
  a.SetDeriv(0,2.0);
  a.SetDeriv(1,3.0);
  a.SetDeriv(2,4.0);
  double x = 1.1;

  c = pow(a,x);

  //test value
  EXPECT_DOUBLE_EQ(std::pow(2.0,1.1), c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(1.1 * std::pow(2.0,1.1-1.0) * 2.0, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(1.1 * std::pow(2.0,1.1-1.0) * 3.0, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(1.1 * std::pow(2.0,1.1-1.0) * 4.0, c.GetDeriv(2));
}
