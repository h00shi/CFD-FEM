#include "gtest/gtest.h"
#include "Surreal/Surreal.h"

//------------------------------ SURREAL POW -----------------------------------
//--->This tests the SurrealPow class of the Surreal package
//-->TEST 1: pow(Surreal, Surreal)
TEST(SurrealPow, SurrealRaisedToSurreal) {
  Surreal<double, 3> a,b,c;
  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  b.Value()  = 3.0;
  b.Deriv(0) = 3.0;
  b.Deriv(1) = 4.0;
  b.Deriv(2) = 5.0;

  c = pow(a,b);
  double y = std::pow(2.0,3.0);
  double dy1 = y*(b.Deriv(0)*std::log(a.Value()) +
                  b.Value()/a.Value()*a.Deriv(0));
  double dy2 = y*(b.Deriv(1)*std::log(a.Value()) +
                  b.Value()/a.Value()*a.Deriv(1));
  double dy3 = y*(b.Deriv(2)*std::log(a.Value()) +
                  b.Value()/a.Value()*a.Deriv(2));

  //test value
  EXPECT_DOUBLE_EQ(std::pow(2.0,3.0), c.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(dy1, c.Deriv(0));
  EXPECT_DOUBLE_EQ(dy2, c.Deriv(1));
  EXPECT_DOUBLE_EQ(dy3, c.Deriv(2));
}

//-->TEST 2: pow(double, Surreal)
TEST(SurrealPow, DoubleRaisedToSurreal) {
  Surreal<double, 3> b,c;
  double x = 1.1;

  b.Value()  = 3.0;
  b.Deriv(0) = 3.0;
  b.Deriv(1) = 4.0;
  b.Deriv(2) = 5.0;

  c = pow(x,b);
  //test value
  EXPECT_DOUBLE_EQ(std::pow(1.1,3.0), c.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(3.0*std::log(x) * 3.0, c.Deriv(0));
  EXPECT_DOUBLE_EQ(3.0*std::log(x) * 4.0, c.Deriv(1));
  EXPECT_DOUBLE_EQ(3.0*std::log(x) * 5.0, c.Deriv(2));
}

//-->TEST 3: pow(Surreal, double)
TEST(SurrealPow, SurrealRaisedToDouble) {
  Surreal<double, 3> a,c;
  double x = 1.1;

  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  c = pow(a,x);
  //test value
  EXPECT_DOUBLE_EQ(std::pow(2.0,1.1), c.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(1.1 * std::pow(2.0,1.1-1.0) * 2.0, c.Deriv(0));
  EXPECT_DOUBLE_EQ(1.1 * std::pow(2.0,1.1-1.0) * 3.0, c.Deriv(1));
  EXPECT_DOUBLE_EQ(1.1 * std::pow(2.0,1.1-1.0) * 4.0, c.Deriv(2));
}
