#include "gtest/gtest.h"
#include "Surreal/Reverse/SurrealR.h"

//------------------------------ SURREAL MULTIPLY ------------------------------
//--->This tests the SurrealMultiply class of the Surreal package
//-->TEST 1.0: Surreal lvalue * Surreal lvalue
TEST(SurrealMultiply, SurrealLvalueTimesSurrealLvalue) {
  SurrealR<double, 3> a,b,c;
  a.SetValue(2.0);
  a.SetDeriv(0,2.0);
  a.SetDeriv(1,3.0);
  a.SetDeriv(2,4.0);

  b.SetValue(3.0);
  b.SetDeriv(0,3.0);
  b.SetDeriv(1,4.0);
  b.SetDeriv(2,5.0);

  c = a * b;
  //test value
  EXPECT_DOUBLE_EQ(6.0, c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(12.0, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(17.0, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(22.0, c.GetDeriv(2));
}

//-->TEST 1.1: Surreal rvalue * Surreal rvalue
TEST(SurrealMultiply, SurrealRvalueTimesSurrealRvalue) {
  SurrealR<double, 3> a,b,c;
  double x = 1.1;
  a.SetValue(2.0);
  a.SetDeriv(0,2.0);
  a.SetDeriv(1,3.0);
  a.SetDeriv(2,4.0);

  b.SetValue(3.0);
  b.SetDeriv(0,3.0);
  b.SetDeriv(1,4.0);
  b.SetDeriv(2,5.0);

  c = (a*x) * (b*x);
  //test value
  EXPECT_DOUBLE_EQ(7.26, c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(14.52, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(20.57, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(26.62, c.GetDeriv(2));
}

//-->TEST 1.2: Surreal lvalue * Surreal rvalue
TEST(SurrealMultiply, SurrealLvalueTimesSurrealRvalue) {
  SurrealR<double, 3> a,b,c;
  double x = 1.1;
  a.SetValue(2.0);
  a.SetDeriv(0,2.0);
  a.SetDeriv(1,3.0);
  a.SetDeriv(2,4.0);

  b.SetValue(3.0);
  b.SetDeriv(0,3.0);
  b.SetDeriv(1,4.0);
  b.SetDeriv(2,5.0);

  c = a * (b*x);
  //test value
  EXPECT_DOUBLE_EQ(6.6, c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(13.20, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(18.70, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(24.20, c.GetDeriv(2));
}

//-->TEST 1.3: Surreal rvalue * Surreal lvalue
TEST(SurrealMultiply, SurrealRvalueTimesSurrealLvalue) {
  SurrealR<double, 3> a,b,c;
  double x = 1.1;
  a.SetValue(2.0);
  a.SetDeriv(0,2.0);
  a.SetDeriv(1,3.0);
  a.SetDeriv(2,4.0);

  b.SetValue(3.0);
  b.SetDeriv(0,3.0);
  b.SetDeriv(1,4.0);
  b.SetDeriv(2,5.0);

  c = (a*x) * b;
  //test value
  EXPECT_DOUBLE_EQ(6.6, c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(13.20, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(18.70, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(24.20, c.GetDeriv(2));
}


//-->TEST 2.0: double * Surreal lvalue
TEST(SurrealMultiply, DoubleTimesSurrealLvalue) {
  SurrealR<double, 3> b,c;
  double x = 1.1;
  b.SetValue(3.0);
  b.SetDeriv(0, 3.0);
  b.SetDeriv(1, 4.0);
  b.SetDeriv(2, 5.0);

  c = x * b;
  //test value
  EXPECT_DOUBLE_EQ(3.3, c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(3.3, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(4.4, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(5.5, c.GetDeriv(2));
}

//-->TEST 2.1: double * Surreal rvalue
TEST(SurrealMultiply, DoubleTimesSurrealRvalue) {
  SurrealR<double, 3> b,c;
  double x = 1.1;
  b.SetValue(3.0);
  b.SetDeriv(0, 3.0);
  b.SetDeriv(1, 4.0);
  b.SetDeriv(2, 5.0);

  c = x * (b*x);
  //test value
  EXPECT_DOUBLE_EQ(3.63, c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(3.63, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(4.84, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(6.05, c.GetDeriv(2));
}

//-->TEST 3.0: Surreal lvalue * double
TEST(SurrealMultiply, SurrealLvalueTimesDouble) {
  SurrealR<double, 3> a,c;
  double x = 1.1;
  a.SetValue(2.0);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);

  c = a * x;
  //test value
  EXPECT_DOUBLE_EQ(2.2, c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(2.2, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(3.3, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(4.4, c.GetDeriv(2));
}

//-->TEST 3.1: Surreal rvalue * double
TEST(SurrealMultiply, SurrealRvalueTimesDouble) {
  SurrealR<double, 3> a,c;
  double x = 1.1;
  a.SetValue(2.0);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);

  c = (a*x) * x;
  //test value
  EXPECT_DOUBLE_EQ(2.42, c.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(2.42, c.GetDeriv(0));
  EXPECT_DOUBLE_EQ(3.63, c.GetDeriv(1));
  EXPECT_DOUBLE_EQ(4.84, c.GetDeriv(2));
}
