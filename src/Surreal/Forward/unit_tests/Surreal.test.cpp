#include "gtest/gtest.h"
#include "Surreal/Forward/Surreal.h"

//------------------------ CONSTRUCTORS ----------------------------------------
TEST(SurrealConstructors, DefaultConstructor) {
  Surreal<double, 3> a;

  //test value
  EXPECT_DOUBLE_EQ(0.0, a.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(0.0, a.Deriv(0));
  EXPECT_DOUBLE_EQ(0.0, a.Deriv(1));
  EXPECT_DOUBLE_EQ(0.0, a.Deriv(2));
}

TEST(SurrealConstructors, ConstructorFromRealNumber) {
  double x = 1.1;
  Surreal<double, 3> a = x;

  //test value
  EXPECT_DOUBLE_EQ(1.1, a.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(0.0, a.Deriv(0));
  EXPECT_DOUBLE_EQ(0.0, a.Deriv(1));
  EXPECT_DOUBLE_EQ(0.0, a.Deriv(2));
}

TEST(SurrealConstructors, CopyConstructor) {
  Surreal<double, 3> a;
  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  Surreal<double, 3> b = a;

  //test value
  EXPECT_DOUBLE_EQ(2.0, b.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(2.0, b.Deriv(0));
  EXPECT_DOUBLE_EQ(3.0, b.Deriv(1));
  EXPECT_DOUBLE_EQ(4.0, b.Deriv(2));
}

//------------------------ SURREAL ASSIGNMENT ----------------------------------
//--> TESTING = real
TEST(SurrealAssignment, AssignToReal) {
  Surreal<double, 3> a;
  double x = 1.1;

  a = x;
  //test value
  EXPECT_DOUBLE_EQ(1.1, a.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(0.0, a.Deriv(0));
  EXPECT_DOUBLE_EQ(0.0, a.Deriv(1));
  EXPECT_DOUBLE_EQ(0.0, a.Deriv(2));
}

//--> TESTING = Surreal
TEST(SurrealAssignment, AssignToSurreal) {
  Surreal<double, 3> a,b;
  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  b = a;

  //test value
  EXPECT_DOUBLE_EQ(2.0, b.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(2.0, b.Deriv(0));
  EXPECT_DOUBLE_EQ(3.0, b.Deriv(1));
  EXPECT_DOUBLE_EQ(4.0, b.Deriv(2));
}

//------------------------ SURREAL COMPOUND ASSIGNMENT -------------------------
//--> TESTING +=
TEST(SurrealCompoundAssignment, AddSurreal) {
  //--->This tests the += feature of the Surreal class
  Surreal<double, 3> a,b;
  Surreal<double,2> c(2.0);
  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  b.Value()  = 3.0;
  b.Deriv(0) = 3.0;
  b.Deriv(1) = 4.0;
  b.Deriv(2) = 5.0;

  //-->TEST 1: += Surreal
  a +=b;
  //test value
  EXPECT_DOUBLE_EQ(5.0, a.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(5.0, a.Deriv(0));
  EXPECT_DOUBLE_EQ(7.0, a.Deriv(1));
  EXPECT_DOUBLE_EQ(9.0, a.Deriv(2));
}

TEST(SurrealCompoundAssignment, AddDouble) {
  //--->This tests the += feature of the Surreal class
  Surreal<double, 3> a;
  double x = 1.1;
  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  //-->TEST 2: += double
  a +=x;
  //test value
  EXPECT_DOUBLE_EQ(3.1, a.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(2.0, a.Deriv(0));
  EXPECT_DOUBLE_EQ(3.0, a.Deriv(1));
  EXPECT_DOUBLE_EQ(4.0, a.Deriv(2));
}

//--> TESTING -=
TEST(SurrealCompoundAssignment, SubtractSurreal) {
  //--->This tests the -= feature of the Surreal class
  Surreal<double, 3> a,b;
  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  b.Value()  = 3.0;
  b.Deriv(0) = 3.0;
  b.Deriv(1) = 5.0;
  b.Deriv(2) = 7.0;

  //-->TEST 1: -= Surreal
  a -=b;
  //test value
  EXPECT_DOUBLE_EQ(-1.0, a.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(-1.0, a.Deriv(0));
  EXPECT_DOUBLE_EQ(-2.0, a.Deriv(1));
  EXPECT_DOUBLE_EQ(-3.0, a.Deriv(2));
}

TEST(SurrealCompoundAssignment, SubtractDouble) {
  //--->This tests the -= feature of the Surreal class
  Surreal<double, 3> a;
  double x = 1.1;
  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  //-->TEST 2: -= double
  a -=x;
  //test value
  EXPECT_DOUBLE_EQ(0.9, a.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(2.0, a.Deriv(0));
  EXPECT_DOUBLE_EQ(3.0, a.Deriv(1));
  EXPECT_DOUBLE_EQ(4.0, a.Deriv(2));
}

//--> TESTING *=
TEST(SurrealCompoundAssignment, MultiplySurreal) {
  //--->This tests the *= feature of the Surreal class
  Surreal<double, 3> a,b;
  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  b.Value()  = 3.0;
  b.Deriv(0) = 3.0;
  b.Deriv(1) = 4.0;
  b.Deriv(2) = 5.0;

  //-->TEST 1: *= Surreal
  a *=b;
  //test value
  EXPECT_DOUBLE_EQ(6.0, a.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(12.0, a.Deriv(0));
  EXPECT_DOUBLE_EQ(17.0, a.Deriv(1));
  EXPECT_DOUBLE_EQ(22.0, a.Deriv(2));
}

TEST(SurrealCompoundAssignment, MultiplyDouble) {
  //--->This tests the *= feature of the Surreal class
  Surreal<double, 3> a;
  double x = 1.1;

  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  //-->TEST 2: *= Double
  a *=x;
  //test value
  EXPECT_DOUBLE_EQ(2.2, a.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(2.2, a.Deriv(0));
  EXPECT_DOUBLE_EQ(3.3, a.Deriv(1));
  EXPECT_DOUBLE_EQ(4.4, a.Deriv(2));
}

//--> TESTING /=
TEST(SurrealCompoundAssignment, DivideBySurreal) {
  //--->This tests the /= feature of the Surreal class
  Surreal<double, 3> a,b;
  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  b.Value()  = 3.0;
  b.Deriv(0) = 3.0;
  b.Deriv(1) = 4.0;
  b.Deriv(2) = 5.0;

  //-->TEST 1: /= Surreal
  a /= b;
  //test value
  EXPECT_DOUBLE_EQ(2.0/3.0, a.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(0.0,     a.Deriv(0));
  EXPECT_DOUBLE_EQ(1.0/9.0, a.Deriv(1));
  EXPECT_DOUBLE_EQ(2.0/9.0, a.Deriv(2));
}

TEST(SurrealCompoundAssignment, DivideByDouble) {
  //--->This tests the /= feature of the Surreal class
  Surreal<double, 3> a;
  double x = 1.1;

  a.Value()  = 2.0;
  a.Deriv(0) = 2.0;
  a.Deriv(1) = 3.0;
  a.Deriv(2) = 4.0;

  //-->TEST 2: /= Double
  a /= x;
  //test value
  EXPECT_DOUBLE_EQ(20.0/11.0, a.Value());
  //test derivatives
  EXPECT_DOUBLE_EQ(20.0/11.0, a.Deriv(0));
  EXPECT_DOUBLE_EQ(30.0/11.0, a.Deriv(1));
  EXPECT_DOUBLE_EQ(40.0/11.0, a.Deriv(2));
}

TEST(SurrealCast, AssignToDouble) {
  Surreal<double,5> y(6.3);
  double x = static_cast<double>(y);
  
  EXPECT_NEAR(x, y.Value(), 5.0e-15);


}

TEST(SurrealCast, AssignToSurrealValue) {
  Surreal<double,5> y(6.3);
  Surreal<double,1> x;
  x.Value() = static_cast<double>(y);
  EXPECT_NEAR(x.Value(), y.Value(), 5.0e-15);


}
