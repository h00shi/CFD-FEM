#include "gtest/gtest.h"
#include "Surreal/Reverse/SurrealR.h"
#include <stdexcept>

//------------------------ CONSTRUCTORS ----------------------------------------
TEST(SurrealConstructors, DefaultConstructor) {
//--->everything should be initialized to zero!
  SurrealR<double, 3> a;
  //test value
  EXPECT_DOUBLE_EQ(0.0, a.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(0.0, a.GetDeriv(0));
  EXPECT_DOUBLE_EQ(0.0, a.GetDeriv(1));
  EXPECT_DOUBLE_EQ(0.0, a.GetDeriv(2));
}

TEST(SurrealConstructors, ConstructorFromRealNumber) {
//--->set value to input, gradient to zero
  double x = 1.1;
  SurrealR<double, 3> a(x);

  //test value
  EXPECT_DOUBLE_EQ(1.1, a.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(0.0, a.GetDeriv(0));
  EXPECT_DOUBLE_EQ(0.0, a.GetDeriv(1));
  EXPECT_DOUBLE_EQ(0.0, a.GetDeriv(2));
}

TEST(SurrealConstructors, DefaultCopyConstructor) {
//--->copy over values and gradient
  SurrealR<double, 4> a;
  a.SetValue(2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(2, 4.0);

  SurrealR<double, 4> b(a);

  //test value
  EXPECT_DOUBLE_EQ(2.0, b.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(2.0, b.GetDeriv(0));
  EXPECT_DOUBLE_EQ(3.0, b.GetDeriv(1));
  EXPECT_DOUBLE_EQ(4.0, b.GetDeriv(2));
  EXPECT_DOUBLE_EQ(0.0, b.GetDeriv(3));
}

TEST(SurrealConstructors, CopyConstructorFromSurreal) {
//--->copy over values and gradient
  SurrealR<double, 4> a;
  a.SetValue(2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(2, 4.0);

  SurrealRAdd<SurrealR<double, 4>,SurrealR<double, 4>> addition = a+a;
  SurrealR<double, 4> b(addition);

  //test value
  EXPECT_DOUBLE_EQ(4.0, b.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(4.0, b.GetDeriv(0));
  EXPECT_DOUBLE_EQ(6.0, b.GetDeriv(1));
  EXPECT_DOUBLE_EQ(8.0, b.GetDeriv(2));
  EXPECT_DOUBLE_EQ(0.0, b.GetDeriv(3));
}


//--------------------------- ACCESSORS ----------------------------------------
TEST(SurrealAccessors, SetGetValue) {
//--->set and get the value
  SurrealR<double, 3> a;
  a.SetValue(5.0);
  //test value
  EXPECT_DOUBLE_EQ(5.0, a.GetValue());
}

TEST(SurrealAccessors, SetGetDeriv) {
//--->set and get the value
  SurrealR<double, 3> a;
  a.SetDeriv(2, 5.0);
  //test derivative values
  EXPECT_DOUBLE_EQ(5.0, a.GetDeriv(2)); //initialized
  EXPECT_DOUBLE_EQ(0.0, a.GetDeriv(1)); //not initialized
}

TEST(SurrealAccessors, SetGetDeriv_UsingDerivAndValue) {
//--->set and get the value
  SurrealR<double, 3> a;
  SurrealR<double, 3> const & b = a;
  a.Deriv(2) = 5.0;

  EXPECT_DOUBLE_EQ(0.0, a.Value());
  EXPECT_DOUBLE_EQ(0.0, b.Value());

  //test derivative values
  EXPECT_DOUBLE_EQ(5.0, b.Deriv(2)); //initialized
  EXPECT_DOUBLE_EQ(0.0, b.Deriv(1)); //not initialized
}

//------------------------ SURREAL ASSIGNMENT ----------------------------------
//--> TESTING = real
TEST(SurrealAssignment, AssignToReal) {
  SurrealR<double, 3> a;
  double x = 1.1;
  //apply operator=
  a.operator=(x);

  //test value
  EXPECT_DOUBLE_EQ(1.1, a.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(0.0, a.GetDeriv(0));
  EXPECT_DOUBLE_EQ(0.0, a.GetDeriv(1));
  EXPECT_DOUBLE_EQ(0.0, a.GetDeriv(2));
}

//--> TESTING = Surreal
TEST(SurrealAssignment, AssignToSurreal) {
  SurrealR<double, 3> a,b;
  a.SetValue(2.0);
  a.SetDeriv(0,2.0);
  a.SetDeriv(1,3.0);
  a.SetDeriv(2,4.0);

  //assignment to a surreal lvalue
  b = a;

  //test value
  EXPECT_DOUBLE_EQ(2.0, b.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(2.0, b.GetDeriv(0));
  EXPECT_DOUBLE_EQ(3.0, b.GetDeriv(1));
  EXPECT_DOUBLE_EQ(4.0, b.GetDeriv(2));
}

//--> TESTING = Surreal Expression
TEST(SurrealAssignment, AssignToSurrealExpression) {
  SurrealR<double, 3> a,b,c;
  b.SetValue(2.0);
  b.SetDeriv(0,2.0);
  b.SetDeriv(2,4.0);
  b.SetDeriv(1,3.0);

  c.SetValue(3.0);
  c.SetDeriv(1,4.0);
  c.SetDeriv(0,3.0);
  //assignment to a surreal expression
  SurrealRMultiply<SurrealR<double, 3>,SurrealR<double, 3>> expression
    = b*c;
  a = expression;

  //test value
  EXPECT_DOUBLE_EQ(6.0, a.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(12.0, a.GetDeriv(0));
  EXPECT_DOUBLE_EQ(17.0, a.GetDeriv(1));
  EXPECT_DOUBLE_EQ(12.0, a.GetDeriv(2));
}

//------------------------ SURREAL COMPOUND ASSIGNMENT -------------------------
//-----------> operator +=
//--> TESTING += Surreal
TEST(SurrealCompoundAssignment, AddSurreal) {
  //--->This tests the += feature of the Surreal class
  SurrealR<double, 3> a,b;
  a.SetValue(2.0);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);

  b.SetValue(3.0);
  b.SetDeriv(0, 3.0);
  b.SetDeriv(1, 4.0);
  b.SetDeriv(2, 5.0);

  a +=b;
  //test value
  EXPECT_DOUBLE_EQ(5.0, a.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(5.0, a.GetDeriv(0));
  EXPECT_DOUBLE_EQ(7.0, a.GetDeriv(1));
  EXPECT_DOUBLE_EQ(9.0, a.GetDeriv(2));
}

//--> TESTING += Double
TEST(SurrealCompoundAssignment, AddDouble) {
  //--->This tests the += feature of the Surreal class
  SurrealR<double, 3> a;
  double x = 1.1;

  a.SetValue(2.0);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);


  a += x;
  //test value
  EXPECT_DOUBLE_EQ(3.1, a.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(2.0, a.GetDeriv(0));
  EXPECT_DOUBLE_EQ(3.0, a.GetDeriv(1));
  EXPECT_DOUBLE_EQ(4.0, a.GetDeriv(2));
}

//-----------> operator -=
//--> TESTING -= Surreal
TEST(SurrealCompoundAssignment, SubtractSurreal) {
  SurrealR<double, 3> a,b;
  a.SetValue(2.0);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);

  b.SetValue(3.0);
  b.SetDeriv(0, 3.0);
  b.SetDeriv(1, 5.0);
  b.SetDeriv(2, 8.0);

  a -=b;
  //test value
  EXPECT_DOUBLE_EQ(-1.0, a.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(-1.0, a.GetDeriv(0));
  EXPECT_DOUBLE_EQ(-2.0, a.GetDeriv(1));
  EXPECT_DOUBLE_EQ(-4.0, a.GetDeriv(2));
}

//--> TESTING -= Double
TEST(SurrealCompoundAssignment, SubtractDouble) {
  //--->This tests the += feature of the Surreal class
  SurrealR<double, 3> a;
  double x = 1.1;

  a.SetValue(2.0);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);


  a -= x;
  //test value
  EXPECT_DOUBLE_EQ(0.9, a.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(2.0, a.GetDeriv(0));
  EXPECT_DOUBLE_EQ(3.0, a.GetDeriv(1));
  EXPECT_DOUBLE_EQ(4.0, a.GetDeriv(2));
}

//-----------> operator *=
//--> TESTING *= Surreal
TEST(SurrealCompoundAssignment, MultiplySurreal) {
  //--->This tests the *= feature of the Surreal class
  SurrealR<double, 3> a,b;
  a.SetValue(2.0);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);

  b.SetValue(3.0);
  b.SetDeriv(0, 3.0);
  b.SetDeriv(1, 4.0);
  b.SetDeriv(2, 5.0);

  a *=b;
  //test value
  EXPECT_DOUBLE_EQ(6.0, a.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(12.0, a.GetDeriv(0));
  EXPECT_DOUBLE_EQ(17.0, a.GetDeriv(1));
  EXPECT_DOUBLE_EQ(22.0, a.GetDeriv(2));
}

//--> TESTING *= double
TEST(SurrealCompoundAssignment, MultiplyDouble) {
  //--->This tests the *= feature of the Surreal class
  SurrealR<double, 3> a;
  double x = 1.1;

  a.SetValue(2.0);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);

  a *=x;

  //test value
  EXPECT_DOUBLE_EQ(2.2, a.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(2.2, a.GetDeriv(0));
  EXPECT_DOUBLE_EQ(3.3, a.GetDeriv(1));
  EXPECT_DOUBLE_EQ(4.4, a.GetDeriv(2));
}
//-----------> operator /=
//--> TESTING /= Surreal
TEST(SurrealCompoundAssignment, DivideSurreal) {
  //--->This tests the *= feature of the Surreal class
  SurrealR<double, 3> a,b;
  a.SetValue(2.0);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);

  b.SetValue(3.0);
  b.SetDeriv(0, 3.0);
  b.SetDeriv(1, 4.0);
  b.SetDeriv(2, 5.0);

  a /=b;
  //test value
  EXPECT_DOUBLE_EQ(2.0/3.0, a.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(0.0, a.GetDeriv(0));
  EXPECT_DOUBLE_EQ(1.0/9.0, a.GetDeriv(1));
  EXPECT_DOUBLE_EQ(2.0/9.0, a.GetDeriv(2));
}

//--> TESTING /= double
TEST(SurrealCompoundAssignment, DivideDouble) {
  SurrealR<double, 3> a;
  double x = 1.1;
  a.SetValue(2.0);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);

  a /= x;
  //test value
  EXPECT_DOUBLE_EQ(2.0/1.1, a.GetValue());
  //test derivatives
  EXPECT_DOUBLE_EQ(2.0/1.1, a.GetDeriv(0));
  EXPECT_DOUBLE_EQ(3.0/1.1, a.GetDeriv(1));
  EXPECT_DOUBLE_EQ(4.0/1.1, a.GetDeriv(2));
}


#ifdef DEV_DEBUG
//--->This tests the self reference error (the same surreal cannot be on
//    both the left and right hand sides of the equation)
TEST(SurrealRExceptionTest, SelfReference){
  SurrealR<double, 3> a;
  double x = 1.1;

  a.SetValue(2.0);
  a.SetDeriv(0, 2.0);
  a.SetDeriv(1, 3.0);
  a.SetDeriv(2, 4.0);

  EXPECT_THROW( a = a*x;, std::invalid_argument);
}
#endif
