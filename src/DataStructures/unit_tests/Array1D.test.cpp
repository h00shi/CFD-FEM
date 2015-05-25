#include "Array1D.h"
#include "gtest/gtest.h"
#include "my_incl.h"
#include <algorithm>

//--------------------------- ARRAY 1D Tests -----------------------------------
//--->This tests the move constructor of the Array1D class
TEST(Array1D, moveconstructor) {
  Array1D<double> x(5);
  x(0) = 4.0;
  x(1) = 3.0;
  x(2) = 2.0;
  x(3) = 1.0;
  x(4) = -1.0;

  double memory_used = x.get_mem();
  unsigned size      = x.get_size(0);
  double* start      = x.begin();

  //call the move constructor by explicitly making x into a rvalue
  Array1D<double> y(std::move(x));

  //--->see if the y is equal to what x was
  //--> check memory
  EXPECT_DOUBLE_EQ(memory_used, y.get_mem());

  //--> check size
  EXPECT_EQ(size, y.get_size(0));

  //--> check pointer of data
  EXPECT_EQ(start, y.begin());

  //--> check data of y
  EXPECT_DOUBLE_EQ( 4.0, y(0));
  EXPECT_DOUBLE_EQ( 3.0, y(1));
  EXPECT_DOUBLE_EQ( 2.0, y(2));
  EXPECT_DOUBLE_EQ( 1.0, y(3));
  EXPECT_DOUBLE_EQ(-1.0, y(4));

  //--->see if the x is set to a default state
  //--> check memory
  EXPECT_DOUBLE_EQ(0.0, x.get_mem());

  //--> check size
  EXPECT_EQ(0, x.get_size(0));

  //--> check pointer of data
  EXPECT_EQ(nullptr, x.begin());
}

//--->This tests the move assignment operator of the Array1D class
TEST(Array1D, moveassignmentoperator) {
  Array1D<double> x(5);
  x(0) = 4.0;
  x(1) = 3.0;
  x(2) = 2.0;
  x(3) = 1.0;
  x(4) = -1.0;

  double memory_used = x.get_mem();
  unsigned size      = x.get_size(0);
  double* start      = x.begin();

  //make a y
  Array1D<double> y(6);
  y(0) = 4.0;
  y(1) = 3.0;
  y(2) = 2.0;
  y(3) = 1.0;
  y(4) = -1.0;
  y(5) = -1.0;

  //call the move assignment operator by explicitly making x into a rvalue
  y = (std::move(x));

  //--->see if the y is equal to what x was
  //--> check memory
  EXPECT_DOUBLE_EQ(memory_used, y.get_mem());

  //--> check size
  EXPECT_EQ(size, y.get_size(0));

  //--> check pointer of data
  EXPECT_EQ(start, y.begin());

  //--> check data of y
  EXPECT_DOUBLE_EQ( 4.0, y(0));
  EXPECT_DOUBLE_EQ( 3.0, y(1));
  EXPECT_DOUBLE_EQ( 2.0, y(2));
  EXPECT_DOUBLE_EQ( 1.0, y(3));
  EXPECT_DOUBLE_EQ(-1.0, y(4));

  //--->see if the x is set to a default state
  //--> check memory
  EXPECT_DOUBLE_EQ(0.0, x.get_mem());

  //--> check size
  EXPECT_EQ(0, x.get_size(0));

  //--> check pointer of data
  EXPECT_EQ(nullptr, x.begin());
}

//--->This tests the use of std::sort on a Array1D object
TEST(Array1D, sort) {
  Array1D<double> x(5);
  x(0) = 4.0;
  x(1) = 3.0;
  x(2) = 2.0;
  x(3) = 1.0;
  x(4) = -1.0;

  std::sort(x.begin(), x.end());

  EXPECT_EQ(x(0), -1.0);
  EXPECT_EQ(x(1),  1.0);
  EXPECT_EQ(x(2),  2.0);
  EXPECT_EQ(x(3),  3.0);
  EXPECT_EQ(x(4),  4.0);
}

//--->This tests the creation and accessing of an Array1D object
TEST(Array1D, access) {
  Array1D<double> x(5);
  x(0) = 4.0;
  x(1) = 3.0;
  x(2) = 2.0;
  x(3) = 1.0;
  x(4) = -1.0;

  EXPECT_EQ(x(0), 4.0);
  EXPECT_EQ(x(1), 3.0);
  EXPECT_EQ(x(2), 2.0);
  EXPECT_EQ(x(3), 1.0);
  EXPECT_EQ(x(4),-1.0);
}

#ifdef DEV_DEBUG
#ifndef HAVE_MPI
//--->This tests the out of bounds death test of an Array1D object
TEST(Array1DDeathTest, outOfBounds){
  //death test to make sure program fails when I try to access
  //an out of bounds index for a 1D Array
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  Array1D<double> x(5);
  Array1D<double> const & y = x;

  // //non constant operator() test
  // EXPECT_EXIT( x(5);, ::testing::ExitedWithCode(EXIT_FAILURE),
  //               "ERROR:  In Array1D.h - Over bounds on 1st index. "
  //               "Accessing Array1D::(5). Size of Array1D is 5");

//   //non constant get_ptr() test
//   EXPECT_EXIT( x.get_ptr(5);, ::testing::ExitedWithCode(EXIT_FAILURE),
//               "ERROR:  In Array1D.h - "
//                "Attempting to access 1st index (thisone): 5"
//                ". Size of 1st dimension is: 5" );

//   //constant operator() test
//   EXPECT_EXIT( y(5);, ::testing::ExitedWithCode(EXIT_FAILURE),
//                "ERROR:  In Array1D.h - Attempting to access ith index: 5. "
//                "Size of ith dimension is: 5");

//   //constant get_ptr() test
//   EXPECT_EXIT( y.get_ptr(5);, ::testing::ExitedWithCode(EXIT_FAILURE),
//                "ERROR:  In Array1D.h - Attempting to access ith index: 5. "
//                "Size of ith dimension is: 5");
}
#endif
#endif
