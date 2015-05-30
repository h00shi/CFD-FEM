#include "Array2D.h"
#include "gtest/gtest.h"
#include "my_incl.h"
#include <algorithm>
//--------------------------- ARRAY 2D Tests -----------------------------------
//--->This tests the move constructor of the Array2D class
TEST(Array2D, moveconstructor) {
  Array2D<double> x(5,5);
  x(0,0) = 4.0;
  x(0,1) = 3.0;
  x(0,2) = 2.0;
  x(0,3) = 1.0;
  x(0,4) = -1.0;

  x(1,0) = 1.0;
  x(1,1) = 2.0;
  x(1,2) = 3.0;
  x(1,3) = 4.0;
  x(1,4) = 5.0;

  double memory_used = x.get_mem();
  unsigned size1     = x.get_size(0);
  unsigned size2     = x.get_size(1);
  double* begin      = x.begin();
  double* end        = x.end();

  //call the move constructor by explicitly making x into a rvalue
  Array2D<double> y(std::move(x));

  //--->see if the y is equal to what x was
  //--> check memory
  EXPECT_DOUBLE_EQ(memory_used, y.get_mem());

  //--> check size
  EXPECT_EQ(size1, y.get_size(0));
  EXPECT_EQ(size2, y.get_size(1));

  //--> check pointer of data
  EXPECT_EQ(begin, y.begin());
  EXPECT_EQ(  end, y.end()  );

  //--> check data of y
  EXPECT_DOUBLE_EQ( 4.0, y(0,0));
  EXPECT_DOUBLE_EQ( 3.0, y(0,1));
  EXPECT_DOUBLE_EQ( 2.0, y(0,2));
  EXPECT_DOUBLE_EQ( 1.0, y(0,3));
  EXPECT_DOUBLE_EQ(-1.0, y(0,4));

  EXPECT_DOUBLE_EQ( 1.0, y(1,0));
  EXPECT_DOUBLE_EQ( 2.0, y(1,1));
  EXPECT_DOUBLE_EQ( 3.0, y(1,2));
  EXPECT_DOUBLE_EQ( 4.0, y(1,3));
  EXPECT_DOUBLE_EQ( 5.0, y(1,4));

  //--->see if the x is set to a default state
  //--> check memory
  EXPECT_DOUBLE_EQ(0.0, x.get_mem());

  //--> check size
  EXPECT_EQ(0, x.get_size(0));
  EXPECT_EQ(0, x.get_size(1));

  //--> make sure data is null
  EXPECT_EQ(nullptr, x.begin());
  EXPECT_EQ(nullptr, x.end()  );
}

//--->This tests the move assignment operator of the Array2D class
TEST(Array2D, moveassignmentoperator) {
  Array2D<double> x(5,5);
  x(0,0) = 4.0;
  x(0,1) = 3.0;
  x(0,2) = 2.0;
  x(0,3) = 1.0;
  x(0,4) = -1.0;

  x(1,0) = 1.0;
  x(1,1) = 2.0;
  x(1,2) = 3.0;
  x(1,3) = 4.0;
  x(1,4) = 5.0;

  double memory_used = x.get_mem();
  unsigned size1     = x.get_size(0);
  unsigned size2     = x.get_size(1);
  double* begin      = x.begin();
  double* end        = x.end();

  //make a y
  Array2D<double> y(2,2);
  y(0,0) = 4.0;
  y(0,1) = 3.0;
  y(1,0) = 2.0;
  y(1,1) = 1.0;

  //call the move assignment operator by explicitly making x into a rvalue
  y = (std::move(x));

  //--->see if the y is equal to what x was
  //--> check memory
  EXPECT_DOUBLE_EQ(memory_used, y.get_mem());

  //--> check size
  EXPECT_EQ(size1, y.get_size(0));
  EXPECT_EQ(size2, y.get_size(1));

  //--> check pointer of data
  EXPECT_EQ(begin, y.begin());
  EXPECT_EQ(end  , y.end()  );

  //--> check data of y
  EXPECT_DOUBLE_EQ( 4.0, y(0,0));
  EXPECT_DOUBLE_EQ( 3.0, y(0,1));
  EXPECT_DOUBLE_EQ( 2.0, y(0,2));
  EXPECT_DOUBLE_EQ( 1.0, y(0,3));
  EXPECT_DOUBLE_EQ(-1.0, y(0,4));

  EXPECT_DOUBLE_EQ( 1.0, y(1,0));
  EXPECT_DOUBLE_EQ( 2.0, y(1,1));
  EXPECT_DOUBLE_EQ( 3.0, y(1,2));
  EXPECT_DOUBLE_EQ( 4.0, y(1,3));
  EXPECT_DOUBLE_EQ( 5.0, y(1,4));

  //--->see if the x is set to a default state
  //--> check memory
  EXPECT_DOUBLE_EQ(0.0, x.get_mem());

  //--> check size
  EXPECT_EQ(0, x.get_size(0));
  EXPECT_EQ(0, x.get_size(1));

  //--> make sure data is null
  EXPECT_EQ(nullptr, x.begin());
  EXPECT_EQ(nullptr, x.end()  );
}

//--->This tests the access of a Array2D object
TEST(Array2D, access) {
  Array2D<double> x(5,5);

  x(0,0) = 1.0;
  x(4,3) = 4;

  EXPECT_EQ(1.0, x(0,0));
  EXPECT_EQ(4.0, x(4,3));
}

//--->This tests the use of std::sort on the rows of a Array2D object
TEST(Array2D, sorteachrow) {
  Array2D<double> x(3,3);
  x(0,0) = 4.0;
  x(0,1) = 3.0;
  x(0,2) = 2.0;
  x(1,0) = -3.0;
  x(1,1) = -2.0;
  x(1,2) = -1.0;
  x(2,0) = 1.0;
  x(2,1) = 2.0;
  x(2,2) = 3.0;

  //sort each row independently
  std::sort(x.get_ptr(0,0), x.get_ptr(0,2)+1);
  std::sort(x.get_ptr(1,0), x.get_ptr(1,2)+1);
  std::sort(x.get_ptr(2,0), x.get_ptr(2,2)+1);

  EXPECT_DOUBLE_EQ( 2.0, x(0,0));
  EXPECT_DOUBLE_EQ( 3.0, x(0,1));
  EXPECT_DOUBLE_EQ( 4.0, x(0,2));
  EXPECT_DOUBLE_EQ(-3.0, x(1,0));
  EXPECT_DOUBLE_EQ(-2.0, x(1,1));
  EXPECT_DOUBLE_EQ(-1.0, x(1,2));
  EXPECT_DOUBLE_EQ( 1.0, x(2,0));
  EXPECT_DOUBLE_EQ( 2.0, x(2,1));
  EXPECT_DOUBLE_EQ( 3.0, x(2,2));
}

//--->This tests the use of std::sort on a Array2D object
TEST(Array2D, sortentire) {
  Array2D<double> x(3,3);
  x(0,0) = 4.0;
  x(0,1) = 3.0;
  x(0,2) = 2.0;
  x(1,0) = -3.0;
  x(1,1) = -2.0;
  x(1,2) = -1.0;
  x(2,0) = 1.0;
  x(2,1) = 2.0;
  x(2,2) = 3.0;

  std::sort(x.begin(), x.end());
  EXPECT_DOUBLE_EQ(-3.0, x(0,0));
  EXPECT_DOUBLE_EQ(-2.0, x(0,1));
  EXPECT_DOUBLE_EQ(-1.0, x(0,2));
  EXPECT_DOUBLE_EQ( 1.0, x(1,0));
  EXPECT_DOUBLE_EQ( 2.0, x(1,1));
  EXPECT_DOUBLE_EQ( 2.0, x(1,2));
  EXPECT_DOUBLE_EQ( 3.0, x(2,0));
  EXPECT_DOUBLE_EQ( 3.0, x(2,1));
  EXPECT_DOUBLE_EQ( 4.0, x(2,2));
}

//out of bounds checks are only performed in compiled in debug mode
#ifdef DEV_DEBUG
#ifndef HAVE_MPI
TEST(Array2DDeathTest, outOfBounds){
  //death test to make sure program fails when I try to access
  //an out of bounds index for a 2D array
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  Array2D<double> x(4,6);
  Array2D<double> const & y = x;

  //non constant operator() test
  EXPECT_EXIT( x(5,5);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array2D.h - ");

  EXPECT_EXIT( x(3,7);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array2D.h - ");

  //non constant get_ptr() test
  EXPECT_EXIT( x.get_ptr(16,3);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array2D.h - ");

  EXPECT_EXIT( x.get_ptr(3,17);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array2D.h - ");

  //constant operator() test
  EXPECT_EXIT( y(18,3);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array2D.h - ");

  EXPECT_EXIT( y(3,19);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array2D.h - ");

  //non constant get_ptr() test
  EXPECT_EXIT( y.get_ptr(16,3);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array2D.h - ");

  EXPECT_EXIT( y.get_ptr(3,17);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array2D.h - ");
}
#endif
#endif
