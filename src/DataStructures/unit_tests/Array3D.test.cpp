#include "Array3D.h"
#include "gtest/gtest.h"
#include "my_incl.h"
#include <algorithm>
//--------------------------- ARRAY 2D Tests -----------------------------------
//--->This tests the move constructor of the Array2D class
TEST(Array3D, moveconstructor) {
  Array3D<double> x(5,5,5);
  for(intT i = 0; i < 5; i++){
    for(intT j = 0; j < 5; j++){
      for(intT k = 0; k < 5; k++){
	x(i,j,k) = static_cast<double>(i+j+k);
      }
    }
  }

  double memory_used = x.get_mem();
  intT size1     = x.get_size(0);
  intT size2     = x.get_size(1);
  intT size3     = x.get_size(2);
  double* begin      = x.begin();
  double* end        = x.end();

  //call the move constructor by explicitly making x into a rvalue
  Array3D<double> y(std::move(x));

  //--->see if the y is equal to what x was
  //--> check memory
  EXPECT_DOUBLE_EQ(memory_used, y.get_mem());

  //--> check size
  EXPECT_EQ(size1, y.get_size(0));
  EXPECT_EQ(size2, y.get_size(1));
  EXPECT_EQ(size3, y.get_size(2));
    
  //--> check pointer of data
  EXPECT_EQ(begin, y.begin());
  EXPECT_EQ(  end, y.end()  );

  //--> check data of y
  for(intT i = 0; i < 5; i++){
    for(intT j = 0; j < 5; j++){
      for(intT k = 0; k < 5; k++){
	EXPECT_DOUBLE_EQ( static_cast<double>(i+j+k), y(i,j,k));
      }
    }
  }
  
  //--->see if the x is set to a default state
  //--> check memory
  EXPECT_DOUBLE_EQ(0.0, x.get_mem());

  //--> check size
  EXPECT_EQ(0, x.get_size(0));
  EXPECT_EQ(0, x.get_size(1));
  EXPECT_EQ(0, x.get_size(2));

  //--> make sure data is null
  EXPECT_EQ(nullptr, x.begin());
  EXPECT_EQ(nullptr, x.end()  );
}

//--->This tests the move assignment operator of the Array2D class
TEST(Array3D, moveassignmentoperator) {
  Array3D<double> x(5,5,5);
  for(intT i = 0; i < 5; i++){
    for(intT j = 0; j < 5; j++){
      for(intT k = 0; k < 5; k++){
	x(i,j,k) = static_cast<double>(i+j+k);
      }
    }
  }
  
  double memory_used = x.get_mem();
  intT size1     = x.get_size(0);
  intT size2     = x.get_size(1);
  intT size3     = x.get_size(1);
  double* begin      = x.begin();
  double* end        = x.end();

  //make a y
  Array3D<double> y(2,2,2);
  y(0,0,0) = 4.0;
  y(0,0,1) = 3.0;
  y(0,1,0) = 2.0;
  y(0,1,1) = 1.0;

  //call the move assignment operator by explicitly making x into a rvalue
  y = (std::move(x));

  //--->see if the y is equal to what x was
  //--> check memory
  EXPECT_DOUBLE_EQ(memory_used, y.get_mem());

  //--> check size
  EXPECT_EQ(size1, y.get_size(0));
  EXPECT_EQ(size2, y.get_size(1));
  EXPECT_EQ(size3, y.get_size(2));
  //--> check pointer of data
  EXPECT_EQ(begin, y.begin());
  EXPECT_EQ(end  , y.end()  );
 
  //--> check data of y
  for(intT i = 0; i < 5; i++){
    for(intT j = 0; j < 5; j++){
      for(intT k = 0; k < 5; k++){
	EXPECT_DOUBLE_EQ( static_cast<double>(i+j+k), y(i,j,k));
      }
    }
  }
 

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

TEST(Array3D, access) {
  Array3D<double> x(5,5,5);

  x(0,0,0) = 1.0;
  x(4,3,1) = 4;
 
  EXPECT_EQ(x(0,0,0),1.0);
  EXPECT_EQ(x(4,3,1),4.0);
 
}

//--->This tests the use of std::sort on a Array2D object
TEST(Array3D, sortentire) {
  Array3D<double> x(1,3,3);
  x(0,0,0) = 4.0;
  x(0,0,1) = 3.0;
  x(0,0,2) = 2.0;
  x(0,1,0) = -3.0;
  x(0,1,1) = -2.0;
  x(0,1,2) = -1.0;
  x(0,2,0) = 1.0;
  x(0,2,1) = 2.0;
  x(0,2,2) = 3.0;

  std::sort(x.begin(), x.end());
  EXPECT_DOUBLE_EQ(-3.0, x(0,0,0));
  EXPECT_DOUBLE_EQ(-2.0, x(0,0,1));
  EXPECT_DOUBLE_EQ(-1.0, x(0,0,2));
  EXPECT_DOUBLE_EQ( 1.0, x(0,1,0));
  EXPECT_DOUBLE_EQ( 2.0, x(0,1,1));
  EXPECT_DOUBLE_EQ( 2.0, x(0,1,2));
  EXPECT_DOUBLE_EQ( 3.0, x(0,2,0));
  EXPECT_DOUBLE_EQ( 3.0, x(0,2,1));
  EXPECT_DOUBLE_EQ( 4.0, x(0,2,2));
}

//out of bounds checks are only performed in compiled in debug mode
#ifdef DEV_DEBUG
#ifndef HAVE_MPI
TEST(Array3DDeathTest, outOfBounds){
  //death test to make sure program fails when I try to access
  //an out of bounds index for a 2D array
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  Array3D<double> x(4,6,3);
  Array3D<double> const & y = x;

  //non constant operator() test
  EXPECT_EXIT( x(5,5,2);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array3D.h - ");

  EXPECT_EXIT( x(3,7,1);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array3D.h - ");
  
  EXPECT_EXIT( x(3,3,5);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array3D.h - ");
  //non constant get_ptr() test
  EXPECT_EXIT( x.get_ptr(16,3,1);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array3D.h - ");

  EXPECT_EXIT( x.get_ptr(3,7,0);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array3D.h - ");
  
  EXPECT_EXIT( x.get_ptr(3,0,20);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array3D.h - ");

  //constant operator() test
  EXPECT_EXIT( y(18,3,0);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array3D.h - ");

  EXPECT_EXIT( y(3,19,1);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array3D.h - ");
  
  EXPECT_EXIT( y(3,1,205);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array3D.h - ");
  
  //non constant get_ptr() test
  EXPECT_EXIT( y.get_ptr(16,3,1);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array3D.h - ");

  EXPECT_EXIT( y.get_ptr(3,7,0);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array3D.h - ");
  
  EXPECT_EXIT( y.get_ptr(3,0,20);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array3D.h - ");
  
}
#endif
#endif
