#include "Array4D.h"
#include "gtest/gtest.h"
#include "my_incl.h"
#include <algorithm>
//--------------------------- ARRAY 2D Tests -----------------------------------
//--->This tests the move constructor of the Array2D class
TEST(Array4D, moveconstructor) {
  Array4D<double> x(5,5,5,5);
  for(intT i = 0; i < 5; i++){
    for(intT j = 0; j < 5; j++){
      for(intT k = 0; k < 5; k++){
	for(intT l = 0; l < 5; l++){
	  x(i,j,k,l) = static_cast<double>(i + j + k + l);
	}
      }
    }
  }
  
  double memory_used = x.get_mem();
  intT size1     = x.get_size(0);
  intT size2     = x.get_size(1);
  intT size3     = x.get_size(2);
  intT size4     = x.get_size(3);
  
  double* begin      = x.begin();
  double* end        = x.end();

  //call the move constructor by explicitly making x into a rvalue
  Array4D<double> y(std::move(x));

  //--->see if the y is equal to what x was
  //--> check memory
  EXPECT_DOUBLE_EQ(memory_used, y.get_mem());

  //--> check size
  EXPECT_EQ(size1, y.get_size(0));
  EXPECT_EQ(size2, y.get_size(1));
  EXPECT_EQ(size3, y.get_size(2));
  EXPECT_EQ(size4, y.get_size(3));

  //--> check pointer of data
  EXPECT_EQ(begin, y.begin());
  EXPECT_EQ(  end, y.end()  );

  //--> check data of y
  for(intT i = 0; i < 5; i++){
    for(intT j = 0; j < 5; j++){
      for(intT k = 0; k < 5; k++){
	for(intT l = 0; l < 5; l++){
	  EXPECT_DOUBLE_EQ( static_cast<double>(i+j+k+l), y(i,j,k,l));
	}
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
  EXPECT_EQ(0, x.get_size(3));

  //--> make sure data is null
  EXPECT_EQ(nullptr, x.begin());
  EXPECT_EQ(nullptr, x.end()  );
}

//--->This tests the move assignment operator of the Array2D class
TEST(Array4D, moveassignmentoperator) {
  Array4D<double> x(5,5,5,5);
  for(intT i = 0; i < 5; i++){
    for(intT j = 0; j < 5; j++){
      for(intT k = 0; k < 5; k++){
	for(intT l = 0; l < 5; l++){
	  x(i,j,k,l) = static_cast<double>(i + j + k + l);
	}
      }
    }
  }
  
  
  double memory_used = x.get_mem();
  intT size1     = x.get_size(0);
  intT size2     = x.get_size(1);
  intT size3     = x.get_size(2);
  intT size4     = x.get_size(3);
  double* begin  = x.begin();
  double* end    = x.end();

  //make a y
  Array4D<double> y(2,2,2,2);
  y(0,0,0,0) = 4.0;
  y(0,0,0,1) = 3.0;
  y(0,0,1,0) = 2.0;
  y(0,0,1,1) = 1.0;

  //call the move assignment operator by explicitly making x into a rvalue
  y = (std::move(x));

  //--->see if the y is equal to what x was
  //--> check memory
  EXPECT_DOUBLE_EQ(memory_used, y.get_mem());

  //--> check size
  EXPECT_EQ(size1, y.get_size(0));
  EXPECT_EQ(size2, y.get_size(1));
  EXPECT_EQ(size3, y.get_size(2));
  EXPECT_EQ(size4, y.get_size(3));
  
  //--> check pointer of data
  EXPECT_EQ(begin, y.begin());
  EXPECT_EQ(end  , y.end()  );

  //--> check data of y
  for(intT i = 0; i < 5; i++){
    for(intT j = 0; j < 5; j++){
      for(intT k = 0; k < 5; k++){
	for(intT l = 0; l < 5; l++){
	  EXPECT_DOUBLE_EQ( static_cast<double>(i+j+k+l), y(i,j,k,l));
	}
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
  EXPECT_EQ(0, x.get_size(3));
  
  //--> make sure data is null
  EXPECT_EQ(nullptr, x.begin());
  EXPECT_EQ(nullptr, x.end()  );
}

TEST(Array4D, access) {
  Array4D<double> x(5,5,5,5);

  x(0,0,0,0) = 1.0;
  x(4,3,1,0) = 4;
 
  EXPECT_EQ(x(0,0,0,0),1.0);
  EXPECT_EQ(x(4,3,1,0),4.0);
  
}

//--->This tests the use of std::sort on a Array2D object
TEST(Array4D, sortentire) {
  Array4D<double> x(1,1,3,3);
  x(0,0,0,0) = 4.0;
  x(0,0,0,1) = 3.0;
  x(0,0,0,2) = 2.0;
  x(0,0,1,0) = -3.0;
  x(0,0,1,1) = -2.0;
  x(0,0,1,2) = -1.0;
  x(0,0,2,0) = 1.0;
  x(0,0,2,1) = 2.0;
  x(0,0,2,2) = 3.0;

  std::sort(x.begin(), x.end());
  EXPECT_DOUBLE_EQ(-3.0, x(0,0,0,0));
  EXPECT_DOUBLE_EQ(-2.0, x(0,0,0,1));
  EXPECT_DOUBLE_EQ(-1.0, x(0,0,0,2));
  EXPECT_DOUBLE_EQ( 1.0, x(0,0,1,0));
  EXPECT_DOUBLE_EQ( 2.0, x(0,0,1,1));
  EXPECT_DOUBLE_EQ( 2.0, x(0,0,1,2));
  EXPECT_DOUBLE_EQ( 3.0, x(0,0,2,0));
  EXPECT_DOUBLE_EQ( 3.0, x(0,0,2,1));
  EXPECT_DOUBLE_EQ( 4.0, x(0,0,2,2));
}

//out of bounds checks are only performed in compiled in debug mode
#ifdef DEV_DEBUG
#ifndef HAVE_MPI
TEST(Array4DDeathTest, outOfBounds){
  //death test to make sure program fails when I try to access
  //an out of bounds index for a 2D array
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  
  Array4D<double> x(4,6,3,5);
  Array4D<double> const & y = x;

  //non constant operator() test
  EXPECT_EXIT( x(5,5,2,0);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array4D.h - ");

  EXPECT_EXIT( x(3,7,1,2);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array4D.h - ");
  
  EXPECT_EXIT( x(3,3,5,0);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array4D.h - ");
  
  EXPECT_EXIT( x(2,3,0,10);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array4D.h - ");

  //non constant get_ptr() test
  EXPECT_EXIT( x.get_ptr(5,5,2,0);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array4D.h - ");

  EXPECT_EXIT( x.get_ptr(3,7,1,2);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array4D.h - ");
  
  EXPECT_EXIT( x.get_ptr(3,3,5,0);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array4D.h - ");
  
  EXPECT_EXIT( x.get_ptr(2,3,0,10);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array4D.h - ");

  
  //constant operator() test
  EXPECT_EXIT( y(18,3,0,0);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array4D.h - ");

  EXPECT_EXIT( y(3,19,1,2);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array4D.h - ");
  
  EXPECT_EXIT( y(3,1,205,0);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array4D.h - ");

  EXPECT_EXIT( y(2,3,0,10);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array4D.h - ");
  
  //non constant get_ptr() test
  EXPECT_EXIT( y.get_ptr(16,3,1,0);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array4D.h - ");

  EXPECT_EXIT( y.get_ptr(3,7,0,1);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array4D.h - ");
  
  EXPECT_EXIT( y.get_ptr(3,0,20,2);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array4D.h - ");
  
  EXPECT_EXIT( y.get_ptr(2,3,0,10);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In Array4D.h - ");

}
#endif
#endif
