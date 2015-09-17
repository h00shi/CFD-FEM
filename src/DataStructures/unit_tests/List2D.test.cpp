#include "DataStructures/List2D.h"
#include "my_incl.h"
#include "gtest/gtest.h"
#include <algorithm>
//--------------------------- LIST 2D Tests -----------------------------------
//--->This tests the move constructor of the List2D class
TEST(Array2D, moveconstructor) {
  List2D<double> x(4,11);
  Array1D<int> ncol(x.get_lead_size());

  // 1
  // 3 2 1
  // 6 1
  // 1 2 3 5 -1
  ncol(0) = 1;
  ncol(1) = 3;
  ncol(2) = 2;
  ncol(3) = 5;
  x.set_ncol(ncol);

  x.set_value(0.0);
  x(0,0) = 1.0;

  x(1,0) = 3.0;
  x(1,1) = 2.0;
  x(1,2) = 1.0;

  x(2,0) = 6.0;
  x(2,1) = 1.0;

  x(3,0) = 1.0;
  x(3,1) = 2.0;
  x(3,2) = 3.0;
  x(3,3) = 5.0;
  x(3,4) =-1.0;

  double memory_used = x.get_mem();
  unsigned size1     = x.get_lead_size();
  unsigned size      = x.get_total_size();
  double* begin      = x.begin();
  double* end        = x.end();

  //call the move constructor by explicitly making x into a rvalue
  List2D<double> y(std::move(x));

  //--->see if the y is equal to what x was
  //--> check memory
  EXPECT_DOUBLE_EQ(memory_used, y.get_mem());

  //--> check size
  EXPECT_EQ(size1, y.get_lead_size());
  EXPECT_EQ(size,  y.get_total_size());

  //--> check pointer of data
  EXPECT_EQ(begin, y.begin());
  EXPECT_EQ(end,   y.end()  );

  //--> check data of y
  EXPECT_DOUBLE_EQ( 1.0, y(0,0));
  EXPECT_DOUBLE_EQ( 3.0, y(1,0));
  EXPECT_DOUBLE_EQ( 2.0, y(1,1));
  EXPECT_DOUBLE_EQ( 1.0, y(1,2));
  EXPECT_DOUBLE_EQ( 6.0, y(2,0));
  EXPECT_DOUBLE_EQ( 1.0, y(2,1));
  EXPECT_DOUBLE_EQ( 1.0, y(3,0));
  EXPECT_DOUBLE_EQ( 2.0, y(3,1));
  EXPECT_DOUBLE_EQ( 3.0, y(3,2));
  EXPECT_DOUBLE_EQ( 5.0, y(3,3));
  EXPECT_DOUBLE_EQ(-1.0, y(3,4));

  //--->see if the x is set to a default state
  //--> check memory
  EXPECT_DOUBLE_EQ(0.0, x.get_mem());

  //--> check size
  EXPECT_EQ(0, x.get_lead_size());
  EXPECT_EQ(0, x.get_total_size());

  //--> make sure data is null
  EXPECT_EQ(nullptr, x.begin());
  EXPECT_EQ(nullptr, x.end()  );
}

//--->This tests the move assignment operator of the List2D class
TEST(Array2D, moveassignmentoperator) {
  List2D<double> x(4,11);
  List2D<double> y(2,3);

  Array1D<int> ncol(x.get_lead_size());

  // 1
  // 3 2 1
  // 6 1
  // 1 2 3 5 -1
  ncol(0) = 1;
  ncol(1) = 3;
  ncol(2) = 2;
  ncol(3) = 5;
  x.set_ncol(ncol);

  x.set_value(0.0);
  x(0,0) = 1.0;

  x(1,0) = 3.0;
  x(1,1) = 2.0;
  x(1,2) = 1.0;

  x(2,0) = 6.0;
  x(2,1) = 1.0;

  x(3,0) = 1.0;
  x(3,1) = 2.0;
  x(3,2) = 3.0;
  x(3,3) = 5.0;
  x(3,4) =-1.0;

  double memory_used = x.get_mem();
  unsigned size1     = x.get_lead_size();
  unsigned size      = x.get_total_size();
  double* begin      = x.begin();
  double* end        = x.end();

  //call the move assignment operator by explicitly making x into a rvalue
  y = (std::move(x));

  //--->see if the y is equal to what x was
  //--> check memory
  EXPECT_DOUBLE_EQ(memory_used, y.get_mem());

  //--> check size
  EXPECT_EQ(size1, y.get_lead_size());
  EXPECT_EQ(size,  y.get_total_size());

  //--> check pointer of data
  EXPECT_EQ(begin, y.begin());
  EXPECT_EQ(  end, y.end()  );

  //--> check data of y
  EXPECT_DOUBLE_EQ( 1.0, y(0,0));
  EXPECT_DOUBLE_EQ( 3.0, y(1,0));
  EXPECT_DOUBLE_EQ( 2.0, y(1,1));
  EXPECT_DOUBLE_EQ( 1.0, y(1,2));
  EXPECT_DOUBLE_EQ( 6.0, y(2,0));
  EXPECT_DOUBLE_EQ( 1.0, y(2,1));
  EXPECT_DOUBLE_EQ( 1.0, y(3,0));
  EXPECT_DOUBLE_EQ( 2.0, y(3,1));
  EXPECT_DOUBLE_EQ( 3.0, y(3,2));
  EXPECT_DOUBLE_EQ( 5.0, y(3,3));
  EXPECT_DOUBLE_EQ(-1.0, y(3,4));

  //--->see if the x is set to a default state
  //--> check memory
  EXPECT_DOUBLE_EQ(0.0, x.get_mem());

  //--> check size
  EXPECT_EQ(0, x.get_lead_size());
  EXPECT_EQ(0, x.get_total_size());

  //--> make sure data is null
  EXPECT_EQ(nullptr, x.begin());
  EXPECT_EQ(nullptr, x.end()  );
}

//--->This tests all the possible data access of a List2D object with
//    size construction
TEST(List2D, access) {
  List2D<double> x(6,22);
  Array1D<int> ncol(x.get_lead_size());
  List2D<double> const & y = x;

  // x
  // x x x
  // x x
  // x x x x x
  // x x x x x x x
  // x x x x
  ncol(0) = 1;
  ncol(1) = 3;
  ncol(2) = 2;
  ncol(3) = 5;
  ncol(4) = 7;
  ncol(5) = 4;
  x.set_ncol(ncol);

  x.set_value(0.0);
  x(0,0) = 1.0;
  x(4,3) = 4.0;

  //(i,j) indexing
  EXPECT_DOUBLE_EQ(1.0, x(0,0));
  EXPECT_DOUBLE_EQ(4.0, x(4,3));

  //(i) indexing
  EXPECT_DOUBLE_EQ(1.0, x(0));
  EXPECT_DOUBLE_EQ(4.0, x(14));

  //(i,j) indexing - constant
  EXPECT_DOUBLE_EQ(1.0, y(0,0));
  EXPECT_DOUBLE_EQ(4.0, y(4,3));

  //(i) indexing - constant
  EXPECT_DOUBLE_EQ(1.0, y(0));
  EXPECT_DOUBLE_EQ(4.0, y(14));
}

//--->This tests all the possible data access of a List2D object with
//    default construction
TEST(List2D, access2) {
  List2D<double> x;
  x.initialize(6,22);
  Array1D<int> ncol(x.get_lead_size());
  List2D<double> const & y = x;

  ncol(0) = 1;
  ncol(1) = 3;
  ncol(2) = 2;
  ncol(3) = 5;
  ncol(4) = 7;
  ncol(5) = 4;
  x.set_ncol(ncol);

  x.set_value(0.0);
  x(0,0) = 1.0;
  x(4,3) = 4.0;

  //(i,j) indexing
  EXPECT_DOUBLE_EQ(1.0, x(0,0));
  EXPECT_DOUBLE_EQ(4.0, x(4,3));

  //(i) indexing
  EXPECT_DOUBLE_EQ(1.0, x(0));
  EXPECT_DOUBLE_EQ(4.0, x(14));

  //(i,j) indexing - constant
  EXPECT_DOUBLE_EQ(1.0, y(0,0));
  EXPECT_DOUBLE_EQ(4.0, y(4,3));

  //(i) indexing - constant
  EXPECT_DOUBLE_EQ(1.0, y(0));
  EXPECT_DOUBLE_EQ(4.0, y(14));
}

//--->This tests thse use of std::Sort on the rows of List2D
TEST(List2D, sortRows) {
  List2D<double> x(4,11);
  Array1D<int> ncol(x.get_lead_size());

  // 1
  // 3 2 1
  // 6 1
  // 1 2 3 5 -1
  ncol(0) = 1;
  ncol(1) = 3;
  ncol(2) = 2;
  ncol(3) = 5;
  x.set_ncol(ncol);

  x.set_value(0.0);
  x(0,0) = 1.0;

  x(1,0) = 3.0;
  x(1,1) = 2.0;
  x(1,2) = 1.0;

  x(2,0) = 6.0;
  x(2,1) = 1.0;

  x(3,0) = 1.0;
  x(3,1) = 2.0;
  x(3,2) = 3.0;
  x(3,3) = 5.0;
  x(3,4) =-1.0;

  //sort each row independently
  std::sort(x.get_ptr(0,0), x.get_ptr(0,0)+1);
  std::sort(x.get_ptr(1,0), x.get_ptr(1,2)+1);
  std::sort(x.get_ptr(2,0), x.get_ptr(2,1)+1);
  std::sort(x.get_ptr(3,0), x.get_ptr(3,4)+1);

  //--->post row sort
  //row0
  EXPECT_DOUBLE_EQ(1.0, x(0,0));

  //row1
  EXPECT_DOUBLE_EQ(1.0, x(1,0));
  EXPECT_DOUBLE_EQ(2.0, x(1,1));
  EXPECT_DOUBLE_EQ(3.0, x(1,2));

  //row2
  EXPECT_DOUBLE_EQ(1.0, x(2,0));
  EXPECT_DOUBLE_EQ(6.0, x(2,1));

  //row3
  EXPECT_DOUBLE_EQ(-1.0, x(3,0));
  EXPECT_DOUBLE_EQ( 1.0, x(3,1));
  EXPECT_DOUBLE_EQ( 2.0, x(3,2));
  EXPECT_DOUBLE_EQ( 3.0, x(3,3));
  EXPECT_DOUBLE_EQ( 5.0, x(3,4));
}

//--->This tests using begin and end to sort an entire list2d
TEST(List2D, sort_entire_list) {
  List2D<double> x(4,11);
  Array1D<int> ncol(x.get_lead_size());

  // 1
  // 3 2 1
  // 6 1
  // 1 2 3 5 -1
  ncol(0) = 1;
  ncol(1) = 3;
  ncol(2) = 2;
  ncol(3) = 5;
  x.set_ncol(ncol);

  x.set_value(0.0);
  x(0,0) = 1.0;

  x(1,0) = 3.0;
  x(1,1) = 2.0;
  x(1,2) = 1.0;

  x(2,0) = 6.0;
  x(2,1) = 1.0;

  x(3,0) = 1.0;
  x(3,1) = 2.0;
  x(3,2) = 3.0;
  x(3,3) = 5.0;
  x(3,4) =-1.0;

  //sort entire list
  std::sort(x.begin(), x.end());

  //--->post row sort
  EXPECT_DOUBLE_EQ(-1.0, x(0,0));
  EXPECT_DOUBLE_EQ( 1.0, x(1,0));
  EXPECT_DOUBLE_EQ( 1.0, x(1,1));
  EXPECT_DOUBLE_EQ( 1.0, x(1,2));
  EXPECT_DOUBLE_EQ( 1.0, x(2,0));
  EXPECT_DOUBLE_EQ( 2.0, x(2,1));
  EXPECT_DOUBLE_EQ( 2.0, x(3,0));
  EXPECT_DOUBLE_EQ( 3.0, x(3,1));
  EXPECT_DOUBLE_EQ( 3.0, x(3,2));
  EXPECT_DOUBLE_EQ( 5.0, x(3,3));
  EXPECT_DOUBLE_EQ( 6.0, x(3,4));
}

//--->This tests pointer access
TEST(List2D, pointer_return){
  List2D<double> x(6,22);
  Array1D<int> ncol(x.get_lead_size());

  ncol(0) = 1;
  ncol(1) = 3;
  ncol(2) = 2;
  ncol(3) = 5;
  ncol(4) = 7;
  ncol(5) = 4;
  x.set_ncol(ncol);

  x.set_value(0.0);
  x(0,0) = 1.0;
  x(4,3) = 4.0;
  const double* p = x.get_ptr(3,2);

  EXPECT_DOUBLE_EQ(p[0], x(3,2));
  EXPECT_DOUBLE_EQ(p[1], x(3,3));
  EXPECT_DOUBLE_EQ(p[2], x(3,4));
}

//--->This tests pointer access
TEST(List2D, pointer_return2){
  List2D<double> x;
  x.initialize(6,22);
  Array1D<int> ncol(x.get_lead_size());

  ncol(0) = 1;
  ncol(1) = 3;
  ncol(2) = 2;
  ncol(3) = 5;
  ncol(4) = 7;
  ncol(5) = 4;
  x.set_ncol(ncol);

  x.set_value(0.0);
  x(0,0) = 1.0;
  x(4,3) = 4.0;
  const double* p = x.get_ptr(3,2);

  EXPECT_DOUBLE_EQ(x(3,2), p[0]);
  EXPECT_DOUBLE_EQ(x(3,3), p[1]);
  EXPECT_DOUBLE_EQ(x(3,4), p[2]);
}

#ifdef DEV_DEBUG
TEST(List2DDeathTest, outOfBounds){
  //death test to make sure program fails when I try to access
  //an out of bounds index for a List2D object
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  List2D<double> x(6,22);
  List2D<double> const & y = x;
  Array1D<int> ncol(x.get_lead_size());

  // x
  // x x x
  // x x
  // x x x x x
  // x x x x x x x
  // x x x x
  ncol(0) = 1;
  ncol(1) = 3;
  ncol(2) = 2;
  ncol(3) = 5;
  ncol(4) = 7;
  ncol(5) = 4;
  x.set_ncol(ncol);

  //---> () access tests
  //non constant access function test
  EXPECT_EXIT( x(8,1);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");

  EXPECT_EXIT( x(1,4);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");

  EXPECT_EXIT( x(22);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");

  //constant access function test
  EXPECT_EXIT( y(15,1);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");
  EXPECT_EXIT( y(4,14);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");
  EXPECT_EXIT( y(23);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");

  //---> get_index tests
  //non constant get_index test
  EXPECT_EXIT( x.get_index(6,1);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");

  EXPECT_EXIT( x.get_index(0,3);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");

  //constant get_index test
  EXPECT_EXIT( y.get_index(55,231);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");

  EXPECT_EXIT( y.get_index(2,3);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");

  //--->get_ncol tests
  //non constant get_ncol test
  EXPECT_EXIT( x.get_ncol(321);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");

  //constant get_ncol test
  EXPECT_EXIT( y.get_ncol(321);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");

  //--->pointer test
  //non constant
  EXPECT_EXIT( x.get_ptr(7,15);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");

  EXPECT_EXIT( x.get_ptr(3,17);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "In List2D.h - ");

  //constant
  EXPECT_EXIT( y.get_ptr(7,15);, ::testing::ExitedWithCode(EXIT_FAILURE),
              "ERROR: In List2D.h - " );

  EXPECT_EXIT( y.get_ptr(3,17);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");

  //--->index pointer access test
  //non constant
  EXPECT_EXIT( x.get_index_ptr(123);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");
  //constant
  EXPECT_EXIT( y.get_index_ptr(123);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");

}

TEST(List2DDeathTest, listTooLarge){
  //death test to make sure program fails when I try to prescribe column
  // sizes that make the entire structure larger than the preallocated size
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  List2D<double> x(6,22);
  Array1D<int> ncol(x.get_lead_size());

  //initialize x such that it is at maximum size
  // x
  // x x x
  // x x
  // x x x x x
  // x x x x x x x
  // x x x x
  ncol(0) = 1;
  ncol(1) = 3;
  ncol(2) = 2;
  ncol(3) = 5;
  ncol(4) = 7;
  ncol(5) = 4;
  x.set_ncol(ncol); //list2d is at maximum size already

  ncol(5) = ncol(5) + 1; //increase request ncol of row 5 by one

  //increase ncol tests
  EXPECT_EXIT( x.set_ncol(ncol);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "In List2D.h - ERROR: Requesting too many elements.");

  ncol(5) = 4;   //revert back to old request
  x.set_ncol(ncol); //this should work now.

  ncol(0) = 100; //request something too large again
  EXPECT_EXIT( x.set_ncol(ncol);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "In List2D.h - ERROR: Requesting too many elements.");

  ncol(0) = 1; //revert back
  x.set_ncol(ncol); //this should work now.
  ncol(1) = 4; //request something too large again
  EXPECT_EXIT( x.set_ncol(ncol);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "In List2D.h - ERROR: Requesting too many elements.");
}

TEST(List2DDeathTest, ncolMismatch){
  //death test to make sure program fails when I try to prescribe
  // too few column sizes
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  List2D<double> x(6,50);
  Array1D<int> ncolTooFew(x.get_lead_size()-1);


  ncolTooFew(0) = 1;
  ncolTooFew(1) = 3;
  ncolTooFew(2) = 2;
  ncolTooFew(3) = 5;
  ncolTooFew(4) = 7;

  //too few ncol specified
  EXPECT_EXIT( x.set_ncol(ncolTooFew);,
               ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: Mismatch of num rows and num of col specification.");

}
#endif
