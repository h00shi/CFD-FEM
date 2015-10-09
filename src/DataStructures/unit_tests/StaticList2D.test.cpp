#include "DataStructures/StaticList2D.h"
#include "my_incl.h"
#include "gtest/gtest.h"
#include <algorithm>
TEST(StaticList2D, Initalize) {
  intT ncol[3] = {2,1,2};
  StaticList2D<intT, 3, 5> x(ncol);
  EXPECT_EQ(x.get_lead_size(), 3);
  EXPECT_EQ(x.get_total_size(), 5);
  EXPECT_EQ(x.get_ncol(0), 2);
  EXPECT_EQ(x.get_ncol(1), 1);
  EXPECT_EQ(x.get_ncol(2), 2);
}

TEST(StaticList2D, Initalize_list) {
 
  StaticList2D<intT, 3, 5> x({2,1,2});
  EXPECT_EQ(x.get_lead_size(), 3);
  EXPECT_EQ(x.get_total_size(), 5);
  EXPECT_EQ(x.get_ncol(0), 2);
  EXPECT_EQ(x.get_ncol(1), 1);
  EXPECT_EQ(x.get_ncol(2), 2);
}

TEST(StaticList2D, Initalize_list2) {
 
  StaticList2D<intT, 3, 5> x({2,1,2},{2,10,20,25,30});
  EXPECT_EQ(x.get_lead_size(), 3);
  EXPECT_EQ(x.get_total_size(), 5);
  EXPECT_EQ(x.get_ncol(0), 2);
  EXPECT_EQ(x.get_ncol(1), 1);
  EXPECT_EQ(x.get_ncol(2), 2);
  
  EXPECT_EQ(x(0,0), 2);
  EXPECT_EQ(x(0,1), 10);
  EXPECT_EQ(x(1,0), 20);
  EXPECT_EQ(x(2,0), 25);
  EXPECT_EQ(x(2,1), 30);
}

TEST(StaticList2D, Access) {
  intT ncol[3] = {2,1,2};
  StaticList2D<intT, 3, 5> x(ncol);
  x(0,1) = 3;
  x(1,0) = 2;
  x(2,0) = 5;
  x(2,1) = 8;

  EXPECT_EQ(3,x(0,1));
  EXPECT_EQ(2,x(1,0));
  EXPECT_EQ(5,x(2,0));
  EXPECT_EQ(8,x(2,1));   
}

TEST(StaticList2D, access2) {
  intT ncol[6] = {1,3,2,5,7,4};
  StaticList2D<double,6,22> x(ncol);
  StaticList2D<double,6,22> const & y = x;

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
TEST(StaticList2D, sortRows) {
  intT ncol[4];
  ncol[0] = 1;
  ncol[1] = 3;
  ncol[2] = 2;
  ncol[3] = 5;
  StaticList2D<double, 4, 11> x(ncol);
  
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

#ifdef DEV_DEBUG
TEST(List2DDeathTest, outOfBounds){
  //death test to make sure program fails when I try to access
  //an out of bounds index for a List2D object
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  intT ncol[6]={1,3,2,5,7,4};
  StaticList2D<double,6,22> x(ncol);
 
  //---> () access tests
  //non constant access function test
  EXPECT_EXIT( x(8,1);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");

  EXPECT_EXIT( x(1,4);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");

  EXPECT_EXIT( x(22);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");

   //---> get_index tests
  //non constant get_index test
  EXPECT_EXIT( x.get_index(6,1);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");

  EXPECT_EXIT( x.get_index(0,3);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");

 
  //--->get_ncol tests
  //non constant get_ncol test
  EXPECT_EXIT( x.get_ncol(321);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");

 

  //--->pointer test
  //non constant
  EXPECT_EXIT( x.get_ptr(7,15);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");

  EXPECT_EXIT( x.get_ptr(3,17);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "In List2D.h - ");

 

  //--->index pointer access test
  //non constant
  EXPECT_EXIT( x.get_index_ptr(123);, ::testing::ExitedWithCode(EXIT_FAILURE),
               "ERROR: In List2D.h - ");
 

}

#endif
