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
//--->This tests the move assignment operator of the List2D class
TEST(List2D, moveassignmentoperator) {
  intT ncol[4] = {1,3,2,5};
  StaticList2D<double,4,11> x(ncol);
  intT ncol2[2] = {0,0};
  StaticList2D<double,4,11> y(ncol);

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
