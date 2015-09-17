#include "ParallelComm/Communication.h"
#include "gtest/gtest.h"

TEST(Communication, SerialRank){
  EXPECT_EQ(Communication::GetCommRank(), 0);
  EXPECT_EQ(Communication::GetCommSize(), 1);

}
