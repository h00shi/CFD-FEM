#include "ParallelComm/Communication.h"
#include "gtest/gtest.h"

TEST(Communication, SerialRank){
  Communication::Initialize();
  EXPECT_EQ(Communication::GetCommRank(), 0);
  EXPECT_EQ(Communication::GetCommSize(), 1);
  Communication::Finalize();
}
