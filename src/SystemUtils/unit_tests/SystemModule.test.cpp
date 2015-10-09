#include "SystemUtils/SystemModule.h"
#include "gtest/gtest.h"

TEST(SystemModule, pout)
{
  SystemModule::cout << "Test" << std::endl;
  EXPECT_EQ(SystemModule::cout.GetOutputRank(), 0);

  EXPECT_EQ(SystemModule::pout(2).GetOutputRank(),2); 
}
