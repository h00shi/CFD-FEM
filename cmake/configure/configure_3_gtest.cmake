##****************************************************************************80
##! 
##! configure gtest library
##! 
##! 
##****************************************************************************80

FIND_PACKAGE(GTEST) 

IF(GTEST_FOUND)
  INCLUDE_DIRECTORIES(${GTEST_INCLUDE_DIR})
ELSE()
  MESSAGE(FATAL_ERROR "CMAKE could not find Google Test (GTEST) installed on your system.")
ENDIF()

