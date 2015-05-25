##****************************************************************************80
##! 
##! try to find package gtest
##! \qiaolx
##! exports 
##!   GTEST
##!   GTEST_MAIN
##!   GTEST_INCLUDE_DIR
##****************************************************************************80

FIND_PATH(GTEST_INCLUDE_DIR
  NAMES gtest.h
  PATHS 
    $ENV{LD_LIBRARY_PATH}/include/gtest
    $ENV{GTEST_ROOT}/include/gtest
  )

FIND_LIBRARY(GTEST_LIB
  NAMES gtest
  PATHS
    $ENV{GTEST_ROOT}/lib
  )

FIND_LIBRARY(GTEST_MAIN_LIB
  NAMES gtest_main
  PATHS
    $ENV{GTEST_ROOT}/lib
  )

MARK_AS_ADVANCED(
  GTEST_INCLUDE_DIR
  GTEST
  GTEST_MAIN
  )


SET(GTEST_FOUND FALSE)
IF(GTEST_INCLUDE_DIR AND GTEST_LIB AND GTEST_MAIN_LIB)
  SET(GTEST_FOUND TRUE)
  MESSAGE(STATUS "Found GTEST Library: ${GTEST_LIB}")
  MESSAGE(STATUS "Found GTEST_MAIN Library: ${GTEST_MAIN_LIB}")
  MESSAGE(STATUS "Found GTEST Include: ${GTEST_INCLUDE_DIR}")
ENDIF()

