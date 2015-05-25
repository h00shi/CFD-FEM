##****************************************************************************80
##! 
##! configure SuperLU library
##! 
##! 
##****************************************************************************80

FIND_PACKAGE(SUPERLU) 

IF(SUPERLU_FOUND)
  INCLUDE_DIRECTORIES(${SUPERLU_INCLUDE_DIR})
ELSE()
  MESSAGE(FATAL_ERROR "CMAKE was not able to find SuperLU in it's search paths.")
ENDIF()

