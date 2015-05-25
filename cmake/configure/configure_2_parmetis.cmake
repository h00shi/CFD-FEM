##****************************************************************************80
##! 
##! configure metis and parmetis libraries
##! 
##! 
##****************************************************************************80
FIND_PACKAGE(PARMETIS)

IF(PARMETIS_FOUND)
  INCLUDE_DIRECTORIES(${PARMETIS_INCLUDE_DIR})
ELSE()
  MESSAGE(FATAL_ERROR "CMAKE could not find parmetis on your system.")
ENDIF()


