##****************************************************************************80
##! 
##! configure cgns libraries
##! 
##! 
##****************************************************************************80
FIND_PACKAGE(CGNS)

IF(CGNS_FOUND)
  INCLUDE_DIRECTORIES(${CGNS_INCLUDE_DIR})
ELSE()
  MESSAGE(FATAL_ERROR "CMAKE could not find cgns on your system.")
ENDIF()


