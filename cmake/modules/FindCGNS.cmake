##****************************************************************************80
##! 
##! try to find packages parmetis and metis
##! \nick
##! exports 
##!   CGNS_LIB
##!   CGNS_INCLUDE_DIR
##****************************************************************************80

FIND_PATH(CGNS_INCLUDE_DIR
  NAMES cgnslib.h
  PATHS 
    $ENV{CGNS_HOME}/include
    ${CMAKE_SOURCE_DIR}/externals/cgns/include
  )

FIND_LIBRARY(CGNS_LIB 
  NAMES cgns
  PATHS
    $ENV{CGNS_HOME}/lib
    ${CMAKE_SOURCE_DIR}/externals/cngs/lib
  )

MARK_AS_ADVANCED(
  CGNS_INCLUDE_DIR
  CGNS_LIB
  )

SET(CGNS_FOUND FALSE)
IF(CGNS_INCLUDE_DIR AND CGNS_LIB)
  SET(CGNS_FOUND TRUE)
  MESSAGE(STATUS "CMAKE Found CGNS Library: ${CGNS_LIB}")
  MESSAGE(STATUS "CMAKE Found CGNS Include: ${CGNS_INCLUDE_DIR}")
ENDIF()
