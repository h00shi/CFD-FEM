##****************************************************************************80
##! 
##! try to find packages parmetis and metis
##! \nick
##! exports 
##!   METIS
##!   PARMETIS
##!   PARMETIS_INCLUDE_DIR
##****************************************************************************80

FIND_PATH(PARMETIS_INCLUDE_DIR
  NAMES parmetis.h
  PATHS 
    $ENV{PARMETIS_HOME}/include
    ${CMAKE_SOURCE_DIR}/externals/parmetis/include
  )

FIND_PATH(METIS_INCLUDE_DIR
  NAMES metis.h
  PATHS
   $ENV{PARMETIS_HOME}/include
   ${CMAKE_SOURCE_DIR}/externals/parmetis/include
  )

FIND_LIBRARY(PARMETIS_LIB 
  NAMES parmetis
  PATHS
    $ENV{PARMETIS_HOME}/lib
    ${CMAKE_SOURCE_DIR}/externals/parmetis/lib
  )

FIND_LIBRARY(METIS_LIB
  NAMES metis
  PATHS
    $ENV{PARMETIS_HOME}/lib
    ${CMAKE_SOURCE_DIR}/externals/parmetis/lib
  )

MARK_AS_ADVANCED(
  PARMETIS_INCLUDE_DIR
  PARMETIS
  METIS
  )

SET(PARMETIS_FOUND FALSE)
IF(PARMETIS_INCLUDE_DIR AND PARMETIS_LIB AND METIS_LIB)
  SET(PARMETIS_FOUND TRUE)
  MESSAGE(STATUS "CMAKE Found METIS Library: ${METIS_LIB}")
  MESSAGE(STATUS "CMAKE Found METIS Include: ${METIS_INCLUDE_DIR}")
  MESSAGE(STATUS "CMAKE Found PARMETIS Library: ${PARMETIS_LIB}")
  MESSAGE(STATUS "CMAKE Found PARMETIS Include: ${PARMETIS_INCLUDE_DIR}")
ENDIF()
