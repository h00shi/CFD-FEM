##****************************************************************************80
##! 
##! Find the SuperLU Package
##! \nick
##! exports 
##!   SUPERLU
##!   SUPERLU_INCLUDE_DIR
##****************************************************************************80
SET(SUPERLU_SEARCH_PATHS ${CMAKE_SOURCE_DIR}/externals/SuperLU/lib
    $ENV{SUPERLU_HOME}/lib
)

FIND_PATH(SUPERLU_INCLUDE_DIR
  NAMES slu_ddefs.h
  PATHS 
    ${CMAKE_SOURCE_DIR}/externals/SuperLU/include
    $ENV{SUPERLU_HOME}/include
  )

FIND_LIBRARY(SUPERLU_LIB 
  NAMES superlu_4.3
  PATHS
    ${CMAKE_SOURCE_DIR}/externals/SuperLU/lib
    $ENV{SUPERLU_HOME}/lib
  )


MARK_AS_ADVANCED(
  SUPERLU_INCLUDE_DIR
  SUPERLU_LIB
  )

MESSAGE(STATUS "CMAKE is looking for SuperLU in default locations plus: ${SUPERLU_SEARCH_PATH}")

SET(SUPERLU_FOUND FALSE)
IF(SUPERLU_INCLUDE_DIR AND SUPERLU_LIB)
  SET(SUPERLU_FOUND TRUE)
  MESSAGE(STATUS "CMAKE Found SuperLU Library: ${SUPERLU_LIB}")
  MESSAGE(STATUS "CMAKE Found SuperLU Header: ${SUPERLU_INCLUDE_DIR}")
ENDIF()

