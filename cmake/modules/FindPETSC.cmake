##****************************************************************************80
##! 
##! try to find package petsc
##! \xdliang
##! exports 
##!   PETSCLIB
##!   PETSC_INCLUDE_DIR
##****************************************************************************80

FIND_PATH(PETSC_INCLUDE_DIR1
  NAMES petsc.h
  PATHS 
    $ENV{PETSC_HOME}/include
  )

FIND_PATH(PETSC_INCLUDE_DIR2
  NAMES petscconf.h
  PATHS 
    $ENV{PETSC_HOME}/include
  )

FIND_LIBRARY(PETSC_LIB 
  NAMES  libpetsc.so libpetsc.a
  PATHS
    $ENV{PETSC_HOME}/lib
  )


MARK_AS_ADVANCED(
  PETSC_INCLUDE_DIR
  PETSCLIB
  )


SET(PETSC_FOUND FALSE)
IF(PETSC_INCLUDE_DIR1 AND PETSC_INCLUDE_DIR2 AND PETSC_LIB)
  SET(PETSC_FOUND TRUE)
  SET(PETSC_INCLUDE_DIR ${PETSC_INCLUDE_DIR1} ${PETSC_INCLUDE_DIR2})
  MESSAGE(STATUS "Found PETSCLIB: ${PETSC_LIB}")  
  MESSAGE(STATUS "Found PETSC_INCLUDE_DIR: ${PETSC_INCLUDE_DIR}")  
ENDIF()

