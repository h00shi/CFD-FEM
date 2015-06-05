##****************************************************************************80
##! 
##! try to find hdf5 from PATH set by module load command
##! \nick
##! exports 
##!   HDF5_LIB
##!   HDF5_INCLUDE_DIR
##****************************************************************************80
FIND_PATH(HDF5_INCLUDE_DIR
  NAMES hdf5.h
  PATHS 
    $ENV{HDF5_HOME}/include
    ${CMAKE_SOURCE_DIR}/externals/hdf5/include
  )

FIND_LIBRARY(HDF5_LIB 
  NAMES hdf5
  PATHS
    $ENV{HDF5_HOME}/lib
    ${CMAKE_SOURCE_DIR}/externals/hdf5/lib
  )

MARK_AS_ADVANCED(
  HDF5_INCLUDE_DIR
  HDF5_LIB
  )

SET(HDF5_FOUND FALSE)
IF(HDF5_INCLUDE_DIR AND HDF5_LIB)
  SET(HDF5_FOUND TRUE)
  MESSAGE(STATUS "CMAKE Found HDF5 Library: ${HDF5_LIB}")
  MESSAGE(STATUS "CMAKE Found HDF5 Include: ${HDF5_INCLUDE_DIR}")
ENDIF()

