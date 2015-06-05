##****************************************************************************80
##! 
##! configure hdf5 libraries
##! 
##! 
##****************************************************************************80
FIND_PACKAGE(HDF5)

IF(HDF5_FOUND)
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
ELSE()
  MESSAGE(FATAL_ERROR "CMAKE could not find hdf5 on your system.")
ENDIF()


