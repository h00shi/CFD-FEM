##****************************************************************************80
##! 
##! try to find package mkl
##! \qiaolx
##! exports 
##!   MKL_INCLUDE_DIR
##!   MKL_LIB_DIR
##!   MKL_LIBS
##****************************************************************************80
FIND_PATH(MKL_INCLUDE_DIR
  NAMES mkl.h
  PATHS
    /opt/intel/mkl/include
  )

FIND_PATH(MKL_LIB_DIR
  NAMES libmkl_core.a
  PATHS
    /opt/intel/mkl/lib/intel64
  )

IF(MKL_INCLUDE_DIR)
  message(STATUS "Found MKL_INCLUDE_DIR: ${MKL_INCLUDE_DIR}")
ELSE()
  message(FATAL_ERROR "Error: can not find MKL_INCLUDE_DIR")
ENDIF()

IF(MKL_LIB_DIR)
  SET(MKL_LIBS
    mkl_intel_lp64
    mkl_core
    mkl_sequential  
    pthread
    m
    )

  message(STATUS "Found MKL_LIBS: ${MKL_LIBS}")
ELSE()
  message(FATAL_ERROR "Error: can not found MKL_LIB_DIR")
ENDIF()

MARK_AS_ADVANCED(
  MKL_INCLUDE_DIR
  MKL_LIB_DIR
  MKL_LIBS
  )

SET(MKL_FOUND TRUE)
