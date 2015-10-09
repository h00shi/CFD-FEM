#****************************************************************************80
# 
# exports 
#   ICC_LIB_DIR
#   ICC_LIBS
##****************************************************************************80

FIND_PATH(ICC_LIB_DIR
  NAMES libimf.so
  PATHS
  /opt/intel/lib/intel64
  )

IF(ICC_LIB_DIR)
  SET(ICC_LIBS 
    imf
    intlc
    irng
    svml)
  
  message(STATUS "Found ICC_LIBS: ${ICC_LIBS}")
ELSE()
  message(FATAL_ERROR "Error: can not found ICC_LIB_DIR")
ENDIF()

MARK_AS_ADVANCED(
  ICC_LIB_DIR
  ICC_LIBS
  )

SET(ICCLIBS_FOUND TRUE)
