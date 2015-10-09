##****************************************************************************80
##! 
##! try to find package doxygen using module
##! \nick
##! exports 
##!   DOXYGEN_BIN_DIR
##!   DOXYGEN_EXECUTABLE
##!   DOXYGEN_FOUND
##****************************************************************************80
FIND_PATH(DOXYGEN_BIN_DIR
  NAMES doxygen
  PATHS 
    $ENV{DOXYGEN_HOME}/bin
    ${CMAKE_SOURCE_DIR}/externals/doxygen/bin
    $ENV{PATH}
  )


MARK_AS_ADVANCED(
  DOXYGEN_BIN_DIR
  )


SET(DOXYGEN_FOUND FALSE)
IF(DOXYGEN_BIN_DIR)
  SET(DOXYGEN_FOUND TRUE)
  SET(DOXYGEN_EXECUTABLE ${DOXYGEN_BIN_DIR}/doxygen)
  MESSAGE(STATUS "Found Doxygen Executable: ${DOXYGEN_EXECUTABLE}")
ENDIF()

