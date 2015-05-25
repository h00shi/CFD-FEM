#****************************************************************************80
##! 
##! Setup Compiler
##! \nick
##! 
##****************************************************************************80

SET(_KNOWN_COMPILER FALSE)

IF (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  INCLUDE(setup_compiler_flags_intel)
  SET(_KNOWN_COMPILER TRUE)
ENDIF ()

IF (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  INCLUDE(setup_compiler_flags_gnu)
  SET(_KNOWN_COMPILER TRUE)
ENDIF ()

IF (NOT _KNOWN_COMPILER)
  MESSAGE(FATAL_ERROR "\nUnknown compiler ${CMAKE_CXX_COMPILER_ID}!\n\n")
ELSE ()
# add user defined flags
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_DEFS}")
ENDIF ()
