##****************************************************************************80
##! 
##! set GNU compiler flags
##! \nick
##! 
##****************************************************************************80

IF (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8" )
  MESSAGE(WARNING "\nYou're using an old version of GNU C++ Compiler!\n")
ENDIF ()

# set common flags
SET (CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wunused-variable -Wextra")

# colorizing diagnostics emitted by gcc when outputting to terminals
ENABLE_FLAG_IF_SUPPORTED(CMAKE_CXX_FLAGS "-fdiagnostics-color=auto")

MESSAGE(STATUS "Currently some Intel compiler libraries are required, "
 "as third party libaries are compiled by Intel compiler.")

FIND_PACKAGE(ICCLIBS)
IF (ICCLIBS_FOUND)
  LINK_DIRECTORIES(${ICC_LIB_DIR})
  SET(SABRES_LIBS_OTHER 
    ${SABRES_LIBS_OTHER}
    ${ICC_LIBS})
ELSE ()
  MESSAGE(FATAL_ERROR "Cannot find necessary Intel compiler libraries."
    "Please load Intel compiler module.")
ENDIF ()

FIND_PACKAGE(MKL)
IF(MKL_FOUND)
  INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIR})
  LINK_DIRECTORIES(${MKL_LIB_DIR})
  SET(SABRES_LIBS_OTHER 
    ${SABRES_LIBS_OTHER}
    ${MKL_LIBS})
ENDIF()

# set flags for debug mode
IF (CMAKE_BUILD_TYPE MATCHES "Debug")
  SET(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
  IF (WITH_COVERAGE) 
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-inline --coverage")
  ELSE ()
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -finline")
  ENDIF ()
ENDIF ()


# set flags for release mode

IF (CMAKE_BUILD_TYPE MATCHES "Release")
  SET(CMAKE_CXX_FLAGS_RELEASE "-O3")

# may produce code optimized for the local machine
  ENABLE_FLAG_IF_SUPPORTED(CMAKE_CXX_FLAGS_RELEASE "-march=native")
  
ENDIF ()



