##***************************************************************************80
##! 
##! setup cached variables prior to PROJECT call 
##! \nick
##! 
##***************************************************************************80
# setup cmake_build_type
SET(CMAKE_BUILD_TYPE
  "Debug"
  CACHE STRING
  "Choose the CMAKE build type, options are: Debug and Release."
  )

IF( (NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Release") AND
    (NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Debug") )
  MESSAGE(FATAL_ERROR
    "INVALID CMAKE_BUILD_TYPE does match either Debug or Release!"
    )
ENDIF()

#
# setup options
#

# enable verbose output from cmake
OPTION(CMAKE_VERBOSE_MAKEFILE "Generate verbose CMake output. Default is on." ON)

OPTION(ICC_WITH_STRICT_ANSI "ICC build with flag -strict_ansi" ON)

# suppress remarks
#   383, reference to temporary value
#   981, operand with unspecified order
#   2547, duplicated include directories
OPTION(ICC_WITH_REMARKS_SUPPRESSION "ICC build with flag -diag-disable 383,981,2547" ON)

OPTION(BUILD_SHARED_LIBS "Create shared libraries. Default is OFF." OFF)

OPTION(BUILD_WITH_LATEX "Build with latex support. Default is ON." ON) 

OPTION(BUILD_WITH_PETSC "Build with petsc library. Default is OFF." OFF)
