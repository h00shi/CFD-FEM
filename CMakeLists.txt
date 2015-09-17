CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")

LIST(APPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/cmake/
  ${CMAKE_SOURCE_DIR}/cmake/modules/
  )

# Load CMAKE macro files
INCLUDE(setup_external_macros)
# load my defined macros
FILE(GLOB my_macro_files "cmake/macros/*.cmake")
FOREACH(ifile ${my_macro_files})
  INCLUDE (${ifile})
ENDFOREACH()

# initialize some cached variables.  This need to be done before project is defined
INCLUDE(setup_cached_variables)

PROJECT(CFD-FEM CXX)
SET(PROJECT_VERSION 1)

#---> Get the Compiler 
INCLUDE(setup_compiler_flags)

IF( CMAKE_BUILD_TYPE MATCHES "Debug")
  MESSAGE(STATUS "C++ Compiler Flags = ${CMAKE_CXX_FLAGS_DEBUG}")
ELSEIF( CMAKE_BUILD_TYPE MATCHES "Release")
   MESSAGE(STATUS "C++ Compiler Flags = ${CMAKE_CXX_FLAGS_RELEASE}")
ELSE()
  MESSAGE(STATUS "WARNING: non-standard build type.  We don't know what you're trying to do so the compiler flags are being set with spartan minimal default values for your compiler.  Please use the Debug or Release CMAKE_BUILD_TYPES to get more predictable behavior")
  MESSAGE(STATUS "C++ Compiler Flags = ${CMAKE_CXX_FLAGS}")
ENDIF()

# setup external libraries
FILE(GLOB my_config_files "cmake/configure/configure_*.cmake")
LIST(SORT my_config_files)
FOREACH(ifile ${my_config_files})
  MESSAGE(STATUS "")
  MESSAGE(STATUS "Include ${ifile}")
  INCLUDE(${ifile})
ENDFOREACH()

INCLUDE(setup_ctest_variables)
INCLUDE(CTest)

########################### INSTALLATION ENVIRONMENT ##########################
SET(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/install")

########################### INCLUDE DIRECTORIES ###############################
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/Include)
#HEADER_DIRECTORIES("${CMAKE_SOURCE_DIR}/src" src_inc_list)
#INCLUDE_DIRECTORIES(${src_inc_list})
#INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src)

ADD_SUBDIRECTORY(src)


############################## PRINT ENVIRONMENT
MESSAGE(STATUS "")
MESSAGE(STATUS "Configuring Build System ${PROJECT_NAME}${PROJECT_VERSION}")
MESSAGE(STATUS "Project source directory: ${PROJECT_SOURCE_DIR}")
MESSAGE(STATUS "Build environment:")
MESSAGE(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
MESSAGE(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")

IF( CMAKE_BUILD_TYPE MATCHES "Debug")
  MESSAGE(STATUS "C++ Compiler Flags = ${CMAKE_CXX_FLAGS_DEBUG}")
ELSEIF( CMAKE_BUILD_TYPE MATCHES "Release")
   MESSAGE(STATUS "C++ Compiler Flags = ${CMAKE_CXX_FLAGS_RELEASE}")
ELSE()
  MESSAGE(STATUS "WARNING: non-standard build type.  We don't know what you're trying to do so the compiler flags are being set with spartan minimal default values for your compiler.  Please use the Debug or Release CMAKE_BUILD_TYPES to get more predictable behavior")
  MESSAGE(STATUS "C++ Compiler Flags = ${CMAKE_CXX_FLAGS}")
ENDIF()

MESSAGE(STATUS "Install directory: ${CMAKE_INSTALL_PREFIX}")
GET_PROPERTY(inc_dirs DIRECTORY ${CMAKE_CURRENT_DIR} PROPERTY INCLUDE_DIRECTORIES)
MESSAGE(STATUS "INCLUDE_DIRECTORIES:")
FOREACH(inc_dir ${inc_dirs})
  MESSAGE(STATUS "  ${inc_dir}")
ENDFOREACH()
