##****************************************************************************80
##!   
##! coverage test script to check coverage of all tests
##! \qiaolx
##! 
##****************************************************************************80

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# CTEST_SOURCE_DIRECTORY:
IF ("${CTEST_SOURCE_DIRECTORY}" STREQUAL "")
  # assume that this script resides under cmake/scipts in the source directory
  GET_FILENAME_COMPONENT(_path "${CMAKE_CURRENT_LIST_DIR}" PATH)
  GET_FILENAME_COMPONENT(CTEST_SOURCE_DIRECTORY "${_path}" PATH)

  IF (NOT EXISTS ${CTEST_SOURCE_DIRECTORY}/CMakeLists.txt)
    MESSAGE(FATAL_ERROR "Cann't find CTEST_SOURCE_DIRECTORYdirectory.")
  ENDIF ()
ENDIF ()
MESSAGE ("-- CTEST_SOURCE_DIRECTORY: ${CTEST_SOURCE_DIRECTORY}")

# CTEST_BINARY_DIRECTORY:
IF ("${CTEST_BINARY_DIRECTORY}" STREQUAL "")
  # If CTEST_BINARY_DIRECTORY is not set we just use the current directory
  SET(CTEST_BINARY_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
ENDIF ()
MESSAGE("-- CTEST_BINARY_DIRECTORY: ${CTEST_BINARY_DIRECTORY}")

SET(CTEST_CMAKE_GENERATOR "Unix Makefiles")
MESSAGE("-- CTEST_CMAKE_GENERATOR: ${CTEST_CMAKE_GENERATOR}")


#INCLUDE(${CTEST_SOURCE_DIRECTORY}/CTestConfig.cmake)
SET(CTEST_PROJECT_NAME "PDE++")

SET(WITH_COVERAGE TRUE)

FIND_PROGRAM(CTEST_COVERAGE_COMMAND NAMES gcov)
FIND_PROGRAM(LCOV_COMMAND NAMES lcov)
FIND_PROGRAM(GENHTML_COMMAND NAMES genhtml)
IF (CTEST_COVERAGE_COMMAND)
  MESSAGE("-- coverage_command: ${CTEST_COVERAGE_COMMAND}")
ELSE ()
  MESSAGE(FATAL_ERROR "Failed in configuring gcov")
ENDIF ()

SET(cxx_defs "-DDEV_DEBUG\\ -DMKL_DSS")

SET(CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug")
SET(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} -DCXX_DEFS:String=${cxx_defs} -DWITH_COVERAGE:bool=TRUE ${CTEST_SOURCE_DIRECTORY}")

ctest_start("Experimetal")

ctest_configure()

MESSAGE("-- configuration: \n ${CTEST_CONFIGURE_COMMAND}")

ctest_build()

IF (WITH_BRANCH)
  SET(branch_arg "1")
ELSE()
  SET(branch_arg "0")
ENDIF ()


IF(WITH_COVERAGE)
  # run lcov before running tests to obtain static info as a base
  FILE(STRINGS  ${CTEST_BINARY_DIRECTORY}/CMakeFiles/TargetDirectories.txt myDIRS)
  FOREACH(_dir ${myDIRS})
    MESSAGE("${_dir}")
    FILE(GLOB _var ${_dir}/*.gcno)
    IF("${_var}" STREQUAL "")
      MESSAGE("--skip--")
    ELSE()
      EXECUTE_PROCESS(COMMAND ${LCOV_COMMAND} --capture --initial --directory "${_dir}"  --rc lcov_branch_coverage=${branch_arg}  --output "${_dir}.base.info")
    ENDIF()
  ENDFOREACH()
ENDIF()

ctest_test()

IF (WITH_COVERAGE)
  # run lcov after running tests
  ctest_coverage()
  FILE(STRINGS  ${CTEST_BINARY_DIRECTORY}/CMakeFiles/TargetDirectories.txt myDIRS)
  SET(INIT_INFO TRUE) 
  FOREACH(_dir ${myDIRS})
    MESSAGE("${_dir}")
    FILE(GLOB _var ${_dir}/*.gcda)
    IF("${_var}" STREQUAL "")
      MESSAGE("--skip--")
    ELSE()
      GET_FILENAME_COMPONENT(_test ${_dir} NAME_WE)
      EXECUTE_PROCESS(COMMAND ${LCOV_COMMAND} --capture --directory "${_dir}" --test-name ${_test} --rc lcov_branch_coverage=${branch_arg} --output "${_dir}.coverage.info")
      EXECUTE_PROCESS(COMMAND ${LCOV_COMMAND} -a "${_dir}.base.info" -a "${_dir}.coverage.info" --rc lcov_branch_coverage=${branch_arg} --output "${_dir}.coverage.info")
      IF(INIT_INFO)
        SET(INIT_INFO FALSE)
        EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy "${_dir}.coverage.info" "${CTEST_BINARY_DIRECTORY}/project_coverage.info")
      ELSE()
        EXECUTE_PROCESS(COMMAND ${LCOV_COMMAND} -a "${_dir}.coverage.info" -a "${CTEST_BINARY_DIRECTORY}/project_coverage.info" --rc lcov_branch_coverage=${branch_arg}  --output "${CTEST_BINARY_DIRECTORY}/project_coverage.info")
      ENDIF()
      # clean the binary folder
      EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove "${_dir}.coverage.info" "${_dir}.base.info")
    ENDIF()
  ENDFOREACH()
  
  EXECUTE_PROCESS(COMMAND ${LCOV_COMMAND} --remove "${CTEST_BINARY_DIRECTORY}/project_coverage.info" /opt/CMT/* --rc lcov_branch_coverage=${branch_arg} --output "${CTEST_BINARY_DIRECTORY}/project_coverage.info")
  EXECUTE_PROCESS(COMMAND ${LCOV_COMMAND} --remove "${CTEST_BINARY_DIRECTORY}/project_coverage.info" externals/* --rc lcov_branch_coverage=${branch_arg} --output "${CTEST_BINARY_DIRECTORY}/project_coverage.info")
  EXECUTE_PROCESS(COMMAND ${GENHTML_COMMAND} --show-details --demangle-cpp --legend --rc lcov_branch_coverage=${branch_arg} "${CTEST_BINARY_DIRECTORY}/project_coverage.info" --output-directory "${CTEST_BINARY_DIRECTORY}/coverage")
  # clean the binary folder
  EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove "${CTEST_BINARY_DIRECTORY}/project_coverage.info")
ENDIF ()
