##****************************************************************************80
##!   
##! memory check script to check memory leaks of all tests
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

SET(WITH_MEMCHECK TRUE)

ctest_start("Experimetal")

ctest_configure()

ctest_build()

# get the list of tests and store it in variable tests
set(tmp_file ${CTEST_BINARY_DIRECTORY}/listoftests.txt)
EXECUTE_PROCESS(COMMAND ctest -N
                OUTPUT_FILE ${tmp_file})
FILE(STRINGS ${tmp_file} testsNames REGEX .test)
EXECUTE_PROCESS(COMMAND rm ${tmp_file})
unset(tests)
foreach(_testline ${testsNames})
  string(REGEX REPLACE ".*: ([^ ]+)" "\\1" _test "${_testline}")
  list(APPEND tests ${_test}) # tests is a list
endforeach()

IF (WITH_MEMCHECK)
##  skip parallel tests labelled "mpi"
#  ctest_memcheck(EXCLUDE_LABEL "mpi")
  set(std_err_string "")
  set(failure_counter 0)
  list(LENGTH tests _num_tests)
  set(_test_index 0)
  foreach(_test ${tests})
    MATH(EXPR _test_index "${_test_index} + 1")
    MESSAGE("  Start ${_test_index}/${_num_tests}: ${_test} ")
    SET(_filename
      ${CTEST_BINARY_DIRECTORY}/Testing/Temporary/MemoryChecker.${_test_index}.log)
    if (EXISTS ${_filename})
      EXECUTE_PROCESS(COMMAND rm ${_filename})
    endif ()
    # call ctest to run memcheck on the test that is not labelled "mpi"
    EXECUTE_PROCESS(COMMAND ctest -T memcheck -R ${_test} 
                    -LE "mpi"
                    OUTPUT_QUIET
                    ERROR_VARIABLE _errorstring)
    if (NOT EXISTS ${_filename}) # call make to run memcheck on the test labelled "mpi"
      EXECUTE_PROCESS(COMMAND make ${_test}_memcheck
                      OUTPUT_QUIET
                      ERROR_VARIABLE _errorstring)
      FILE(WRITE ${_filename} "${_errorstring}")
    endif ()

    FILE(STRINGS ${_filename}  _filestring)
    if (NOT "_filestring" STREQUAL "")
      MATH(EXPR failure_counter "${failure_counter} + 1")
      set(_err_message "  MemCheck ${_test_index}: ${_test} ... Failed!")
      set(std_err_string "${std_err_string}  ${_err_message}\n")
      MESSAGE("${_err_message}")
    else ()
      MESSAGE("  MemCheck ${_test_index}: ${_test} ... Passed!")
    endif ()
  endforeach()

  message("  \nSummary of MemCheck:\n")
  if (${failure_counter} GREATER 0)    
    message("${std_err_string}\n  ${failure_counter} tests failed. Please see log files ${CTEST_BINARY_DIRECTORY}/Testing/Temporary/MemoryChecker.TestID.log, where TestId is the ID number of each test.\n" )
  else ()
    message("  All tests passed!\n")
  endif ()
ENDIF ()
