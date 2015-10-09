##***************************************************************************80
##! 
##! setup ctest variables:
##!   MEMORYCHECK_SUPPRESSIONS_FILE
##!   MEMORYCHECK_COMMAND_OPTIONS
##! \nick
##! 
##***************************************************************************80
#
# ctest memcheck variables
#
SET(MEMORYCHECK_SUPPRESSIONS_FILE
  "${CMAKE_SOURCE_DIR}/cmake/scripts/valgrind_user.supp"
  CACHE FILEPATH "File that contains suppressions for the memory checker")

SET(MEMORYCHECK_COMMAND_OPTIONS "--quiet --tool=memcheck --leak-check=yes --leak-resolution=high --track-origins=yes --show-reachable=yes --num-callers=50"
  CACHE STRING "Options to be used when invoking valgrind"
  )
MARK_AS_ADVANCED(MEMORYCHECK_COMMAND_OPTIONS)
