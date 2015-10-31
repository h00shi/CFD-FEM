###############################################################################
#
# configure petsc library
# \qiaolx
###############################################################################
IF(BUILD_WITH_PETSC)
  FIND_PACKAGE(PETSC)
  IF(PETSC_FOUND)
    add_definitions(-DPETSC)
    INCLUDE_DIRECTORIES(${PETSC_INCLUDE_DIR})
  ELSE()
    MESSAGE(FATAL_ERROR "Failed in founding petsc. Please module load petsc OR turn option BUILD_WITH_PETSC off.")
  ENDIF()
ELSE()
  MESSAGE(STATUS "Skipped safely as BUILD_WITH_PETSC is off.")
ENDIF()
