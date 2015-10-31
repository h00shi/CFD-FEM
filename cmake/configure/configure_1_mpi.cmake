##***************************************************************************80
##!
##! configure mpi library
##! \qiaolx
##!
##***************************************************************************80
# use the default cmake configuration to find the native mpi libraries.
IF(BUILD_WITH_MPI)
        FIND_PACKAGE(MPI)
        IF (MPI_CXX_FOUND)
           INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
           MESSAGE(STATUS "MPI_CXX_LIBRARIES = ${MPI_CXX_LIBRARIES}")
           MESSAGE(STATUS "MPIEXEC  = ${MPIEXEC}")
           MARK_AS_ADVANCED(CLEAR MPIEXEC)
           ADD_DEFINITIONS(-DHAVE_MPI)
        ELSE ()
        MESSAGE(FATAL_ERROR 
        "CMAKE did not find an MPI distribution. Please make sure that an MPI disct#rubtion is compiled and in your PATH and LD_LIBRARY_PATH environment variables.")
        ENDIF ()
ELSE()
        MESSAGE(STATUS "WARNING: BUILD_WITH_MPI set to OFF...no parallel execution is possible.")
ENDIF()