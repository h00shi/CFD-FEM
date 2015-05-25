##****************************************************************************80
##! 
##! add a mpi test to the test suite and also add a target for memory check.
##! this macro allows variable number of input parameters.
##! definite input parameters:
##!   _test     the name of the test to be added to the test suite
##!   _nproc    the number of processes to be used to run this test
##!   _exe      the name of the executable of this test with absolute path
##! indefinite input (parameters for the executable) will be held in ARGN
##! \qiaolx
##****************************************************************************80

MACRO(ADD_MPI_TEST _test _nproc _exe)
  # first, add the test to the test suite and label it "mpi" 
  SET(_exe_pars_list ${ARGN})
  ADD_TEST(NAME ${_test}
           COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${_nproc} 
             ${_exe} ${_exe_pars_list})
  SET_PROPERTY(TEST ${_test}  PROPERTY LABELS "mpi")
  
  # second, add a custom target for memcheck
  STRING(REPLACE " " ";" _options_list ${MEMORYCHECK_COMMAND_OPTIONS})
 
  ADD_CUSTOM_TARGET(${_test}_memcheck
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${_nproc} valgrind 
      ${_options_list} 
      --suppressions=${MEMORYCHECK_SUPPRESSIONS_FILE}
      ${exe} ${_exe_pars_list}
    VERBATIM)
ENDMACRO()

