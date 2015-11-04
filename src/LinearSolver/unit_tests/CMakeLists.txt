SET(Unit_Tests 
  MKLDSS.test
  PetscSolver.test
)

FOREACH(Unit ${Unit_Tests})
  ADD_EXECUTABLE(${Unit}.unit ${Unit}.cpp )
  TARGET_LINK_LIBRARIES(${Unit}.unit 
    SparseMatrix 
    DataStructures 
    LinearSolver
    IO
    ${GTEST_LIB} 
    ${GTEST_MAIN_LIB})
  INSTALL(TARGETS ${Unit}.unit DESTINATION "{CMAKE_INSTALL_PREFIX}/test")
  ADD_TEST(NAME ${Unit} COMMAND ${Unit}.unit --gtest_death_test_use_fork)
ENDFOREACH()
# Files to copy for unit tests
FILE (COPY "../../Mesh/unit_tests/Square.grid" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
FILE (COPY "petscrc" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
FILE (COPY "Square_quad.msh" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
FILE (COPY "Square_quad.idmap" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
FILE (COPY "Answer-Solve2D_2_MKL.bin" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)