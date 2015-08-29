SET(Unit_Tests 
  Array1D.test
  Array2D.test
  Array3D.test
  Array4D.test
  SquareMatrix.test
  List2D.test
  CSRMatrix.test
)

FOREACH(Unit ${Unit_Tests})
  ADD_EXECUTABLE(${Unit}.unit ${Unit}.cpp)
  TARGET_LINK_LIBRARIES(${Unit}.unit SystemUtils Mesh ${GTEST_LIB} ${GTEST_MAIN_LIB})
  INSTALL(TARGETS ${Unit}.unit DESTINATION "{CMAKE_INSTALL_PREFIX}/test")
  ADD_TEST(NAME ${Unit} COMMAND ${Unit}.unit --gtest_death_test_use_fork)
ENDFOREACH()

# Copy Square.grid
file(COPY Square.grid DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
