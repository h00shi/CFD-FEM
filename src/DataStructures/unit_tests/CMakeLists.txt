SET(Unit_Tests 
  Array1D.test
  Array2D.test
  Array3D.test
  Array4D.test
  SquareMatrix.test
  List2D.test
  Graph.test
  StaticList2D.test
  EdgeSet.test
)

FOREACH(Unit ${Unit_Tests})
  ADD_EXECUTABLE(${Unit}.unit ${Unit}.cpp)
  TARGET_LINK_LIBRARIES(${Unit}.unit DataStructures 
  ${GTEST_LIB} ${GTEST_MAIN_LIB})
  INSTALL(TARGETS ${Unit}.unit DESTINATION "{CMAKE_INSTALL_PREFIX}/test")
  ADD_TEST(NAME ${Unit} COMMAND ${Unit}.unit --gtest_death_test_use_fork)
ENDFOREACH()

