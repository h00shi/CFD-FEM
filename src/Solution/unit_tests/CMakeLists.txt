SET(Unit_Tests
    CGElementField.test
)

FOREACH(Unit ${Unit_Tests})
  ADD_EXECUTABLE(${Unit}.unit ${Unit}.cpp )
  TARGET_LINK_LIBRARIES(${Unit}.unit Solution IO DataStructures ${GTEST_LIB} ${GTEST_MAIN_LIB})
  INSTALL(TARGETS ${Unit}.unit DESTINATION "{CMAKE_INSTALL_PREFIX}/test")
  ADD_TEST(NAME ${Unit} COMMAND ${Unit}.unit --gtest_death_test_use_fork)
ENDFOREACH()
# Files to copy for unit tests
FILE (COPY "../../Mesh/unit_tests/Square.grid" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
FILE (COPY "../../IO/unit_tests/Box.msh-bin" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
FILE (COPY "../../IO/unit_tests/Box.idmap" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

