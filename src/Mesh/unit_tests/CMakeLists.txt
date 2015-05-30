SET(Unit_Tests 
 )

FOREACH(Unit ${Unit_Tests})
  ADD_EXECUTABLE(${Unit}.unit ${Unit}.cpp )
  TARGET_LINK_LIBRARIES(${Unit}.unit SystemUtils Mesh ${GTEST_LIB} ${GTEST_MAIN_LIB})
  INSTALL(TARGETS ${Unit}.unit DESTINATION "{CMAKE_INSTALL_PREFIX}/test")
  ADD_TEST(NAME ${Unit} COMMAND ${Unit}.unit --gtest_death_test_use_fork)
ENDFOREACH()