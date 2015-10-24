SET(TESTS
  SurrealR.test
  SurrealRAdd.test
  SurrealRSubtract.test
  SurrealRMultiply.test
  SurrealRDivide.test
  SurrealRPow.test
  SurrealRLogic.test
  SurrealRUnary.test
)

########################## SERIAL UNIT TESTS ###################################
FOREACH(UNIT ${TESTS})
   ADD_EXECUTABLE(${UNIT}.unit 
                  ${UNIT}.cpp)
   TARGET_LINK_LIBRARIES(${UNIT}.unit ${GTEST_MAIN_LIB} ${GTEST_LIB} pthread)
   ADD_TEST(NAME ${UNIT}
            COMMAND ${UNIT}.unit --gtest_death_test_use_fork)
ENDFOREACH(UNIT)

