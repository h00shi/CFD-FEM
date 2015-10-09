SET(TESTS
  Surreal.test
  SurrealLogic.test
  SurrealAdd.test
  SurrealSubtract.test
  SurrealMultiply.test
  SurrealDivide.test
  SurrealPow.test
  SurrealUnary.test
)

########################## SERIAL UNIT TESTS ###################################
FOREACH(UNIT ${TESTS})
   ADD_EXECUTABLE(${UNIT}.unit 
                  ${UNIT}.cpp)
   TARGET_LINK_LIBRARIES(${UNIT}.unit ${GTEST_MAIN_LIB} ${GTEST_LIB} pthread)
   ADD_TEST(NAME ${UNIT}
            COMMAND ${UNIT}.unit --gtest_death_test_use_fork)
ENDFOREACH(UNIT)

