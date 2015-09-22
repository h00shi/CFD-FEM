SET(EXECS 
  PolyViz
)

FOREACH(Unit ${EXECS})
  ADD_EXECUTABLE(${Unit} ${Unit}.cpp)
  TARGET_LINK_LIBRARIES(${Unit} SystemUtils BasisFunctions Elements
   DataStructures)
  INSTALL(TARGETS ${Unit} DESTINATION "{CMAKE_INSTALL_PREFIX}/test") 
ENDFOREACH()
