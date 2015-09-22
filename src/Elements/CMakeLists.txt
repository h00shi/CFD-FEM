ADD_SUBDIRECTORY(unit_tests)

ADD_LIBRARY(Elements BarElement.cpp TriElement.cpp TetElement.cpp)
TARGET_LINK_LIBRARIES(Elements BasisFunctions DataStructures)
