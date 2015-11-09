ADD_SUBDIRECTORY(unit_tests)

ADD_LIBRARY(BasisFunctions H1BarBasis.cpp Polynomials.cpp)
TARGET_LINK_LIBRARIES(BasisFunctions SystemUtils)