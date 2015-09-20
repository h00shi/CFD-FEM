# Add unit testing directories for Data Structures
ADD_SUBDIRECTORY(unit_tests)

ADD_LIBRARY(Graph Graph.cpp)
TARGET_LINK_LIBRARIES(Graph SystemUtils)
