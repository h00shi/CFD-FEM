# Add unit testing directories for Data Structures
ADD_SUBDIRECTORY(unit_tests)

ADD_LIBRARY(DataStructures 
  Array1D.cpp 
  Array2D.cpp 
  Array3D.cpp 
  Array4D.cpp
  List2D.cpp 
  Heap.cpp
  SquareMatrix.cpp 
  Graph.cpp
  EdgeSet.cpp)

TARGET_LINK_LIBRARIES(DataStructures SystemUtils)