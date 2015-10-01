# Add unit testing directories for Data Structures
ADD_SUBDIRECTORY(unit_tests)

ADD_LIBRARY(Mesh 
  UnstMesh.cpp 
  UnstMeshEdges.cpp 
  ElementTopologyBar.cpp
  ElementTopologyTriangle.cpp
  ElementTopologyQuadrilateral.cpp
  ElementTopologyTetrahedron.cpp
  ElementTopologyPrism.cpp
  ElementTopologyPyramid.cpp
  ElementTopologyHexahedron.cpp
  DummyMesh.cpp)
TARGET_LINK_LIBRARIES(Mesh SystemUtils DataStructures)
