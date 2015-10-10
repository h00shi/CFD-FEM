# Add unit testing directories for Data Structures
ADD_SUBDIRECTORY(unit_tests)

ADD_LIBRARY(Mesh 
  UnstMeshGeom.cpp
  UnstMeshElements.cpp
  UnstMeshBcFaces.cpp
  UnstMeshEdges.cpp 
  UnstMesh.cpp 
  CGMesh.cpp
  ElementTopology.cpp
  ElementTopologyBar.cpp
  ElementTopologyTriangle.cpp
  ElementTopologyQuadrilateral.cpp
  ElementTopologyTetrahedron.cpp
  ElementTopologyPrism.cpp
  ElementTopologyPyramid.cpp
  ElementTopologyHexahedron.cpp
  DummyMesh.cpp)
TARGET_LINK_LIBRARIES(Mesh IO SystemUtils DataStructures)
