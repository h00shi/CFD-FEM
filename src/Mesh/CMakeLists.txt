# Add unit testing directories for Data Structures
ADD_SUBDIRECTORY(unit_tests)

ADD_LIBRARY(Mesh UnstMesh.cpp UnstMeshEdges.cpp ElementTopology.cpp)
TARGET_LINK_LIBRARIES(Mesh SystemUtils DataStructures)
