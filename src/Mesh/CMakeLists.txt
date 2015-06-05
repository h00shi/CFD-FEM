# Add unit testing directories for Data Structures
ADD_SUBDIRECTORY(unit_tests)

ADD_LIBRARY(Mesh UnstMesh.cpp UnstMeshWriterCGNS.cpp)
TARGET_LINK_LIBRARIES(Mesh SystemUtils ${CGNS_LIB} ${HDF5_LIB} ${ZLIB_LIB})