ADD_SUBDIRECTORY(unit_tests)

ADD_LIBRARY(IO UnstMeshWriterCGNS.cpp 
               UnstMeshWriterVTK.cpp 
               UnstMeshReader.cpp 
               UnstMeshReaderNKBGrid.cpp)
TARGET_LINK_LIBRARIES(IO Mesh SystemUtils ParallelComm DataStructures
                          ${VTK_LIB} ${CGNS_LIB} ${HDF5_LIB} ${ZLIB_LIB})
