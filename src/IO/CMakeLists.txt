ADD_SUBDIRECTORY(unit_tests)

ADD_LIBRARY(IO UnstMeshWriterVTK.cpp 
               UnstMeshReader.cpp 
               UnstMeshReaderNKBGrid.cpp
               UnstMeshReaderGMSH.cpp)
TARGET_LINK_LIBRARIES(IO Mesh SystemUtils ParallelComm DataStructures
                          ${VTK_LIB} ${ZLIB_LIB})
