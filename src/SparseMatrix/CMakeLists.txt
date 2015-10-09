# Add unit testing directories for Data Structures
ADD_SUBDIRECTORY(unit_tests)

ADD_LIBRARY(SparseMatrix SparseMatrix.cpp CSRMatrix.cpp)
TARGET_LINK_LIBRARIES(SparseMatrix DataStructures)
