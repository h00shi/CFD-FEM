# Add unit testing directories for Data Structures
ADD_SUBDIRECTORY(unit_tests)

ADD_LIBRARY(LinearSolver MKLDSS.cpp)
TARGET_LINK_LIBRARIES(LinearSolver 
SparseMatrix 
DataStructures 
SystemUtils)
