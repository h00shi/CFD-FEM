# Add unit testing directories for Data Structures
ADD_SUBDIRECTORY(unit_tests)

ADD_LIBRARY(LinearSolver LinearSolver.cpp MKLDSS.cpp PetscSolver.cpp)
TARGET_LINK_LIBRARIES(LinearSolver 
SparseMatrix 
DataStructures 
SystemUtils
${PETSC_LIB})
