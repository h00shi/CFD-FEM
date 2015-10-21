ADD_SUBDIRECTORY(unit_tests)

ADD_LIBRARY(Solution 
			MeshField.cpp
			NodalField.cpp)
TARGET_LINK_LIBRARIES(Solution Mesh SystemUtils ParallelComm DataStructures)
