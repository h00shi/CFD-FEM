ADD_SUBDIRECTORY(unit_tests)

ADD_LIBRARY(Solution 
			Field.cpp
			ElementalField.cpp
			CGElementField.cpp)
TARGET_LINK_LIBRARIES(Solution Mesh SystemUtils ParallelComm DataStructures)
