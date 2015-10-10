# Add unit testing directories for Data Structures
ADD_SUBDIRECTORY(unit_tests)

ADD_LIBRARY(ParallelComm Communication.cpp )
TARGET_LINK_LIBRARIES(ParallelComm)
