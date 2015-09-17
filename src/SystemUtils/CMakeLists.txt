##****************************************************************************80
##! \file CMakeLists.txt
##! \brief setup library "SystemUtils.a"
##! \details Add cpp files to SYS_Source then add library
##! \nick
##****************************************************************************80
# Add unit testing directories for Data Structures
ADD_SUBDIRECTORY(unit_tests)


SET(SYS_Source
  SystemModule.cpp
)

ADD_LIBRARY(SystemUtils ${SYS_Source})

