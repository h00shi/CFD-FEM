##****************************************************************************80
##! \file CMakeLists.txt
##! \brief setup library "SystemUtils.a"
##! \details Add cpp files to SYS_Source then add library
##! \nick
##****************************************************************************80
SET(SYS_Source
  SystemModule.cpp
)

ADD_LIBRARY(SystemUtils ${SYS_Source})

