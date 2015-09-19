##****************************************************************************80
##! 
##! configure VTK libraries
##! 
##! 
##****************************************************************************80
FIND_PACKAGE(VTK 6.0 REQUIRED NO_MODULE)

IF(VTK_FOUND)
  INCLUDE(${VTK_USE_FILE})
 
  SET(VTK_INCLUDE_DIR ${VTK_INCLUDE_DIRS})
  #SET(VTK_LIB vtkHybrid)
  #SET(VTK_LIB ${VTK_LIB} vtkWidgets)
  set(VTK_LIB ${VTK_LIBRARIES})
  INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIR})
  MESSAGE(STATUS "CMAKE Found VTK Library: ${VTK_LIB}")
  MESSAGE(STATUS "CMAKE Found VTK Include: ${VTK_INCLUDE_DIR}")
ELSE()
  MESSAGE(FATAL_ERROR "CMAKE could not find VTK on your system.")
ENDIF()


