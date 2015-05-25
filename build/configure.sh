#!/bin/bash
rm -rf CMakeCache.txt CMakeFiles cmake_install.cmake CTestTestfile.cmake Makefile Testing 
cmake -D CMAKE_BUILD_TYPE=Debug -D CXX_DEFS:String="-DDEV_DEBUG -DMKL_DSS" -D CMAKE_CXX_COMPILER=mpicxx ../../
