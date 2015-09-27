##****************************************************************************80
##! 
##! Set Intel compiler flags
##! \nick
##! 
##****************************************************************************80

IF (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "15.0")
  MESSAGE(WARNING "\nYou're using an old version of Intel C++ Compiler!\n")
ENDIF ()

# set common flags
SET(CMAKE_CXX_FLAGS "-std=c++11 -mkl=sequential")

IF (ICC_WITH_STRICT_ANSI)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -strict-ansi")
ENDIF ()

# set flags for debug mode
IF (CMAKE_BUILD_TYPE MATCHES "Debug")
  SET (CMAKE_CXX_FLAGS_DEBUG "-g -O0 -finline -w3 -diag-error 1599,3280")
  IF (ICC_WITH_REMARKS_SUPPRESSION)
    SET(CMAKE_CXX_FLAGS_DEBUG  
      "${CMAKE_CXX_FLAGS_DEBUG} -diag-disable 383,981,2547,424")
  ENDIF ()
# explicitly turn on warnings
# SET(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -finline -Wunused-variable "
#     " -Wunused-function -Wuninitialized -Wreorder -Wreturn-type "
#     " -Wsign-compare -Wshadow -Wpointer-arith -Wnon-virtual-dtor "
#     " -Wformat -Wextra-tokens -Wcheck")
 
#  ENABLE_FLAG_IF_SUPPORTED(CMAKE_CXX_FLAGS_DEBUG "-check=uninit") 
#  ENABLE_FLAG_IF_SUPPORTED(CMAKE_CXX_FLAGS_DEBUG "-early-template-check") 
  SET(CMAKE_CXX_FLAGS_DEUBG "${CMAKE_CXX_FLAGS_DEBUG} -check=uninit")
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -early-template-check")
ENDIF ()


# set flags for release mode
IF (CMAKE_BUILD_TYPE MATCHES "Release")
  SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -finline-limit=300 -no-inline-min-size -no-inline-max-per-compile")
# may generate intel advanced vector extensions 2 instructions
  ENABLE_FLAG_IF_SUPPORTED(CMAKE_CXX_FLAGS_RELEASE "-axCORE-AVX2")
ENDIF ()
