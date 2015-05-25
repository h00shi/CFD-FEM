###############################################################################
#
# configure latex library
# \nick
###############################################################################
# use the default cmake configuration to find the native latex libraries.
IF (BUILD_WITH_LATEX)
  FIND_PACKAGE(LATEX)
  IF ("${PDFLATEX_COMPILER}" STREQUAL "")
    MESSAGE(FATAL_ERROR "CMAKE failed to find a latex installation with pdflatex. Please make sure a latex distrubution including pdflatex is installed correctly. Or set OPTION BUILD_WITH_LATEX off.")
  ELSE()
    MESSAGE(STATUS "Found PDFLATEX_COMPILER: ${PDFLATEX_COMPILER}")   
  ENDIF()
  
  IF ("${LATEX2HTML_CONVERTER}" STREQUAL "")
    MESSAGE(FATAL_ERROR "CMAKE failed to find a latex installation withlatex2html. Please make sure a latex distribution including latex2html is installed correctly. Or set OPTION BUILD_WITH_LATEX off.")
  ELSE()
    MESSAGE(STATUS "Found LATEX2HTML_CONVERTER: ${LATEX2HTML_CONVERTER}")
  ENDIF()
ELSE ()
  MESSAGE(STATUS "Will not build latex documentation because the value of BUILD_WITH_LATEX = off.")
  MESSAGE(STATUS "BUILD_WITH_LATEX=${BUILD_WITH_LATEX}")
ENDIF()
