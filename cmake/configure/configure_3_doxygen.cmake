###############################################################################
#
# configure doxygen 
#
###############################################################################

FIND_PACKAGE(Doxygen)
IF(DOXYGEN_FOUND)
ELSE()
  MESSAGE(FATAL_ERROR "Failed in configuring doxygen.")
ENDIF()
