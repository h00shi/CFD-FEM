##****************************************************************************80
##! 
##! test compiler flag using cmake macro CHECK_CXX_COMPILER_FLAG
##! 
##! 
##****************************************************************************80

MACRO(ENABLE_FLAG_IF_SUPPORTED _var _flag)
  STRING(STRIP "${_flag}" _flag_stripped)
  IF(NOT "${_flag_stripped}" STREQUAL "")
    STRING(REGEX REPLACE "^-" "" _flag_name "${_flag_stripped}")
    STRING(REPLACE "++" "__" _flag_name "${_flag_name}")
    CHECK_CXX_COMPILER_FLAG(
      "${_flag_stripped}"
      ${_flag_name}
      )
    IF(${_flag_name})
      SET(${_var} "${${_var}} ${_flag_stripped}")
      STRING(STRIP "${${_var}}" ${_var})
    ENDIF()
  ENDIF()
ENDMACRO()
