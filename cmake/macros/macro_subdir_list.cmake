##****************************************************************************80
##! 
##! macro to find the sub=directories for a given relative path
##! \qiaolx
##! 
##****************************************************************************80

MACRO(SUBDIR_LIST path_name return_list)
  FILE(GLOB new_list RELATIVE ${path_name} ${path_name}/*)
  SET(dir_list "")
  FOREACH(child ${new_list})
    IF(IS_DIRECTORY ${path_name}/${child})
      LIST(APPEND dir_list ${child})
    ENDIF()
  ENDFOREACH()
  SET(${return_list} ${dir_list})
ENDMACRO()

