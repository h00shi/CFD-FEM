##****************************************************************************80
##! 
##! macro to recursively find the sub=directories w/ header files for a given
##! relative path
##! 
##****************************************************************************80

MACRO(HEADER_DIRECTORIES path_name return_list)
  FILE(GLOB_RECURSE new_list ${path_name}/*.h)
  SET(dir_list "")
  FOREACH(file_path ${new_list})
    GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
    SET(dir_list ${dir_list} ${dir_path})
  ENDFOREACH()
  LIST(REMOVE_DUPLICATES dir_list)
  SET(${return_list} ${dir_list})
ENDMACRO()

