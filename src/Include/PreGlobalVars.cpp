//****************************************************************************80
//! \file PreGlobalVars.cpp 
//! \brief Contains variable definitions for pre-processor global Variables.  
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//****************************************************************************80
#include "my_incl.h"
#include <map>
#include "PreGlobalVars.h"
#include "system_module.h"
#include "Utils.h"

////////////////////////////////////////////////////////////////////////////////
//
//+++++++++++++++++++++++++++++++ Member Data Definitions ++++++++++++++++++++++
//
////////////////////////////////////////////////////////////////////////////////
std::string PreGlobalVars::project = "NULL"; 
std::string PreGlobalVars::meshtype = "NULL";
int PreGlobalVars:: pg = 0; 
int PreGlobalVars::nproc = 1; 
bool PreGlobalVars::mergebl = false; 
bool PreGlobalVars::extractline = false;
std::string PreGlobalVars::renumber = "None";

////////////////////////////////////////////////////////////////////////////////
//
//+++++++++++++++++++++++++++++ Member Function Definitions ++++++++++++++++++++
//
////////////////////////////////////////////////////////////////////////////////
    
//****************************************************************************80
//!
//! \brief initialize : Initializes the variables in namespace PreGlobalVars 
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] pre_global_vars_map A STL map that contains the global data 
//****************************************************************************80
void PreGlobalVars::initialize(std::map<std::string,std::string>& 
			       pre_global_vars_map)
{
   
  //---> Tell user what values of things are in the inputs
  std::cout << "Setting Global Variables:" << std::endl;
  std::cout << "PRE_GLOBAL: "<<std::endl;
  
  Utils::extract_var_string(pre_global_vars_map, "project", 
		     PreGlobalVars::project);
  Utils::extract_var_string(pre_global_vars_map, "meshtype", 
		     PreGlobalVars::meshtype);
  Utils::extract_var_int(pre_global_vars_map, "pg", 
		  PreGlobalVars::pg);
  Utils::extract_var_int(pre_global_vars_map, "nproc", 
		  PreGlobalVars::nproc);
  Utils::extract_var_bool(pre_global_vars_map, "mergebl", 
		   PreGlobalVars::mergebl);
  Utils::extract_var_bool(pre_global_vars_map, "extractline", 
		   PreGlobalVars::extractline);
  
  Utils::extract_var_string(pre_global_vars_map, "renumber", 
		     PreGlobalVars::renumber);


  
  std::cout << std::endl;
} // End Function initialize()

