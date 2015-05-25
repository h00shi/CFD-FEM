//****************************************************************************80
//! \file SolGlobalVars.cpp 
//! \brief This file contains implementation for dealing with SolGlobalVars
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! 
//****************************************************************************80
#include "my_incl.h"
#include "SolGlobalVars.h"

////////////////////////////////////////////////////////////////////////////////
//
//+++++++++++++++++++++++++++++++ Member Data Definitions ++++++++++++++++++++++
//
////////////////////////////////////////////////////////////////////////////////
std::string SolGlobalVars::Project = "NULL"; 
std::string SolGlobalVars::MeshType = "NULL";
std::string SolGlobalVars::Domain="NULL";
std::string SolGlobalVars::EquationType="NULL";
int SolGlobalVars:: PG = 0; 
double SolGlobalVars::Fmach;
double SolGlobalVars::Alpha; 
double SolGlobalVars::Beta;
int SolGlobalVars::NTS; 
std::string SolGlobalVars::LinearSolver; 
std::string SolGlobalVars::Preconditioner; 
int SolGlobalVars::IcType; 


//****************************************************************************80
//!
//! \brief initialize : Initializes the variables in namespace PreGlobalVars 
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] sol_global_vars_map A STL map that contains the global data 
//****************************************************************************80
void SolGlobalVars::initialize(std::map<std::string,std::string>& 
			       sol_global_vars_map)

{
  //---> Tell user what values of things are in the inputs
  std::cout << "Setting Global Variables:" << std::endl;
  std::cout << "GLOBAL: "<<std::endl;
  
  Utils::extract_var_string(sol_global_vars_map, "Project", 
		     SolGlobalVars::Project);
  Utils::extract_var_string(sol_global_vars_map, "MeshType", 
		     SolGlobalVars::MeshType);
  Utils::extract_var_string(sol_global_vars_map, "Domain", 
			    SolGlobalVars::Domain);
  Utils::extract_var_string(sol_global_vars_map, "EquationType", 
			    SolGlobalVars::EquationType);
  Utils::extract_var_int(sol_global_vars_map, "PG", 
		  SolGlobalVars::PG);
  Utils::extract_var_real(sol_global_vars_map, "Fmach",
			  SolGlobalVars::Fmach);
  Utils::extract_var_real(sol_global_vars_map, "Alpha",
			  SolGlobalVars::Alpha);
  Utils::extract_var_real(sol_global_vars_map, "Beta",
			  SolGlobalVars::Beta);
  Utils::extract_var_int(sol_global_vars_map, "NTS",
			  SolGlobalVars::NTS);

  
  std::cout << "SOLVER: "<<std::endl;

  Utils::extract_var_string(sol_global_vars_map, "LinearSolver", 
			    SolGlobalVars::LinearSolver);
  Utils::extract_var_string(sol_global_vars_map, "Preconditioner", 
			    SolGlobalVars::Preconditioner);
  
 
  std::cout << "INITIAL_CONDITIONS: "<<std::endl;
  
  Utils::extract_var_int(sol_global_vars_map, "IcType",
			 SolGlobalVars::IcType);
  std::cout << std::endl;
  std::cout << std::endl;
}
