// -*-c++-*-
#ifndef SOL_GLOBALVARS
#define SOL_GLOBALVARS

//****************************************************************************80
//! \file SolGlobalVars.h
//! \brief Contains the global Variables declarations for Solver global 
//!        Variables.
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! 
//****************************************************************************80
#include "my_incl.h"
#include <map>
#include "Utils.h"
namespace SolGlobalVars {
////////////////////////////////////////////////////////////////////////////////
//
//++++++++++++++++++++++++++++++++ Member Data Declarations ++++++++++++++++++++
//
////////////////////////////////////////////////////////////////////////////////
  extern std::string Project; /*!< Project String */
  extern std::string MeshType; /*!< Type of mesh specified in string */
  extern std::string Domain; /*!< String for 1-D, 2-D or 3-D */
  extern std::string EquationType; /*!< Equation Type String */
  extern int PG; /*!< Global Polynomial order */
  extern double Fmach; /*!< Freestream Mach number */
  extern double Alpha; /*!< Angle of attack for freestream flow */
  extern double Beta; /*!<  Angle of side slip for freestream flow */
  extern int NTS; /*!< Number of time-steps to run */
  extern std::string LinearSolver; /*!< Chooses the linear solver */
  extern std::string Preconditioner; /*!< Chooses the preconditioner */
  extern int IcType; /*!< Chooses the initial condition type */

////////////////////////////////////////////////////////////////////////////////
//
//++++++++++++++++++++++++++++++ Member Functions Declarations +++++++++++++++++
//
////////////////////////////////////////////////////////////////////////////////
  void initialize(std::map<std::string,std::string>& );

} // End namespace SolGlobalVars
#endif
