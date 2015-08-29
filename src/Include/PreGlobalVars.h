// -*-c++-*-
#ifndef PRE_GLOBALVARS
#define PRE_GLOBALVARS
//****************************************************************************80
// \file PreGlobalVars.h 
//! \namespace PreGlobalVars
//! \brief Contains variable declarations for pre-processor global Variables.  
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//****************************************************************************80
#include "my_incl.h"
#include <map>

namespace PreGlobalVars {
////////////////////////////////////////////////////////////////////////////////
//
//++++++++++++++++++++++++++++++++ Member Data Declarations ++++++++++++++++++++
//
////////////////////////////////////////////////////////////////////////////////
  extern std::string project; 
  extern std::string meshtype;
  extern int pg; 
  extern int nproc; 
  extern bool mergebl; 
  extern bool extractline;
  extern std::string renumber;
 
////////////////////////////////////////////////////////////////////////////////
//
//++++++++++++++++++++++++++++++ Member Functions Declarations +++++++++++++++++
//
////////////////////////////////////////////////////////////////////////////////
  void initialize(std::map<std::string,std::string>& );
  /*!< Namepspace memberfunction initialize, \refcpp{PreGlobalVars} */
 
} //End namespace PreGlobalVars

#endif
