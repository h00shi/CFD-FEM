//-*-c++-*-
#ifndef SYS_MODULE
#define SYS_MODULE
//****************************************************************************80
//! \file SystemModule.h 
//! \namespace SystemModule
//! \brief Defines namespace called system module.  
//! \details system module is a namespace where I wrap certain system functions
//!          so that inclusion of MPI is more seemless.  Things like exit and
//!          cout.  This file header defines interfaces to these functions 
//!          as well as namepsace.  
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//****************************************************************************80
#include <stdlib.h>
#include "my_incl.h"

namespace SystemModule {
////////////////////////////////////////////////////////////////////////////////
//
//++++++++++++++++++++++++++++++++ Member Data Declarations ++++++++++++++++++++
//
////////////////////////////////////////////////////////////////////////////////
  extern std::stringstream my_cout;

////////////////////////////////////////////////////////////////////////////////
//
//++++++++++++++++++++++++++++ Member Function Declarations ++++++++++++++++++++
//
////////////////////////////////////////////////////////////////////////////////
  
  void my_exit(); /*!< Namespace member function my_exit(), 
		   \refcpp{SystemModule} */
  void writeout();/*!< Namespace member function writeout(), 
		    \refcpp{SystemModule} */ 
  bool string_to_bool(std::string); /*!< Namespace memeber function 
				      string_to_bool, \refcpp{SystemModule} */
  void pause(); /*!< Namespace member function pause, \refcpp{SystemModule} */
 
//****************************************************************************80
//!
//! \brief  alloc_mem : Memory allocation wrapper function 
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \tparam dataT Datatype being allocated
//! \tparam intT Integer data type
//! \tparam realT Real data type
//! \param[out] data The data to be allocated 
//! \param[in] size Size to be allocated
/// \return mem Amount of memory in megabites for this array
//****************************************************************************80
  template<class dataT, class intT, class realT>
  realT alloc_mem(dataT*& data, const intT size);

//****************************************************************************80
  //---> Implementation
  template<typename dataT, typename intT, typename realT>
  realT alloc_mem(dataT*& data, const intT size )
  {
    realT mem;
    intT bytes = sizeof(dataT);
    
    //---> Allocate the memeory
    data = new dataT[size];
   
    mem = bytes*((realT) size)/1000000.0;
    return(mem);
  }

//****************************************************************************80
//! \brief  VarMemoryPrint : For a specified variable print it's size  and 
//!         memory used.  
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] var_name The name of the variable
//! \param[in] var The variable to print the diago
//****************************************************************************80

} //End namespace SystemModule
#endif
