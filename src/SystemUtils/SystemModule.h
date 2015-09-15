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
 
  class pout; /*!< Namespace member class pout for parallel output */
  extern SystemModule::pout my_pout;
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

} //End namespace 

//****************************************************************************80
//!
//! \brief A class to enable parallel cout type output trivially
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
class SystemModule::pout : public std::ostream
{
public:

//****************************************************************************80
//!
//! \brief pout : Default constructor. 
//! \details Default constructor will set output_rank to 0
//! \nick
//! \version $Rev$
//****************************************************************************80
  pout()
  {
    output_rank_ = 0;
  }// end pout

//****************************************************************************80
//!
//! \brief pout : Constructor that sets value of output rank for this 
//!               instance of the class 
//! \details Sets internal output_rank_
//! \nick
//! \version $Rev$
//! \param[in] output_rank The rank for outputting data
//****************************************************************************80
  pout(const intT& output_rank)
  {
    output_rank_ = output_rank;
  }// end pout

//****************************************************************************80
//!
//! \brief operator << : '<<' operator for wrapping std::cout and only printing 
//!                      if rank == specifed value
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
  template<class T>
  inline SystemModule::pout& operator << (const T& output) 
  {
    std::cout << output;
    
    return *this;
  }// End operator << 
//****************************************************************************80
//!
//! \brief operator << : Special instance of operator to handle endl
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
  inline SystemModule::pout& operator 
  << (std::ostream&(*f)(std::ostream&))
  {
    
    if(f == (std::basic_ostream<char>& (*)(std::basic_ostream<char>&)) 
       std::endl ) { // Check that input is endl;
      std::cout << std::endl;
    } // End Check that input is endl           
    
    return *this;
  } // End operator << 
   
private:
  intT output_rank_;

}; //End class pout

//#define pout pout()
#endif
