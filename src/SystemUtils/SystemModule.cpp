
#include "my_incl.h"
#include <stdlib.h>
#include "SystemModule.h"
#include <algorithm>
#include <limits>

////////////////////////////////////////////////////////////////////////////////
//
//+++++++++++++++++++++++++++++++ Member Data Definitions ++++++++++++++++++++++
//
////////////////////////////////////////////////////////////////////////////////
SystemModule::pout SystemModule::my_pout;

////////////////////////////////////////////////////////////////////////////////
//
//++++++++++++++++++++++++++++++ Member Function Definitions +++++++++++++++++++
//
////////////////////////////////////////////////////////////////////////////////

//****************************************************************************80
//! \brief my_exit : This function is a wrapper for internal c++ exit
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//****************************************************************************80
void SystemModule::my_exit () {
#ifdef HAVE_MPI
  
#else
  std::cout << "Exiting program SystemModule::my_exit()" << std::endl;
  exit (EXIT_FAILURE);
#endif  
  
}// End my_exit

//****************************************************************************80
//!
//! \brief string_to_bool : Converts std::string to boolean
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] str The string to be converted
//! \return b The boolena to ver returned
//****************************************************************************80
bool SystemModule::string_to_bool(std::string str)
{
  //---> Return var
  bool b;
  //---> Not exactly sure how this works...solution from stack overflow. 
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);
  std::istringstream is(str);
  is >> std::boolalpha >> b;
  return(b);
} // End string_to_bool

//****************************************************************************80
//!
//! \brief pause : A proper c++ function that causes the code to pause until
//!                enter is pressed.  It's for debugging loops
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! 
//****************************************************************************80
void SystemModule::pause()
{
  //---> Print message to screen
  std::cout << "Presss ENTER to continue..." << std::flush;
  //---> Now use cin to wait for input
  std::cin.ignore( std::numeric_limits <std::streamsize> ::max(), '\n' );
} //End pause
