#include "SystemUtils/SystemModule.h"

////////////////////////////////////////////////////////////////////////////////
//
//+++++++++++++++++++++++++++++++ Member Data Definitions ++++++++++++++++++++++
//
////////////////////////////////////////////////////////////////////////////////
SystemModule::pout SystemModule::cout;

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


