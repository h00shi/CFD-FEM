#include "SystemUtils/SystemModule.h"
#include "ParallelComm/Communication.h"
#include <stdlib.h>
#include <execinfo.h>
#include <cxxabi.h>
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

//****************************************(FILE *out = stderr, unsigned int max_frames = 63************************************80
//! \brief my_exit : This function is a wrapper for internal c++ exit
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//****************************************************************************80
void SystemModule::my_exit () {
  if(Communication::GetCommRank()== 0){SystemModule::PrintStacktrace();}
  SystemModule::cout << "Exiting program SystemModule::my_exit()" << std::endl;
  Communication::Abort();
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

//****************************************************************************80
void SystemModule::PrintStacktrace(FILE *out, unsigned int max_frames)
{
  fprintf(out, "STACK TRACE:\n");

  // storage array for stack trace address data
  void* addrlist[max_frames+1];

  // retrieve current stack addresses
  int addrlen = backtrace(addrlist, sizeof(addrlist) / sizeof(void*));

  if (addrlen == 0) {
    fprintf(out, "  <empty, possibly corrupt>\n");
    return;
  }

  // resolve addresses into strings containing "filename(function+address)",
  // this array must be free()-ed
  char** symbollist = backtrace_symbols(addrlist, addrlen);

  // allocate string which will be filled with the demangled function name
  size_t funcnamesize = 256;
  char* funcname = (char*)malloc(funcnamesize);

  // iterate over the returned symbol lines. skip the first, it is the
  // address of this function.
  for (int i = 1; i < addrlen; i++)
  {
    char *begin_name = 0, *begin_offset = 0, *end_offset = 0;

    // find parentheses and +address offset surrounding the mangled name:
    // ./module(function+0x15c) [0x8048a6d]
    for (char *p = symbollist[i]; *p; ++p)
    {
      if (*p == '(')
        begin_name = p;
      else if (*p == '+')
        begin_offset = p;
      else if (*p == ')' && begin_offset) {
        end_offset = p;
        break;
      }
    }

    if (begin_name && begin_offset && end_offset
        && begin_name < begin_offset)
    {
      *begin_name++ = '\0';
      *begin_offset++ = '\0';
      *end_offset = '\0';

      // mangled name is now in [begin_name, begin_offset) and caller
      // offset in [begin_offset, end_offset). now apply
      // __cxa_demangle():

      int status;
      char* ret = abi::__cxa_demangle(begin_name,
                                      funcname, &funcnamesize, &status);
      if (status == 0) {
        funcname = ret; // use possibly realloc()-ed string
        fprintf(out, "  %s : %s+%s\n",
                symbollist[i], funcname, begin_offset);
      }
      else {
        // demangling failed. Output function name as a C function with
        // no arguments.
        fprintf(out, "  %s : %s()+%s\n",
                symbollist[i], begin_name, begin_offset);
      }
    }
    else
    {
      // couldn't parse the line? print the whole line.
      fprintf(out, "  %s\n", symbollist[i]);
    }
  }

  free(funcname);
  free(symbollist);
}

