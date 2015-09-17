//-*-c++-*-
#ifndef COMMUNICATION_H
#define COMMUNICATION_H
#include "my_incl.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

//****************************************************************************80
//! \brief Communication : A namespace for wrapping communication 
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//****************************************************************************80
namespace Communication {

//****************************************************************************80
//! \brief Initialize : Initializes the communication system
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//****************************************************************************80
  void Initialize();

//****************************************************************************80
//! \brief : Finalize the communication system
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//****************************************************************************80
  void Finalize();

//****************************************************************************80
//! \brief GetCommRank : Gets the rank of the current process
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//****************************************************************************80
  intT GetCommRank();

//****************************************************************************80
//! \brief GetCommSize : Gets the number of rank of the current process
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//****************************************************************************80
  intT GetCommSize();

}// End Namespace Communicator

#endif
