//-*-c++-*-
#ifndef COMMUNICATION_H_
#define COMMUNICATION_H_
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
  inline void Initialize()
  {
#ifdef HAVE_MPI
    MPI::Init();
#endif
    return;
  }// End Communication::Initialize
//****************************************************************************80
//! \brief : Finalize the communication system
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//****************************************************************************80
  inline void Finalize()
  {
  #ifdef HAVE_MPI
    MPI::Finalize();
  #endif
    return;
  }// End Communication::Finalize
//****************************************************************************80
//! \brief Abort the communication system
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//****************************************************************************80
  inline void Abort()
  {
#ifdef HAVE_MPI
    int flag = static_cast<int>(true);
    MPI_Initialized(&flag);
    if (flag){
      int my_rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }else{
      exit(EXIT_FAILURE);
    }
#else
    exit (EXIT_FAILURE);
#endif
  }//End Communication::Abort

//****************************************************************************80
//! \brief GetCommRank : Gets the rank of the current process
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//****************************************************************************80
  inline intT GetCommRank()
  {
#ifdef HAVE_MPI
    return MPI::COMM_WORLD.Get_rank();
#else
    return 0;
#endif
  }// End Communication::GetCommRank
//****************************************************************************80
//! \brief GetCommSize : Gets the number of rank of the current process
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//****************************************************************************80
 inline intT GetCommSize()
 {
#ifdef HAVE_MPI
   return MPI::COMM_WORLD.Get_size();
#else
   return 1;
#endif
 }// End Communication::GetCommRank
//****************************************************************************80
//! \brief AllGather : Method whereby each process gathers all specified
//!                    data from other processes
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//****************************************************************************80
inline void AllGather(const intT* send_data, intT send_count, intT* recv_data,
                      intT recv_count)
{
#ifdef HAVE_MPI
 MPI::Datatype mpi_intT =
    MPI::Datatype::Match_size(MPI_TYPECLASS_INTEGER, sizeof(intT));

MPI::COMM_WORLD.Allgather(send_data, send_count, mpi_intT,
                          recv_data, recv_count, mpi_intT);
#endif
}
inline void AllGather(const realT* send_data, intT send_count, realT* recv_data,
                      intT recv_count)
{
#ifdef HAVE_MPI
 MPI::Datatype mpi_realT =
    MPI::Datatype::Match_size(MPI_TYPECLASS_REAL, sizeof(realT));

MPI::COMM_WORLD.Allgather(send_data, send_count, mpi_realT,
                     recv_data, recv_count, mpi_realT);
#endif
}

};// End Namespace Communicator

#endif
