#include "ParallelComm/Communication.h"

//****************************************************************************80
void Communication::Initialize()
{
#ifdef HAVE_MPI
  MPI::Init(NULL, NULL);
#endif
  return;
}// End Communication::Initialize

//****************************************************************************80
void Communication::Finalize()
{
#ifdef HAVE_MPI
  MPI::Finalize();
#endif
  return;
}// End Communication::Finalize

//****************************************************************************80
intT Communication::GetCommRank()
{
#ifdef HAVE_MPI
  return MPI::COMM_WORLD.Get_rank();
#else
  return 0;
#endif
}// End Communication::GetCommRank

//****************************************************************************80
intT Communication::GetCommSize()
{
#ifdef HAVE_MPI
  return MPI::COMM_WORLD.Get_size();
#else
  return 1;
#endif
}// End Communication::GetCommRank
