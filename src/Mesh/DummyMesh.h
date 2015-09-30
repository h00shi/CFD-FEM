// -*-c++-*-
#ifndef DUMMYMESH_H
#define DUMMYMESH_H
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/List2D.h"
#include "Include/my_incl.h"
#include <limits>

//****************************************************************************80
//! \brief A namespace for creating dummy mesh streams as strings. 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
namespace DummyMesh
{

//****************************************************************************80
//! \brief Builds a hex mesh of a rectangular cubic domain of size 
//!        (0,Lx)X(0,Ly)X(0,Lz);
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] Lx Length in the x-dir
//! \param[in] Ly Length in the y-dir
//! \param[in] Lz Length in the z-dir
//! \param[in] Nx Number of points in the x-dir
//! \param[in] Ny Number of points in the y-dir
//! \param[in] Nz Number of points in the z-dir
//****************************************************************************80
  std::string SetupCubicHexMesh(const realT& Lx, const realT& Ly, 
                                const realT& Lz, const intT& Nx, 
                                const intT& Ny, const intT& Nz);
  
  
};// End Namespace DummyMesh
#endif
