// -*-c++-*-
#ifndef UNSTMESHREADERNKBGRID_H
#define UNSTMESHREDEARNKBGRID_H
#include "my_incl.h"
#include "UnstMeshReader.h"

//****************************************************************************80
//! \brief Reads a mesh of the format NKB-Grid
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
class UnstMeshReaderNKBGrid
{
public:

//****************************************************************************80
//! \brief Default constructor 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  UnstMeshReaderNKBGrid();
//****************************************************************************80
//! \brief  Default destructor
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  ~UnstMeshReaderNKBGrid();
private:
  FILE * file_ptr;
  Array2D<realT> node_coords_;
  List2D<intT> elem2node_;
  Array1D<intT> elem_type_;
  Array1D<intT> elem_region_;
  List2D<intT> bc_face2node_;
  Array1D<intT> bc_id_;
  Array1D<intT> bc_face_type_;
  
}
#endif;
