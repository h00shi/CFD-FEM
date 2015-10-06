// -*-c++-*-
#ifndef UNSTMESHREADERNKBGRID_H
#define UNSTMESHREDEARNKBGRID_H
#include "my_incl.h"
#include "IO/UnstMeshReader.h"
#include "Mesh/ElementTopology.h"
//****************************************************************************80
//! \brief Reads a mesh of the format NKB-Grid
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
class UnstMeshReaderNKBGrid : public UnstMeshReader
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

//****************************************************************************80
//! \brief A method to open the mesh file.  It also reads the file into 
//!        class members.    
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  void Open(const std::string& filename);
  Array2D<realT> ReadNodes();
  // List2D<intT>  ReadElement2Node();
  // Array1D<intT> ReadElementType();
  // Array1D<intT> ReadElementRegion();
  // List2D<intT>  ReadBcFace2Node();
  // Array1D<intT> ReadBcID();
  // Array1D<intT> ReadBcFaceType();
 
private:
  bool is_open_;
  FILE* mesh_file_;
  Array2D<realT> x_;
  List2D<intT> element2node_;
  Array1D<intT> element_type_;
  Array1D<intT> element_region_;
  List2D<intT> bc_face2node_;
  Array1D<intT> bc_face_id_;
  Array1D<intT> bc_face_type_;
  
}; // End UnstMeshReaderNKBGrid
#endif
