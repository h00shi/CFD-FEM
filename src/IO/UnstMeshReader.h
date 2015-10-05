// -*-c++-*-
#ifndef UNSTMESHREADER_H
#define UNSTMESHREDEAR_H
#include "my_incl.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/List2D.h"

//****************************************************************************80
//! \brief This is the base class for reading unstructure meshes.  All classes
//!        for this task must inherit from this one. 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
class UnstMeshReader
{
public:
  
//****************************************************************************80
//! \brief   Default constructor, there is nothing to construct. 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//****************************************************************************80
  UnstMeshReader();

//****************************************************************************80
//! \brief Destructor
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  virtual ~UnstMeshReader();

//****************************************************************************80
//! \brief 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  virtual Open(const std::string& filename) = 0;
  virtual Array2D<realT> ReadNodes() = 0;
  virtual List2D<intT>  ReadElement2Node() = 0;
  virtual Array1D<intT> ReadElementType() = 0;
  virtual Array1D<intT> ElementRegion() = 0;
  virtual List2D<intT>  ReadBcFace2Node() = 0;
  virtual Array1D<intT> ReadBcID() = 0;
  virtual Array1D<intT> ReadBcFaceType() = 0;

private:
 
};// End UnstMeshReader
#endif
