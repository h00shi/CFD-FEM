// -*-c++-*-
#ifndef UNSTMESHREADER_H
#define UNSTMESHREADER_H
#include "my_incl.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/List2D.h"
#include "Mesh/ElementTopology.h"

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
  UnstMeshReader(const std::string& filename);

//****************************************************************************80
//! \brief Destructor
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  virtual ~UnstMeshReader() = 0;

//****************************************************************************80
//! \brief 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  inline intT ReadNnode() const {return nnode_;}
  inline intT ReadNelement()const {return nelement_;}
  inline intT ReadNbcFace() const {return nbc_face_;}
  inline intT ReadNbcID() const {return nbc_id_;}
  virtual Array2D<realT> ReadNodes() = 0;
  virtual List2D<intT>  ReadElement2Node() = 0;
  virtual Array1D<ElementTopology::element_types> ReadElementType() = 0;
  virtual Array1D<intT> ReadElementRegion() = 0;
  virtual List2D<intT>  ReadBcFace2Node() = 0;
  virtual Array1D<intT> ReadBcID() = 0;
  virtual Array1D<ElementTopology::face_types> ReadBcFaceType() = 0;

protected:
  std::string filename_;
  intT nnode_;
  intT nelement_;
  intT nbc_face_;
  intT nbc_id_;
};// End UnstMeshReader
#endif
