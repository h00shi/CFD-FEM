// -*-c++-*-
#ifndef UNST_MESH_ELEMENTS_H
#define UNST_MESH_ELEMENTS_H
#include "my_incl.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/List2D.h"
#include "Mesh/ElementTopology.h"

//****************************************************************************80
//! \brief UnstMeshElements : A class that describes elements in an 
//!                           unstructured Mesh.
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
class UnstMeshElements
{
 public:
//****************************************************************************80
//! \brief UnstMeshElements : 
//!         Constructor that takes in element2node and number of nodes
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  UnstMeshElements(const intT& nnode, const List2D<intT>& element2node);
private:
  intT nelement_; //!< Number of elements in the mesh
  intT nnode_; //!< Number of nodes in the mesh
  List2D<intT> element2node_; //!< Element2node connectivity
  List2D<intT> node2element_; //!< Node2lement connectivity
  
};// UnstMeshElements
#endif
