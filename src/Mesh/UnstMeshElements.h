// -*-c++-*-
#ifndef UNST_MESH_ELEMENTS_H
#define UNST_MESH_ELEMENTS_H
#include "my_incl.h"
#include "SystemUtils/SystemModule.h"
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
  UnstMeshElements(const intT& nnode, 
                   const intT& nelement, 
                   List2D<intT>& element2node);
//****************************************************************************80
//! \brief Destructor 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  ~UnstMeshElements();
private:
  intT nelement_; //!< Number of elements in the mesh
  intT nbar_; //!< Number of 1-D bar elements in the mesh
  intT ntri_; //!< Number of 2-D tri elements in the mesh
  intT nquad_; //!< Number of 2-D quad elements in the mesh
  intT ntet_; //!< Number of 3-D tet elements in the mesh
  intT nprism_; //!< Number of 3-D prism elements in the mesh
  intT npyr_; //!< Number of 3-D pyramid elements in the mesh
  intT nhex_; //!< Number of 3-D hex elements in the mesh
  intT nnode_; //!< Number of nodes in the mesh
  realT mem_;//!< Amount of memory used by this class
  List2D<intT> element2node_; //!< Element2node connectivity
  List2D<intT> node2element_; //!< Node2lement connectivity
  List2D<intT> node2node_;//!< All nodes that share an element with a node
  Array1D<ElementTopology::element_types> element_type_; /*!< The type of each 
                                                          element */
  Array1D<realT> volume_;//!< Volume of each element in the mesh

//****************************************************************************80
//! \brief Default constructor.  DELETED
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  UnstMeshElements() = delete;
  
//****************************************************************************80
//! \brief FormNode2Element : Forms the node2element connecitivity
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  void FormNode2Element();

//****************************************************************************80
//! \brief FormNode2Node : Forms the list of node2node connectivity
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  void FormNode2Node();
  void CountElementTypes();

};// UnstMeshElements
#endif
