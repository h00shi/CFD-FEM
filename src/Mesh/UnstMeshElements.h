// -*-c++-*-
#ifndef UNST_MESH_ELEMENTS_H
#define UNST_MESH_ELEMENTS_H
#include "my_incl.h"
#include "SystemUtils/SystemModule.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/List2D.h"
#include "Mesh/ElementTopology.h"
#include "IO/UnstMeshReader.h"
#include "IO/UnstMeshReaderNKBGrid.h"
#include <vector>
#include <algorithm>
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
  UnstMeshElements(UnstMeshReader& mesh_reader);
//****************************************************************************80
//! \brief Destructor 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  ~UnstMeshElements();

  inline intT get_nelement() const{return nelement_;}
  inline intT get_nbar() const {return(nbar_);}
  inline intT get_ntri() const {return(ntri_);}
  inline intT get_nquad() const {return(nquad_);}
  inline intT get_ntet() const {return(ntet_);}
  inline intT get_nprism() const {return(nprism_);}
  inline intT get_npyr() const {return(npyr_);}
  inline intT get_nhex() const {return(nhex_);}
  inline const List2D<intT>& get_element2node() const {return element2node_;}
  inline const List2D<intT>& get_node2element() const {return node2element_;}
  inline const Array1D<ElementTopology::element_types>& 
  get_element_type() const {return element_type_;}

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
//! \brief Count the types of elements in the mesh
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  void CountElementTypes();
  
};// UnstMeshElements
#endif
