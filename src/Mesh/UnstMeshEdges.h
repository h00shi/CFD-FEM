// -*-c++-*-
#ifndef UNSTMESHEDGES_H
#define UNSTMESHEDGES_H
#include "my_incl.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/List2D.h"
#include "Mesh/ElementTopology.h"
#include "Mesh/UnstMeshElements.h"
#include "DataStructures/EdgeSet.h"
#include <algorithm>
#include <vector>
#include <cmath>

//****************************************************************************80
//! \brief UnstMeshEdges : A class for extracting edges in a mesh
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
class UnstMeshEdges
{
public:
//****************************************************************************80
//! \brief UnstMeshEdges : A class that builds the edges of a mesh 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  UnstMeshEdges(UnstMeshElements& mesh_elements);
  
  inline const Array2D<intT>& get_edge2node()const {return edge2node_;}
  inline const List2D<intT>& get_node2edge() const {return node2edge_;}
  inline const List2D<intT>& get_edge2element() const {return edge2element_;}
  inline const List2D<intT>& get_element2edge() const {return element2edge_;}

private:
  const UnstMeshElements& mesh_elements_;//!< Reference to elements of mesh
  intT nedge_; //!< The number of edges...must be predicted
  intT nnz_node2edge_; //!< Number of non-zero entries n nnz_node2edge
  Array2D<intT> edge2node_;//!< Lists nodes attached to each edge
  List2D<intT> node2edge_;//!< Lists the edges attached to each node
  List2D<intT> edge2element_;//!< Lists elements attached to an edge
  List2D<intT> element2edge_;//!< Lists the edges attached to an element
  
//****************************************************************************80
//! \brief UnstMeshEdges : Default constructor...DELETED.  
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  UnstMeshEdges() = delete;

//****************************************************************************80
//! \brief ComputeNumberOfEdges : Computes the number of edges in a mesh
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  void ComputeNumberOfEdges();

//****************************************************************************80
//! \brief Extract : Driver function for doing extraction.
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  void Extract();

//****************************************************************************80
//! \brief Given an element extract all it's edges
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] elem The element you want to extract the edges of 
//****************************************************************************80
  template<class Topology>
 void ExtractEdgesOfElement(const intT& elem, EdgeSet& edges)
  {
    for(intT e = 0; e < Topology::nEdge; e++){
      intT n1 = mesh_elements_.get_element2node()(elem,Topology::Edges[e][0]);
      intT n2 = mesh_elements_.get_element2node()(elem,Topology::Edges[e][1]);
      //---> Inserting the two nodes returns the edge that is made of these two
      //      nodes.
      element2edge_(elem,e) = edges.InsertEdge(n1,n2);
    }

  } //End UnstMeshEdges::ExtractEdgesOfElement

};// End UnstMeshEdges
#endif
