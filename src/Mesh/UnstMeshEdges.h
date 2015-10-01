// -*-c++-*-
#ifndef UNSTMESHEDGES_H
#define UNSTMESHEDGES_H
#include "my_incl.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/List2D.h"
#include "Mesh/ElementTopology.h"
#include <algorithm>
#include <vector>

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
  UnstMeshEdges(const intT& nnode, const intT& nelement, 
                const List2D<intT>& elem2node, 
                const List2D<intT>& node2elem,
                const Array1D<intT>& element_type);
  
  inline const Array2D<intT>& get_edge2node()const {return edge2node_;}
  inline const List2D<intT>& get_node2edge() const {return node2edge_;}
  inline const List2D<intT>& get_edge2element() const {return edge2element_;}
  inline const List2D<intT>& get_element2edge() const {return element2edge_;}

private:
  const intT& nnode_; //!< Number of nodes from which to build edges
  const intT& nelement_; /*!< Number of elements from which to build edges */
  const List2D<intT>& element2node_; /*!< Reference element2node 
                                       connectivity. */
  const List2D<intT>& node2element_; /*!< Reference to node2lement 
                                        connectivity. */
  const Array1D<intT>& element_type_;//!< Types of elements

  intT nedge_; //!< The number of edges...must be predicted
  intT nnz_node2edge_; //!< Number of non-zero entries n nnz_node2edge
  Array2D<intT> edge2node_;//!< Lists nodes attached to each edge
  List2D<intT> node2edge_;//!< Lists the edges attached to each node
  List2D<intT> edge2element_;//!< Lists elements attached to an edge
  List2D<intT> element2edge_;//!< Lists the edges attached to an element
  Array1D<intT> nelement_surr_edge_;//!< Number of elements around and edge
  Array1D<intT> LastEdgeFromNode_; //!< Last edge created by a node
  Array1D<intT> NextEdgeFromEdge_; /*!< Value is the next edge created after 
                                    current edge.*/
  
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
//! \brief ExtractEdgesBar : Extracts single edge of a bar element
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] elem The element you want to extract the edges of 
//****************************************************************************80
  void ExtractEdgesBar(const intT& elem, intT& edge_count);

//****************************************************************************80
//! \brief ExtractEdgesTri : Extracts edges of a triangle
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] elem The element you want to extract the edges of 
//****************************************************************************80
  void ExtractEdgesTri(const intT& elem, intT& edge_count);

//****************************************************************************80
//! \brief ExtractEdgesQuad : Extracts edges of a Qaud
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] elem The element you want to extract the edges of 
//****************************************************************************80
  void ExtractEdgesQuad(const intT& elem, intT& edge_count);

//****************************************************************************80
//! \brief ExtractEdgesTet : Extracts edges of a Tet
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] elem The element you want to extract the edges of 
//****************************************************************************80
  void ExtractEdgesTet(const intT& elem, intT& edge_count);

//****************************************************************************80
//! \brief ExtractEdgesPrism : Extracts of a Prism
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] elem The element you want to extract the edges of 
//****************************************************************************80
  void ExtractEdgesPrism(const intT& elem, intT& edge_count);

//****************************************************************************80
//! \brief ExtractEdgesPyr: Extracts edges of a Pyramid
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] elem The element you want to extract the edges of 
//****************************************************************************80
  void ExtractEdgesPyr(const intT& elem, intT& edge_count);

//****************************************************************************80
//! \brief ExtractEdgesHex : Extracts edges of a Hex
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] elem The element you want to extract the edges of 
//****************************************************************************80
  void ExtractEdgesHex(const intT& elem, intT& edge_count);

//****************************************************************************80
//! \brief FormEdge : Forms an edge from two nodes
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] node0 One of the nodes on an edge
//! \param[in] node1 The other node on the edge
//! \param[out] edge_count Keeps count of the edges returns by reference
//! \return edge The edge that was found by extraction
//****************************************************************************80
  intT FormEdge(const intT& node0, const intT& node1, intT& edge_count);
//****************************************************************************80
//! \brief GetEdgeTag : Check to see if the two specified nodes 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] node0 One of the nodes we want to match
//! \param[in] node1 The other node we ant to match on the edge
//! \param[in] edge_count The number of edges currently extracted
//! \return tag A tag to decided to create a new edge tag = -1 is new edge
//****************************************************************************80
  intT GetEdgeTag(const intT& node0, const intT& node1, const intT& edge_count);
};// End UnstMeshEdges
#endif
