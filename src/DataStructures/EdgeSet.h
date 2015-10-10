//-*-c++-*-
#ifndef DCEL_H
#define DCEL_H
#include "my_incl.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/List2D.h"
#include <cmath>
#include <algorithm>
//****************************************************************************80
//! \brief A class to represent edges as pairs of nodes in a unique fashion.   
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
class EdgeSet{
public:
//****************************************************************************80
//! \brief Constructor the EdgeSet class.  Requires that number of nodes and
//!        edges are known at construction time.  
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] nnode The number of nodes that build up the edge set
//! \param[in] nedge The number of unique edges in the set. 
//****************************************************************************80
  EdgeSet(const intT& nnode, const intT& nedge);

//****************************************************************************80
//! \brief Takes two nodes that define and edge and inserts these into the set
//!        such that no two edges are duplicated.  
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] node0 First node on the edge
//! \param[in] node1 Second node on the edge
//! \return The unique edge made out of the two input nodes
//****************************************************************************80
  intT InsertEdge(const intT& node0, const intT& node1);
  
//****************************************************************************80
//! \brief Uses the move constructor to remove the edge2node array from the 
//!        class.  Class will no longer function this way and must be deleted. 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  inline Array2D<intT> RemoveEdge2Node() {return std::move(edge2node_);}

//****************************************************************************80
//! \brief FormNode2Edge Given that the edge2node array has not yet been removed
//!        we form the inverse map node2edge and return it to user.  
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
 List2D<intT> FormNode2Edge();

 //***************************************************************************80
 //! \brief Form node adjacency, which is a list of all nodes surrounding a node
 //! including the node itself. This list is sorted by the node numbers of all
 //! neighboring nodes.
 //! \details
 //! \nick
 //! \version $Rev$
 //! \date $Date$
 //!
 //****************************************************************************80
 List2D<intT> FormNodeAdj();
 inline intT get_nedge() const {return nedge_;}
private:
  intT nnode_; //!< Number of nodes in set
  intT nedge_; //!< Number of edges in set
  intT iedge_; //!< Counter of current edge being built/inserted
  Array1D<intT> last_edge_from_node_; /*!< Most recently created edge for 
                                        each node. */
  Array1D<intT> next_edge_from_edge_; /*!< Next edge in creation sequence from 
                                        current edge. */
  Array2D<intT> edge2node_;//!< Edge2node this is the product of the list.

//****************************************************************************80
//! \brief GetEdgeTag : Check to see if the two specified nodes 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] node0 One of the nodes we want to match
//! \param[in] node1 The other node we ant to match on the edge
//! \return tag A tag to decided to create a new edge tag = -1 is new edge
//****************************************************************************80
  intT GetEdgeTag(const intT& node0, const intT& node1);

  EdgeSet() = delete;

};// End EdgeSet

#endif
