//-*-c++-*-
#ifndef DCEL_H
#define DCEL_H
#include "my_incl.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"

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
//****************************************************************************80
  void InsertEdge(const intT& node0, const intT& node1);
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
