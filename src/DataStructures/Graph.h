//-*-c++-*-

#ifndef GRAPH_H
#define GRAPH_H
//****************************************************************************80
//! \class Graph
//! \brief Graph : A class to represent graphs and connectivity of that graph
//! \nick
//! \version $Rev: 5 $
//****************************************************************************80
#include "my_incl.h"
#include "SystemUtils/SystemModule.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/List2D.h"
class Graph
{
public:
 
//****************************************************************************80
//! \brief Graph : Constructor for this class. Note will use move constructor 
//!                of input variables...effectly deleting them. 
//! \nick
//! \version $Rev: 5 $
//! \param[in] adj The adjacency for this graph 
//! \param[in] edge2node The nodes attached to each edge
//****************************************************************************80
  Graph(List2D<intT>& adj, Array2D<intT>& edge2node);	
			
//****************************************************************************80
//! \brief Graph : Constructor for this class. Note will use move constructor 
//!                of input variables...effectly deleting them. 
//! \nick
//! \version $Rev: 5 $
//! \param[in] adj The adjacency for this graph 
//****************************************************************************80
  Graph(List2D<intT>& adj);
//****************************************************************************80
//! \brief ~Graph : Destructor for this class
//! \nick
//! \version $Rev: 5 $
//****************************************************************************80
  ~Graph();

//****************************************************************************80
//! \brief get_nnode : Returns the number of nodes in a graph
//! \nick
//! \version $Rev: 5 $
//! \return nnode_ The number of nodes in the graph
//****************************************************************************80
  inline const intT& get_nnode() const {return nnode_;}

//****************************************************************************80
//! \brief get_nnode : Returns the number of edges in a graph
//! \nick
//! \version $Rev: 5 $
//! \return nedge_ The number of edges in the graph
//****************************************************************************80
  inline const intT& get_nedge() const {return nedge_;}

//****************************************************************************80
//! \brief get_GraphAdj : Returns the adjacency of the graph
//! \nick
//! \version $Rev: 5 $
//! \return adj_ The adjacency of the graph
//****************************************************************************80
  inline const List2D<intT>& get_GraphAdj() const {return adj_;}

//****************************************************************************80
//! \brief get_GraphEdge2Node : Returns the adjacency of the graph
//! \nick
//! \version $Rev: 5 $
//! \return edge2node_ The nodes that are attached to each edge 
//****************************************************************************80
  inline const Array2D<intT>& get_GraphEdge2Node() const {return edge2node_;}
  
//****************************************************************************80
//! \brief get_GraphNode2SelfIndex : Returns a 1-D array that contains the 
//!        self neighbor index of the adj.  Example array(i) = j : adj_(i,j) = i
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  inline const Array1D<intT>& get_GraphNode2SelfIndex() const 
  { 
    return node2self_adj_index_;
  }

//****************************************************************************80
//! \brief get_GraphEdge2AdjIndex : Returns a 2-D array that is 
//!                                 j = edge2adj_index_(e,0)
//!		         	  such that adj(edge2node_(e,0),j) = 
//!				  edg2node_(e,1) and 
//!				  j = edge2adj_index(e,1) such that 
//!				  adj(edge2node(e,1),j) = edge2node(e,0)
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  inline const Array2D<intT>& get_GraphEdge2AdjIndex() const
  {
    return edge2adj_index_;
  }
//****************************************************************************80
//! \brief NeighborBegin : Returns the iterator the begining of the neighbor 
//!                        list of the specified node
//! \nick
//! \version $Rev: 5 $
//! \return n The node who's neighbors you want to access; 
//****************************************************************************8
  inline intT NeighborBegin(const intT& n)
  {
    return 0;
  }// End NeighborBegin

//****************************************************************************80
//! \brief NeighborEnd : Returns the accessor to the end of the neighbor list
//!                      of the specified node
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \return n The node who's neighbors you want to access
//****************************************************************************80
  inline intT NeighborEnd(const intT& n)
  {
    return adj_.get_ncol(n);
  } // End NeighborLast

//****************************************************************************80
//! \brief GetNeighbor : Retreive specified neighbor of specified node
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] n The node who's neighbor you want
//****************************************************************************80
  inline intT& GetNeighbor (const intT& n, const intT& j)
  {
    return adj_(n,j);
  } // End GetNeighbor

protected:
  
  List2D<intT> adj_; //!< Adjacency of nodes in the graph
  Array2D<intT> edge2node_; //!< For each edge give the nodes attached to it
  Array1D<intT> node2self_adj_index_;/*!< For each node i, 
				       j = node2self_adj_index(i), such that  
				       adj_ such that adj(i,j) = i */
  Array2D<intT> edge2adj_index_;/*!< For each edge e, j = edge2adj_index_(e,0)
				  such that adj(edge2node_(e,0),j) = 
				  edg2node_(e,1) and 
				  j = edge2adj_index(e,1) such that 
				  adj(edge2node(e,1),j) = edge2node(e,0)*/
  
  intT nnode_; //!< Number of nodes in the graph
  intT nedge_; //!< Number of edges in the graph
  intT nnz_; //!< Number of non-zero adjacency entries
private:
//****************************************************************************80
//! \brief Graph : Default Constructor for this class, deleted
//! \nick
//! \version $Rev: 5 $
//****************************************************************************80
  Graph() = delete;
//****************************************************************************80
//! \brief 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] node0 One of the nodes we want to match
//! \param[in] node1 The other node we ant to match on the edge
//! \param[in] iedge Most recently constructured edge
//! \param[out] LastEdgeFromNode Last edge created by each node in the mesh
//! \param[out] NextEdgeFromEdge Next edge created after each edge
//! \return tag A tag to decided to create a new edge tag = -1 is new edge
//****************************************************************************80
  intT GetEdgeTag(const intT& node0, const intT& node1, const int& iedge,
                  Array1D<intT>& LastEdgeFromNode,
                  Array1D<intT>& NextEdgeFromEdge);

//****************************************************************************80
//! \brief Forms the edge2node list in the case that it's not input
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  void FormEdge2Node();

//****************************************************************************80
//! \brief 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  void FormGraphIndicies();
}; //End Class Graph
#endif
