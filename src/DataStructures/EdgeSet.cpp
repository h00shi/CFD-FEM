#include "DataStructures/EdgeSet.h"
//****************************************************************************80
EdgeSet::EdgeSet(const intT& nnode, const intT& nedge) : nnode_(nnode), 
                                                         nedge_(nedge)
{
  last_edge_from_node_.initialize(nnode_);
  last_edge_from_node_.set_value(-1);
  
  next_edge_from_edge_.initialize(nedge_);
  next_edge_from_edge_.set_value(-1);
  edge2node_.initialize(nedge_, 2);
  iedge_ = 0;
}// End EdgeSet::EdgeSet
//****************************************************************************80
void EdgeSet::InsertEdge(const intT& node0, const intT& node1)
{
  intT tag = GetEdgeTag(node0, node1);
  if(tag == -1) {
    edge2node_(iedge_,0) = std::max(node0, node1);
    edge2node_(iedge_,1) = std::min(node0, node1);
    iedge_ += 1;
  }
  return;
}// End EdgeSet::InsertEdge
//****************************************************************************80
intT EdgeSet::GetEdgeTag(const intT& node0, const intT& node1)
{
  intT tag = -1;
  intT index_node = std::max(node0, node1);
  //---> Get the last edge made by the index node
  intT edge = last_edge_from_node_(index_node);
   if( edge == -1) {// EdgeTag
    //---> Do nothing we've already done this edge
    last_edge_from_node_(index_node) = iedge_;
  }
  else {
    bool cont = true;
    while(cont){ // Search Loop
      //---> Check to see if edge is an edge connecting node0 and node1
      intT nL = edge2node_(edge,0);
      intT nR = edge2node_(edge,1);
      if(std::max(node0,node1) == nL && std::min(node0,node1) == nR){ // Check Nodes
        //---> Edge edge does connect node0 and node1
        tag = edge;
        cont = false;
      }
      else {
        if( next_edge_from_edge_(edge) != -1 ) {
          edge = next_edge_from_edge_(edge);
        }
        else {
          //---> We've determined that  edge = NextEdgeFromEdge(edge) = -1;
          cont = false;
        }
      }// End Nodes Check
    }// End Search Loop 
    /*---> If after all that searching we haven't found an edge made of node0
     and node1.  Then update that edge's next edge in the list is iedge_. */
    if( tag == -1) {
      next_edge_from_edge_(edge) = iedge_;
    }
  } // End EdgeTag
 
  return tag;
}
//****************************************************************************80
List2D<intT> EdgeSet::FormNode2Edge()
{
  Array1D<intT> nedge_per_node(nnode_);
  for(intT i = 0; i < nedge_; i++){
    intT nL = edge2node_(i,0);
    intT nR = edge2node_(i,1);
    
    nedge_per_node(nL) += 1;
    nedge_per_node(nR) += 1;
  }
  
  intT nnz = 0;
  for(intT i = 0; i < nnode_; i++){
	  nnz += nedge_per_node(i);
  }
  //---> Instantiate
  List2D<intT> node2edge(nnode_, nnz);
  //---> Set number of columns
  node2edge.set_ncol(nedge_per_node);
  //---> Reset to zero to be used as helper/counter
  nedge_per_node.set_value(0);
  for(intT i = 0; i < nedge_; i++)
  {
	  //---> Get Nodes left and right
	  intT nL = edge2node_(i,0);
	  intT nR = edge2node_(i,1);

	  //---> Inform each node that edge i is incident on them.
	  node2edge(nL, nedge_per_node(nL)) = i;
	  node2edge(nR, nedge_per_node(nR)) = i;

	  //---> Increment local counters
	  nedge_per_node(nL) += 1;
	  nedge_per_node(nR) += 1;

  }
  return node2edge;
}// End EdgeSet::FormNode2Edge()
//****************************************************************************80
List2D<intT> EdgeSet::FormNodeAdj()
{
  Array1D<intT> nnode_per_node(nnode_);
  nnode_per_node.set_value(1); //---> Set adjacency
  for(intT i = 0; i < nedge_; i++){
    intT nL = edge2node_(i,0);
    intT nR = edge2node_(i,1);

    nnode_per_node(nL) += 1;
    nnode_per_node(nR) += 1;
  }

  intT nnz = 0;
  for(intT i = 0; i < nnode_; i++){
	  nnz += nnode_per_node(i);
  }
  //---> Instantiate
  List2D<intT> adj(nnode_, nnz);
  //---> Set number of columns
  adj.set_ncol(nnode_per_node);
  //---> Reset to zero to be used as helper/counter
  nnode_per_node.set_value(1);
  //---> Set self adjacency first
  for(intT i = 0; i < nnode_; i++){adj(i,0) = i;}

  for(intT i = 0; i < nedge_; i++)
  {
	  //---> Get Nodes left and right
	  intT nL = edge2node_(i,0);
	  intT nR = edge2node_(i,1);

	  //---> Inform each node that edge i is incident on them.
	  adj(nL, nnode_per_node(nL)) = nR;
	  adj(nR, nnode_per_node(nR)) = nL;

	  //---> Increment local counters
	  nnode_per_node(nL) += 1;
	  nnode_per_node(nR) += 1;

  }

  for(intT i = 0; i < nnode_; i++){
	  intT ncol = adj.get_ncol(i);
	  std::sort(adj.get_ptr(i,0), adj.get_ptr(i,0) + ncol);
  }
  //---> Now loop and sort
  return adj;
}// End EdgeSet::FormNode2Edge()
