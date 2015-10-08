#include "DataStructures/EdgeSet.h"
//****************************************************************************80
EdgeSet::EdgeSet(const intT& nnode, const intT& nedge) : nnode_(nnode), 
                                                         nedge_(nedge)
{
  last_edge_from_node_.initialize(nnode_);
  next_edge_from_edge_.initialize(nedge_);
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
        //---> Edge edge does conntect node0 and node1
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
