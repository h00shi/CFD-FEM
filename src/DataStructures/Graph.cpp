#include "DataStructures/Graph.h"
//****************************************************************************80
Graph::Graph(List2D<intT>& adj, Array2D<intT>& edge2node) :
  adj_(std::move(adj)),
  edge2node_(std::move(edge2node))
{
  /*---> For Robustness Set these values here because order of intializer list
    is not the most reliable thing.*/
  nnode_ = adj_.get_lead_size(); 
  nnz_   = adj_.get_total_size();
  nedge_ = edge2node_.get_size(0); 
  node2self_adj_index_.initialize(nnode_);
  edge2adj_index_.initialize(nedge_,2);
  
  FormGraphIndicies();
}
//****************************************************************************80
Graph::Graph(List2D<intT>& adj) :
  adj_(std::move(adj))
{
  /*---> For Robustness Set these values here because order of intializer list
    is not the most reliable thing.*/
  nnode_ = adj_.get_lead_size(); 
  nnz_   = adj_.get_total_size();
  nedge_ = (nnz_- nnode_)/2; 
  edge2node_.initialize(nedge_,2);
  node2self_adj_index_.initialize(nnode_);
  edge2adj_index_.initialize(nedge_,2);
  
  FormEdge2Node();
  FormGraphIndicies();
}
//****************************************************************************80
Graph::~Graph(){} 
//****************************************************************************80
void Graph::FormEdge2Node() 
{
  
  Array1D<intT> LastEdgeFromNode(nnode_);
  Array1D<intT> NextEdgeFromEdge(nedge_);
  LastEdgeFromNode.set_value(-1);
  NextEdgeFromEdge.set_value(-1);
  
  intT iedge = 0;
  for(intT n = 0; n < nnode_; n++){// Node loop 
    for(int j = 0; j < adj_.get_ncol(n); j++){// Neighbor loop 
      intT node = adj_(n,j);
      if(node == n) {
        //---> Do nothing...we don't need self adj for edges
      }
      else {
        intT tag = GetEdgeTag(n, node, iedge, LastEdgeFromNode, 
                              NextEdgeFromEdge);
        if( tag == -1){
          //---> Create a new edge
          edge2node_(iedge,0) = max(n,node);
          edge2node_(iedge,1) = min(n,node);
          ++iedge;
        }
      } 
    }// End Neighbor loop 
  }// End Node loop 
 
}// End FormEdge2Node

//****************************************************************************80
intT Graph::GetEdgeTag(const intT& node0, const intT& node1, const int& iedge,
                       Array1D<intT>& LastEdgeFromNode,
                       Array1D<intT>& NextEdgeFromEdge)
{
  intT tag = -1;
  intT index_node = max(node0, node1);
  //---> Get the last edge made by the index node
  intT edge = LastEdgeFromNode(index_node);
   if( edge == -1) {// EdgeTag
    //---> Do nothing we've already done this edge
    LastEdgeFromNode(index_node) = iedge;
  }
  else {
    bool cont = true;
    while(cont){ // Search Loop
      //---> Check to see if edge is an edge connecting node0 and node1
      intT nL = edge2node_(edge,0);
      intT nR = edge2node_(edge,1);
      if( max(node0,node1) == nL && min(node0,node1) == nR){ // Check Nodes
        //---> Edge edge does conntect node0 and node1
        tag = edge;
        cont = false;
      }
      else {
        if( NextEdgeFromEdge(edge) != -1 ) {
          edge = NextEdgeFromEdge(edge);
        }
        else {
          //---> We've determined that  edge = NextEdgeFromEdge(edge) = -1;
          cont = false;
        }
      }// End Nodes Check
    }// End Search Loop 
    /*---> If after all that searching we haven't found an edge made of node0
     and node1.  Then update that edge's next edge in the list is edge_count. */
    if( tag == -1) {
      NextEdgeFromEdge(edge) = iedge;
    }
  } // End EdgeTag
 
  return tag;
}// End GetEdgeTag
//****************************************************************************80
void Graph::FormGraphIndicies()
{
  //---> Setup Self adjacency index
  for(intT n = 0; n < nnode_; n++){// Node Loop
    for(intT j = 0; j < adj_.get_ncol(n); j++) { // Neighbor Loop 
      intT node = adj_(n,j);
      if(node == n){node2self_adj_index_(n) = j;break;}
    }// End Neighbor loop
  }
  
  //---> Setup edge adjacency index
  for(intT e = 0; e < nedge_; e++){ // Edge loop
    //---> Left and right nodes
    intT nl = edge2node_(e,0);
    intT nr = edge2node_(e,1);
    
    for(intT j = 0; j < adj_.get_ncol(nl); j++){// left node neighbor idx
      intT node = adj_(nl,j);
      if(node == nr){edge2adj_index_(e,0) = j; break;}
    }// End left node neighbor idx
    
    for(intT j = 0; j < adj_.get_ncol(nr); j++){//right node neighbor idx
      intT node = adj_(nr,j);
      if(node == nl){edge2adj_index_(e,1) = j; break;}
    }//End right node neighbor idx
    
  }// End Edge loop 


}// End FormGraphIndicies
