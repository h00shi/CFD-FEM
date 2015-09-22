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
}
