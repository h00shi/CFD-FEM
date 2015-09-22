//****************************************************************************80
Graph(List2D<intT>& adj, Array2D<intT>& edge2node) :
  nnode_( adj.get_lead_size()), 
  nnz_(adj.get_total_size() ), 
  nedge_(edge2node.get_size(0)), 
  adj_(adj),
  edge2node_(edge2node),
  node2self_adj_index_(nnode_),
  edge2adj_index_(nedge_,2)
{


}
