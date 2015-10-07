#include "Mesh/UnstMeshEdges.h"
//****************************************************************************80
UnstMeshEdges::UnstMeshEdges(const intT& nnode, const intT& nelement, 
                             const List2D<intT>& elem2node, 
                             const List2D<intT>& node2elem,
                             const Array1D<intT>& element_type) :
  nnode_(nnode), nelement_(nelement), element2node_(elem2node), 
  node2element_(node2elem), element_type_(element_type)
{
  ComputeNumberOfEdges();

  intT nnz_elem2edge = 0;
 
  for(intT e = 0; e < nelement_; e++){//Element Loop  
  
    switch (element_type_(e)) {
    case ElementTopology::element_types::BAR :
      nnz_elem2edge += 1;
      break;
    case ElementTopology::element_types::TRI : 
      nnz_elem2edge += 3;
      break;
    case ElementTopology::element_types::QUAD :
      nnz_elem2edge += 4;
      break;
    case ElementTopology::element_types::TET : 
      nnz_elem2edge += 6;
      break;
    case ElementTopology::element_types::PRISM :
      nnz_elem2edge += 9;
      break;
    case ElementTopology::element_types::PYR :
      nnz_elem2edge += 8;
      break;
    case ElementTopology::element_types::HEX :
      nnz_elem2edge += 12;
      break;
    default :
      SystemModule::cout << "ERROR: Invalid Element type detected durning edge "
                         << "exraction.  " << std::endl;
      SystemModule::cout << "Element: " << e << " is of type " 
                         << element_type_(e) << " which is not valid"
                         << "." << std::endl
                         << "Valid types are BAR=0, TRI=1, QUAD = 2, TET = 3,"
                         << " PRISM=4, PYR = 5, HEX = 6" << std::endl;
      SystemModule::my_exit();
        
    }
  }// End Element Loop  
  //---> Initialize Internal members
  edge2node_.initialize(nedge_,2);
  node2edge_.initialize(nnode_, nnz_node2edge_);
  nelement_surr_edge_.initialize(nedge_);
  element2edge_.initialize(nelement_, nnz_elem2edge);
  
  LastEdgeFromNode_.initialize(nnode_);
  LastEdgeFromNode_.set_value(-1);
  NextEdgeFromEdge_.initialize(nedge_);
  NextEdgeFromEdge_.set_value(-1);

  for(intT e = 0; e < nelement_; e++){//Element Loop  
  
    switch (element_type_(e)) {
    case ElementTopology::element_types::BAR :
      element2edge_.set_ncol(e,1);
      break;
    case ElementTopology::element_types::TRI : 
      element2edge_.set_ncol(e,3);
      break;
    case ElementTopology::element_types::QUAD :
      element2edge_.set_ncol(e,4);
      break;
    case ElementTopology::element_types::TET : 
      element2edge_.set_ncol(e,6);
      break;
    case ElementTopology::element_types::PRISM :
      element2edge_.set_ncol(e,9);
      break;
    case ElementTopology::element_types::PYR :
      element2edge_.set_ncol(e,8);
      break;
    case ElementTopology::element_types::HEX :
      element2edge_.set_ncol(e,12);
      break;
    default :
      SystemModule::cout << "ERROR: Invalid Element type detected durning edge "
                         << "exraction.  " << std::endl;
      SystemModule::cout << "Element: " << e << " is of type " 
                         << element_type_(e) << " which is not valid"
                         << "." << std::endl
                         << "Valid types are BAR=0, TRI=1, QUAD = 2, TET = 3,"
                         << " PRISM=4, PYR = 5, HEX = 6" << std::endl;
      SystemModule::my_exit();
        
    }
  }// End Element Loop 



  Extract();

}//End UnstMeshEdges
  
//****************************************************************************80
void UnstMeshEdges::ComputeNumberOfEdges()
{
   
  intT ne = 0;// Edge counter
  //---> Loop over the nodes
  for(intT n = 0; n < nnode_; n++) { // Node_loop
  
    intT nn_tot = 0;
    //---> Loop over the elements that surround this node 
    for(intT i = 0; i < node2element_.get_ncol(n); i++){
      intT e = node2element_(n,i);
      nn_tot += element2node_.get_ncol(e);
    }
    
    //---> Temporary array for nodes attached to node n
    std::vector<int> list_of_neighbors (nn_tot,-1);
    
    //---> Initialize counter
    intT neighbor_count = 0;
    
    //---> Now loop over the elements containing that node
    for( intT i = 0; i < node2element_.get_ncol(n); i++ ) {// element_loop
      //---> Get element containing node n
      intT e = node2element_(n, i);
      
      //---> Loop over nodes attached to elem
      for ( intT j = 0; j < element2node_.get_ncol(e); j++){ // Node_on_element
	
	//---> Get node index of jth node attached to elem
	intT node = element2node_(e, j);
        //---> Store node in temp 
	list_of_neighbors[neighbor_count] = node;
        //---> Increment counter
	neighbor_count += 1;
      } // Node_on_element
      
    } // End element_loop 

    //---> Sort all possible nodes adjacent to node n in lexigraphical order
    std::sort(list_of_neighbors.begin(), list_of_neighbors.end());
    //---> Define variables to use unique function 
    std::vector<int>::iterator itemp;
    
    /*---> Eleminate all duplicate nodes from list temp NOTE: Must call 
      sort first for std::unique to work as expected.  */
    itemp = std::unique(list_of_neighbors.begin(), list_of_neighbors.end());
    //---> Note that self adjacency is in this list...subtract 1 to take this
    //     away
    //---> Resize temp to remove undefined entries
    list_of_neighbors.resize(std::distance( list_of_neighbors.begin(),itemp ) );
    ne += list_of_neighbors.size() - 1;
  }// End Node_loop

  //---> ne now contains sum over nodes of nedge_per_node*2
  nedge_ = ne/2; 
  nnz_node2edge_ = ne;

}//End ComputeNumberOfEdges
//****************************************************************************80
void UnstMeshEdges::Extract()
{
  //---> Loop over elements
  intT edge_count = 0;
  for(intT e = 0; e < nelement_; e++){// Element Loop 
    switch (element_type_(e)){
    case ElementTopology::element_types::BAR : 
      ExtractEdgesBar(e, edge_count);
      break;
    case ElementTopology::element_types::TRI :
      ExtractEdgesTri(e, edge_count);
      break;
    case ElementTopology::element_types::QUAD : 
      ExtractEdgesQuad(e, edge_count);
      break;
    case ElementTopology::element_types::TET :
      ExtractEdgesTet(e, edge_count);
      break;
    case ElementTopology::element_types::PRISM :
      ExtractEdgesPrism(e, edge_count);
      break;
    case ElementTopology::element_types::PYR : 
      ExtractEdgesPyr(e, edge_count);
      break;
    case ElementTopology::element_types::HEX : 
      ExtractEdgesHex(e, edge_count);
      break;
    }
  }// End Element Loop 

}// End UnstMeshEdges::Extract
//****************************************************************************80
void UnstMeshEdges::ExtractEdgesBar(const intT& elem, intT& edge_count)
{
  intT node0 = element2node_(elem,0);
  intT node1 = element2node_(elem,1);
   
  intT edge;
  
  edge = FormEdge(node0, node1, edge_count);
  element2edge_(elem,0) = edge;
  nelement_surr_edge_(edge) += 1;
  
  edge_count += 1;
} // UnstMeshEdges::ExtractEdgesBar
//****************************************************************************80
void UnstMeshEdges::ExtractEdgesTri(const intT& elem, intT& edge_count)
{
  intT node0 = element2node_(elem,0);
  intT node1 = element2node_(elem,1);
  intT node2 = element2node_(elem,2);
  
  intT edge;
 
  //---> Edge 0
  edge = FormEdge(node0, node1, edge_count);
  element2edge_(elem, 0) = edge;
  nelement_surr_edge_(edge) += 1;
 
  //---> Edge 1
  edge = FormEdge(node1, node2, edge_count);
  element2edge_(elem, 1) = edge;
  nelement_surr_edge_(edge) += 1;

  //---> Edge 2
  edge = FormEdge(node2, node0, edge_count);
  element2edge_(elem, 2) = edge;
  nelement_surr_edge_(edge) += 1;
}// End UnstMeshEdges::ExtractEdgeTri
//****************************************************************************80
void UnstMeshEdges::ExtractEdgesQuad(const intT& elem, intT& edge_count)
{
  intT node0 = element2node_(elem,0);
  intT node1 = element2node_(elem,1);
  intT node2 = element2node_(elem,2);
  intT node3 = element2node_(elem,3);

  intT edge;
  //---> Edge 0
  edge = FormEdge(node0, node1, edge_count);
  element2edge_(elem, 0) = edge;
  nelement_surr_edge_(edge) += 1;
  
  //---> Edge 1
  edge = FormEdge(node1, node2, edge_count);
  element2edge_(elem, 1) = edge;
  nelement_surr_edge_(edge) += 1;

  //---> Edge 2
  edge = FormEdge(node2, node3, edge_count);
  element2edge_(elem, 2) = edge;
  nelement_surr_edge_(edge) += 1;

  //---> Edge 3
  edge = FormEdge(node3, node0, edge_count);
  element2edge_(elem, 3) = edge;
  nelement_surr_edge_(edge) += 1;
  
}// End UnstMeshEdges::ExtractEdgesQuad
//****************************************************************************80
void UnstMeshEdges::ExtractEdgesTet(const intT& elem, intT& edge_count)
{
}// End UnstMeshEdges::ExtractEdgesTet(
//****************************************************************************80
void UnstMeshEdges::ExtractEdgesPrism(const intT& elem, intT& edge_count)
{
}// End UnstMeshEdges::ExtractEdgesPrism
//****************************************************************************80
void UnstMeshEdges::ExtractEdgesPyr(const intT& elem, intT& edge_count)
{
} // End UnstMeshEdges::ExtractEdgesPyr
//****************************************************************************80
void UnstMeshEdges::ExtractEdgesHex(const intT& elem, intT& edge_count)
{
}// End UnstMeshEdges::ExtractEdgesHex
//****************************************************************************80
intT UnstMeshEdges::FormEdge(const intT& node0, const intT& node1, 
                             intT& edge_count)
{
  intT next_edge = GetEdgeTag(node0, node1, edge_count);
  intT edge;
  
  if( next_edge == -1) {
    //---> Mark edge2node for a new edge;
    edge2node_(edge_count, 0) = std::max(node0, node1);
    edge2node_(edge_count, 1) = std::min(node0, node1);
    
    //---> Set the return value to be edge
    edge = edge_count;
    //---> Increment the edge count
    edge_count += 1;
  }
  else {
    edge = next_edge;
  }

  return edge;
}// End FormEdge

//****************************************************************************80
intT UnstMeshEdges::GetEdgeTag(const intT& node0, const intT& node1, 
                               const intT& edge_count)
{
  intT tag = -1;
  intT index_node = std::max(node0, node1);
  //---> Get the last edge made by the index node
  intT edge = LastEdgeFromNode_(index_node);
   if( edge == -1) {// EdgeTag
    //---> Do nothing we've already done this edge
    LastEdgeFromNode_(index_node) = edge_count;
  }
  else {
    bool cont = true;
    while(cont){ // Search Loop
      //---> Check to see if edge is an edge connecting node0 and node1
      intT nL = edge2node_(edge,0);
      intT nR = edge2node_(edge,1);
      if( std::max(node0,node1) == nL && std::min(node0,node1) == nR){ // Check Nodes
        //---> Edge edge does conntect node0 and node1
        tag = edge;
        cont = false;
      }
      else {
        if( NextEdgeFromEdge_(edge) != -1 ) {
          edge = NextEdgeFromEdge_(edge);
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
      NextEdgeFromEdge_(edge) = edge_count;
    }
  } // End EdgeTag
 
  return tag;
}// End GetEdgeTag

