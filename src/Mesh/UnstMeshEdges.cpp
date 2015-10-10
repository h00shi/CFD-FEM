#include "Mesh/UnstMeshEdges.h"
#include "DataStructures/EdgeSet.h"
//****************************************************************************80
UnstMeshEdges::UnstMeshEdges(UnstMeshElements& mesh_elements) :
mesh_elements_(mesh_elements)
{
  ComputeNumberOfEdges();
  intT nnz_elem2edge = 0;
 
  for(intT e = 0; e < mesh_elements_.get_nelement(); e++){//Element Loop
  
    switch (mesh_elements_.get_element_type()(e)) {
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
                         << "Extraction.  " << std::endl;
      SystemModule::cout << "Element: " << e << " is of type " 
                         << mesh_elements_.get_element_type()(e) << " which is not valid"
                         << "." << std::endl
                         << "Valid types are BAR=0, TRI=1, QUAD = 2, TET = 3,"
                         << " PRISM=4, PYR = 5, HEX = 6" << std::endl;
      SystemModule::my_exit();
        
    }
  }// End Element Loop  
  element2edge_.initialize(mesh_elements_.get_nelement(), nnz_elem2edge);

  for(intT e = 0; e < mesh_elements_.get_nelement(); e++){//Element Loop
  
    switch (mesh_elements_.get_element_type()(e)) {
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
      SystemModule::cout << "ERROR: Invalid Element type detected during edge "
                         << "Extraction.  " << std::endl;
      SystemModule::cout << "Element: " << e << " is of type " 
                         << mesh_elements.get_element_type()(e) << " which is not valid"
                         << "." << std::endl
                         << "Valid types are BAR=0, TRI=1, QUAD = 2, TET = 3,"
                         << " PRISM=4, PYR = 5, HEX = 6" << std::endl;
      SystemModule::my_exit();
        
    }
  }// End Element Loop */
  //---> Will populate edg2node, node2edge and elem2edge
  Extract();

  //---> Build edge2element
  Array1D<intT> nelem_per_edge(nedge_);

  for(intT e = 0; e < element2edge_.get_lead_size(); e++){
    for(intT i = 0; i < element2edge_.get_ncol(e); i++){
      intT edge = element2edge_(e,i);
      nelem_per_edge(edge)++;
    }
  }

  edge2element_.initialize(nelem_per_edge);
  nelem_per_edge.set_value(0);
  for(intT e = 0; e < element2edge_.get_lead_size(); e++){
     for(intT i = 0; i < element2edge_.get_ncol(e); i++){
       intT edge = element2edge_(e,i);
       intT j = nelem_per_edge(edge);
       edge2element_(edge,j) = e;
       nelem_per_edge(edge)++;
     }
   }

}//End UnstMeshEdges
  
//****************************************************************************80
void UnstMeshEdges::ComputeNumberOfEdges()
{
  const List2D<intT>& node2element = mesh_elements_.get_node2element();
  const List2D<intT>& element2node = mesh_elements_.get_element2node();

  intT ne = 0;// Edge counter
  //---> Loop over the nodes
  for(intT n = 0; n < node2element.get_lead_size(); n++) { // Node_loop
  
    intT nn_tot = 0;
    //---> Loop over the elements that surround this node 
    for(intT i = 0; i < node2element.get_ncol(n); i++){
      intT e = node2element(n,i);
      nn_tot += element2node.get_ncol(e);
    }
    
    //---> Temporary array for nodes attached to node n
    std::vector<int> list_of_neighbors (nn_tot,-1);
    
    //---> Initialize counter
    intT neighbor_count = 0;
    
    //---> Now loop over the elements containing that node
    for( intT i = 0; i < node2element.get_ncol(n); i++ ) {// element_loop
      //---> Get element containing node n
      intT e = node2element(n, i);
      
      //---> Loop over nodes attached to elem
      for ( intT j = 0; j < element2node.get_ncol(e); j++){ // Node_on_element
	
	//---> Get node index of jth node attached to elem
	intT node = element2node(e, j);
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
    
    /*---> Eliminate all duplicate nodes from list temp NOTE: Must call
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
  intT nnode = mesh_elements_.get_node2element().get_lead_size();
  EdgeSet edges(nnode, nedge_);
  //---> Loop over elements
  const Array1D<ElementTopology::element_types>& element_type =
      mesh_elements_.get_element_type();
  for(intT e = 0; e < mesh_elements_.get_nelement(); e++)
  {
    switch(element_type(e)){
      case ElementTopology::element_types::BAR :
        ExtractEdgesOfElement<ElementTopology::Bar>(e, edges);
        break;
      case ElementTopology::element_types::TRI :
        ExtractEdgesOfElement<ElementTopology::Triangle>(e, edges);
        break;
      case ElementTopology::element_types::QUAD :
        ExtractEdgesOfElement<ElementTopology::Quadrilateral>(e, edges);
        break;
      case ElementTopology::element_types::TET:
        ExtractEdgesOfElement<ElementTopology::Tetrahedron>(e, edges);
        break;
      case ElementTopology::element_types::PRISM :
        ExtractEdgesOfElement<ElementTopology::Prism>(e, edges);
        break;
      case ElementTopology::element_types::PYR :
        ExtractEdgesOfElement<ElementTopology::Pyramid>(e, edges);
        break;
      case ElementTopology::element_types::HEX :
        ExtractEdgesOfElement<ElementTopology::Hexahedron>(e, edges);
        break;
    }
  }
  node2edge_ = edges.FormNode2Edge();
  edge2node_ = edges.RemoveEdge2Node();
}
