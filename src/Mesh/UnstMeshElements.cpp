#include "Mesh/UnstMeshElements.h"
//****************************************************************************80
UnstMeshElements::UnstMeshElements(UnstMeshReader& mesh_reader)
{
  element2node_ = mesh_reader.ReadElement2Node();
  element_type_ = mesh_reader.ReadElementType();
  
  FormNode2Element();
  FormNode2Node();
  CountElementTypes();
  
}// UnstMeshElement::UnstMeshElements

//****************************************************************************80
void UnstMeshElements::FormNode2Element()
{
  //---> Local Variables

 
  intT size_node2element; // Size of the array node2element
  SystemModule::cout 
    << "-------------Node 2 Element Connectivity ---------------------"
    << std::endl;
  //---> Initialize the temporary counter array
  Array1D<intT> nodeitemp(nnode_ + 1);
  
  /*---> Loop over elements and count the number of elements attached to a
    node */
  for (intT e = 0; e < nelement_; e++) { // element_loop
    //---> Get number of nodes on this element
    intT nne = element2node_.get_ncol(e);
    //---> Loop over nodes on the element
    for(intT n = 0; n < nne; n++) { // node_loop
      //---> Get node number
      intT node = element2node_(e, n);
      /*---> Add 1 to the count of number of elements attached to a node
        for node */
      nodeitemp(node + 1) += 1;
    } // End node_loop
  } // End element_loop

    /*---> Now add up all the element surrounding each node and setup
      linked list index array */
  size_node2element = 0;
  for (intT n = 0; n < nnode_ ; n++) { // Node_loop
    //---> Size of node2element array computation
    size_node2element += nodeitemp(n + 1);
  } // End Node_loop

  //---> Allocate the memory for node2element data array
  node2element_.initialize(nnode_, size_node2element);
  for(intT n = 0; n < nnode_; n++){ node2element_.set_ncol(n,nodeitemp(n + 1));}
  
  nodeitemp.set_value(0);
  
  //---> Now for the connectivity node2element
  for (intT e = 0; e < nelement_; e++ ) { // Element_loop
    //---> Get number of nodes on this element
    intT nne = element2node_.get_ncol(e);
    //---> Loop over nodes on the element
    for (intT n = 0; n < nne; n++) { // node_loop
      //---> Get node number
      intT node = element2node_(e, n);
      //----> Assign element e to index position of node2element
      node2element_(node, nodeitemp(node)) = e;
      /*----> We've added an element to the node2element for node: node...so
        add 1 to node2elementi(node) to act as a counter for how many
        elements we've added to the node, this acts as a local index now */
      //node2elementi(node) += 1;
      nodeitemp(node) += 1;
    } // End node_loop
  } // End element_loop
 
  SystemModule::cout 
    << "------------- End Node 2 Element Connectivity -----------------"
    << std::endl;
 
  return;
}// End Function FormNode2Element
//****************************************************************************80
void UnstMeshElements::CountElementTypes()
{

  nbar_ = 0;
  ntri_ = 0;
  nquad_ = 0;
  ntet_ = 0;
  nprism_ = 0;
  npyr_ = 0;
  nhex_ = 0;

  for (intT e = 0; e < nelement_; e++) { // Count element types
    switch (element_type_(e)) {
    case ElementTopology::element_types::BAR:
      nbar_ += 1;
      break;
    case ElementTopology::element_types::TRI:
      ntri_ += 1;
      break;
    case ElementTopology::element_types::QUAD:
      nquad_ += 1;
      break;
    case ElementTopology::element_types::TET:
      ntet_ += 1;
      break;
    case ElementTopology::element_types::PRISM:
      nprism_ += 1;
      break;
    case ElementTopology::element_types::PYR:
      npyr_ += 1;
      break;
    case ElementTopology::element_types::HEX:
      nhex_ += 1;
      break;
    }
  }
}// End CountElementTypes
//****************************************************************************80
void UnstMeshElements::FormNode2Node() 
{  
  
  //---> Initialize size of list to zero
  intT nnz_adj = 0;
  //---> Loop over the nodes
  for(intT n = 0; n < nnode_; n++) { // Node_loop
    intT nn_tot = 0;
    //---> Loop over the elements that surround this node 
    for(intT i = 0; i < node2element_.get_ncol(n); i++){
      intT e = node2element_(n,i);
      nn_tot += element2node_.get_ncol(e);
    }
    
    //---> Temporary array for nodes attached to node n
    std::vector<int> temp (nn_tot,-1);
    
    //---> Initialize counter
    intT counter = 0;
    
    //---> Now loop over the elements containing that node
    for( intT i = 0; i < node2element_.get_ncol(n); i++ ) {// element_loop
      //---> Get element containing node n
      intT e = node2element_(n, i);
      
      //---> Loop over nodes attached to elem
      for ( intT j = 0; j < element2node_.get_ncol(e); j++){ // Node_on_element
	
	//---> Get node index of jth node attached to elem
	intT node = element2node_(e, j);	 
	
	//---> Store node in temp 
	temp[counter] = node;
	
	//---> Increment counter
	counter += 1;
	
      } // Node_on_element
      
    } // End element_loop 
    
    //---> Sort all possible nodes adjacent to node n in lexigraphical order
    std::sort (temp.begin(), temp.end());
    
    //---> Define variables to use unique function 
    std::vector<int>::iterator itemp;
    
    /*---> Eleminate all duplicate nodes from list temp NOTE: Must call 
      sort first for std::unique to work as expected.  */
    itemp = std::unique(temp.begin(), temp.end());
    
    //---> Resize temp to remove undefined entries
    temp.resize(std::distance( temp.begin(),itemp ) );
    
    //---> Return size of vector temp which is number of unique adjacent nodes
    intT c = temp.size();
    nnz_adj += c;

  } // End Node_loop 
  //--------------------------------------------------------------------------
  node2node_.initialize(nnode_, nnz_adj);
  mem_ += node2node_.get_mem();
  //--------------------------------------------------------------------------
  
  //---> Build Adjacency List
  //---> Loop over the nodes
  for (intT n = 0; n < nnode_; n++) { // Node_loop
    intT nn_tot = 0;
    //---> Loop over the elements that surround this node 
    for(intT i = 0; i < node2element_.get_ncol(n); i++){
      intT e = node2element_(n,i);
      nn_tot += element2node_.get_ncol(e);
    }
    
    //---> Temporary array for nodes attached to node n
    std::vector<int> temp (nn_tot,-1);
    
    //---> Initialize counter
    intT counter = 0;
    
    //---> Now loop over the elements containing that node
    for( intT i = 0; i < node2element_.get_ncol(n); i++ ) {// element_loop
      //---> Get element containing node n
      intT e = node2element_(n, i);
      
      //---> Loop over nodes attached to elem
      for ( intT j = 0; j < element2node_.get_ncol(e); j++){ // Node_on_element
	
	//---> Get node index of jth node attached to elem
	intT node = element2node_(e, j);	 
	
	//---> Store node in temp 
	temp[counter] = node;
	
	//---> Increment counter
	counter += 1;
	
      } // Node_on_element
      
    } // End element_loop 
    
    //---> Sort all possible nodes adjacent to node n in lexigraphical order
    std::sort (temp.begin(), temp.end());
    
    //---> Define variables to use unique function 
    std::vector<int>::iterator itemp;
    
    /*---> Eleminate all duplicate nodes from list temp NOTE: Must call 
      sort first for std::unique to work as expected.  */
    itemp = std::unique(temp.begin(), temp.end());
    
    //---> Resize temp to remove undefined entries
    temp.resize(std::distance( temp.begin(),itemp ) );
    
    //---> Return size of vector temp which is number of unique adjacent nodes
    intT c = temp.size();
    node2node_.set_ncol(n,c);
 
    for (intT i = 0; i < c; i++){node2node_(n,i) = temp[i];}
    
  } // End Node_loop 
  
  //graph_ = new Graph(adj);
}// End FormAdjacency
