#include "Mesh/UnstMeshElements.h"
//****************************************************************************80
UnstMeshElements::UnstMeshElements(UnstMeshReader& mesh_reader)
{
  nelement_ = mesh_reader.ReadNelement();
  nnode_ = mesh_reader.ReadNnode();

  element2node_ = mesh_reader.ReadElement2Node();
  element_type_ = mesh_reader.ReadElementType();
  
  FormNode2Element();
  CountElementTypes();
  
}// UnstMeshElement::UnstMeshElements
//****************************************************************************80
UnstMeshElements::~UnstMeshElements(){}
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

