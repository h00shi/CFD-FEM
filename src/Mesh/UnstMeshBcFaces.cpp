/*
 * UnstMeshBcFaces.cpp
 *
 *  Created on: Oct 9, 2015
 *      Author: rabbit
 */
#include "Mesh/UnstMeshBcFaces.h"
//****************************************************************************80
UnstMeshBcFaces::UnstMeshBcFaces(UnstMeshReader& mesh_reader)
{
  //---> Read from mesh file
  nbc_face_ = mesh_reader.ReadNbcFace();
  nbc_id_ = mesh_reader.ReadNbcID();
  nbc_node_ = 0;
  nbc_bar_ = 0;
  nbc_tri_ = 0;
  nbc_quad_ = 0;

  bc_face2node_ = mesh_reader.ReadBcFace2Node();
  bc_face_id_   = mesh_reader.ReadBcID();
  bc_face_type_ = mesh_reader.ReadBcFaceType();

  //---> Initialize Internal members;
  nface_per_bcid_.initialize(nbc_id_);

  for(intT f = 0; f < nbc_face_; f++){
    intT id = bc_face_id_(f);
    nface_per_bcid_(id) += 1;

    switch(bc_face_type_(f)) {
      case ElementTopology::face_types::FACE_NODE :
        nbc_node_ += 1;
        break;
      case ElementTopology::face_types::FACE_BAR :
        nbc_bar_ += 1;
        break;
      case ElementTopology::face_types::FACE_TRI :
        nbc_tri_ += 1;
        break;
      case ElementTopology::face_types::FACE_QUAD :
        nbc_quad_ += 1;
        break;
    }
  }

}// end UnstMeshBcFaces::UnstMeshBcFaces
//****************************************************************************80
UnstMeshBcFaces::UnstMeshBcFaces(UnstMeshReader& mesh_reader,
                                 UnstMeshElements& mesh_elements) :
  UnstMeshBcFaces(mesh_reader)
{
 FormBcFace2Element(mesh_elements);
}
//****************************************************************************80
UnstMeshBcFaces::~UnstMeshBcFaces(){}
//****************************************************************************80
void UnstMeshBcFaces::FormBcFace2Element(UnstMeshElements& mesh_elements)
{
  bc_face2elem_.initialize(nbc_face_);
  bc_local_face_.initialize(nbc_face_);
  for(intT f = 0; f < nbc_face_; f++){
   bc_face2elem_(f) = FindBcElement(f, mesh_elements.get_node2element());
  }

  //---> Now form the local face numbers;
  for(intT f = 0; f < nbc_face_; f++){
    intT e = bc_face2elem_(f);

    switch(mesh_elements.get_element_type()(e)){
      case ElementTopology::element_types::BAR :
        bc_local_face_(f) = ExtractLocalFace
          <ElementTopology::Bar>(f, mesh_elements.get_element2node());
        break;
      case ElementTopology::element_types::TRI :
        bc_local_face_(f) = ExtractLocalFace
          <ElementTopology::Triangle>(f, mesh_elements.get_element2node());
        break;
      case ElementTopology::element_types::QUAD :
        bc_local_face_(f) = ExtractLocalFace
          <ElementTopology::Quadrilateral>(f, mesh_elements.get_element2node());
        break;
      case ElementTopology::element_types::TET :
        bc_local_face_(f) = ExtractLocalFace
          <ElementTopology::Tetrahedron>(f, mesh_elements.get_element2node());
        break;
      case ElementTopology::element_types::PRISM :
        bc_local_face_(f) = ExtractLocalFace
          <ElementTopology::Prism>(f,mesh_elements.get_element2node());
        break;
      case ElementTopology::element_types::PYR :
        bc_local_face_(f) = ExtractLocalFace
          <ElementTopology::Pyramid>(f, mesh_elements.get_element2node());
        break;
      case ElementTopology::element_types::HEX :
        bc_local_face_(f) = ExtractLocalFace
          <ElementTopology::Hexahedron>(f, mesh_elements.get_element2node());
        break;
    }
  }

}// End UnstMeshBcFaces::FormBcFace2Element
//****************************************************************************80
intT UnstMeshBcFaces::FindBcElement(const intT & f, const List2D<intT>& node2element)
{

  intT count = 0; // Counter variable counting

  //---> Counter number of boundary nodes
  intT nbnode = bc_face2node_.get_ncol(f); //bc_face2nodei(f + 1) - bc_face2nodei(f);

  //---> Get first node on boundary face
  intT bnode1 = bc_face2node_(f,0);
  intT e = -1;

  //---> Loop over elements attached to node
  for(intT i = 0; i < node2element.get_ncol(bnode1); i++) {
    // Element on node loop

    //---> Get the ith element attached to node1
    intT elem = node2element(bnode1, i);
    count = 1;

    //---> Loop over the rest of the nodes on bc_face f
    for ( intT j = 1; j < bc_face2node_.get_ncol(f); j ++){
      // bc_node_loop
      intT node = bc_face2node_(f,j);
      //---> Now loop over all the elements of node
      for(intT k = 0; k < node2element.get_ncol(node); k++){
        // node loop
        // Get the elements containing node
        intT elem2 = node2element(node, k);
        /*---> Check to see if the element containing node 2 is the same
          as elem */
        if(elem2 == elem) { //end if elem check
          /*---> Since elem2 is the same as elem...then elem also contains *
            node node.  Therefore elem contains count + 1 nodes of the
            face f */
          count += 1;
          /*---> Leave loop because we don't need to check anymore elements
            containing the node denoted node */
          break;
        } // End if elem check

      } // End node loop

    }// End bc node loop

    if(count == nbnode) { // check count
      e = elem;
      break;
    }// End count check

  } // End element on node loop

    //---> Some error reporting
  if(e == -1) { // End if error
    std::cout << "ERROR: Could not find an element that contained all "
              <<"nodes of the boundary face " << f
              << ".  There ie probably something wrong with the "
              << "connectivity of the mesh.  " << std::endl;
#ifdef DEV_DEBUG
    std::cout << " The face "  << f << " contains the following nodes:"
              << std::endl;
    for(intT i = 0 ; i < bc_face2node_.get_ncol(f); i++){
      std::cout << " Node " << i << ": "
                << bc_face2node_(f,i) << std::endl;
    }
    SystemModule::pause();
#endif

  } // End if error

    //---> Return the element containing face f
  return e;

} //end FindBcElemen

