#include "UnstMesh.h"

//****************************************************************************80
UnstMesh::UnstMesh(std::string const & filename, std::string const & file_type)
{
  

  if(file_type.compare("Umesh2D-Binary") == 0 ) {
    
  }
  else if(file_type.compare("Grid-NKB") == 0 ){
    ReadGridFile(filename);
  }
  else
    std::cout << "ERROR: Invalid mesh file type specified." 
	      << "Mesh cannot be constructed." << std::endl
	      << "The supplied mesh file type is " << file_type << std::endl
	      << "Valid file types are Umesh2D-Binary, Grid-NKB" << std::endl;
  
#ifdef DEV_DEBUG
  std::cout << "Error from UnstMesh constructor at line 21 of UnstMesh.cpp" 
	    << std::endl;
#endif
  

}

//****************************************************************************80
void UnstMesh::MemoryDiagnostic() 
{


}

//****************************************************************************80
void UnstMesh::AllocateMemory(const intT& nnz_element2node, 
			 const intT& nnz_bc_face2node) {

  //----------------------------------------------------------------------------
  nodeitemp.initialize(nnode + 1);
  grid_mem += nodeitemp.get_mem();
  //----------------------------------------------------------------------------
  element2node.initialize(nnode, nnz_element2node);
  grid_mem += element2node.get_mem();
  //----------------------------------------------------------------------------
  bc_face2elem.initialize(nbc_face);
  grid_mem += bc_face2elem.get_mem();
  //----------------------------------------------------------------------------
  bc_face2node.initialize(nbc_face, nnz_bc_face2node);
  grid_mem += bc_face2node.get_mem();
  //----------------------------------------------------------------------------
  element_type.initialize(nelement);
  grid_mem += element_type.get_mem();
  //----------------------------------------------------------------------------
  bc_face_id.initialize(nbc_face);
  grid_mem += bc_face_id.get_mem();		       
  //----------------------------------------------------------------------------
  bcid_type.initialize(nbc_id);
  grid_mem += bcid_type.get_mem();
  //----------------------------------------------------------------------------
  bc_local_face.initialize(nbc_face);
  grid_mem += bc_local_face.get_mem();
  //----------------------------------------------------------------------------
  nface_per_bcid.initialize(nbc_id);
  grid_mem += nface_per_bcid.get_mem();
  //----------------------------------------------------------------------------
}

//****************************************************************************80
intT UnstMesh::FindBcElement(const intT & f)
{
  //---> Local Variables
  intT i; // Looping index
  intT j; // Looping index
  intT k; // Looping index
  intT count = 0; // Counter variable counting
  intT e; // Return value

  //---> Counter number of boundary nodes
  intT nbnode = bc_face2node.get_ncol(f); //bc_face2nodei(f + 1) - bc_face2nodei(f);

  //---> Get first node on boundary face
  intT bnode1 = bc_face2node(f,0);
  e = -1;

  //---> Loop over elements attached to node
  for( i = 0; i < node2element.get_ncol(bnode1); i++) {
    // Element on node loop

    //---> Get the ith element attached to node1
    intT elem = node2element(bnode1, i);
    count = 1;

    //---> Loop over the rest of the nodes on bc_face f
    for ( j = 1; j < bc_face2node.get_ncol(f); j ++){
      // bc_node_loop
      intT node = bc_face2node(f,j);
      //---> Now loop over all the elements of node
      for( k = 0; k < node2element.get_ncol(node); k++){
        // node loop
        // Get the elements containing node
        intT elem2 = node2element(node, k);
        /*---> Check to see if the element containing node 2 is the same
          as elem */
        if(elem2 == elem) { //end if elem check
          /*---> Since elem2 is the same as elem...then elem also contains *
            node node.  Thefefore elem contains count + 1 nodes of the
            face f */
          count += 1;
          /*---> Leave loop because we don't need to check anymore elements
            containg the node denoted node */
          break;
        } // End if elem check

      } // End node loop

    }// End bc node loop

    if( count == nbnode) { // check count
      e = elem;
      break;
    }// End count check

  } // End element on node loop

    //---> Some error reporting
  if( e == -1 ) { // End if error
    std::cout << "ERROR: Could not find an element that contained all "
              <<"nodes of the boundary face " << f
              << ".  There ie probably something wrong with the "
              << "connectivity of the mesh.  " << std::endl;
#ifdef DEV_DEBUG
    std::cout << " The face "  << f << " contains the following nodes:"
              << std::endl;
    for( i = 0 ; i < bc_face2node.get_ncol(f); i++){
      std::cout << " Node " << i << ": "
                << bc_face2node(f,i) << std::endl;
    }
    SystemModule::pause();
#endif

  } // End if error

    //---> Return the element containing face f
  return(e);

} //end FindBcElement

//****************************************************************************80
intT UnstMesh::GetLocalFace1D(const Array1D<intT>& face_nodes,
			      const Array1D<intT>& elem_nodes)
{
  //---> Return variable
  intT loc_face = -99; //---> Initialized to - 1 to throw error if it's wrong

  //---> Since it's a 1-D element it's very easy, face nodes is only 1 node
  if( face_nodes(0) == elem_nodes(0) ) loc_face = 0;
  if( face_nodes(0) == elem_nodes(1) ) loc_face = 1;

  return(loc_face);

} // End GetLocalFace1D

//****************************************************************************80
intT UnstMesh::GetLocalFaceTri(const Array1D<intT>& face_nodes,
			       const Array1D<intT>& elem_nodes)
{
  //---> Return variable
  intT loc_face = -99; //---> Initialized to -99 to throw error if it's wrong
  Array1D<intT> loc_nodes(2);

  for(intT i = 0; i < elem_nodes.get_size(0); i++) { // elem node loop
    for(intT j = 0; j < face_nodes.get_size(0); j++) {
      if( face_nodes(j) == elem_nodes(i) ) {
        loc_nodes(j) = i;
      }
    }
  } // End elem node loop
    //---> Local Node numbers sume to unique value
  intT sum = 0;
  for (intT i = 0; i < 2; i++){
    sum += loc_nodes(i);
  }

  if( sum == 1) loc_face = 0;
  if( sum == 3) loc_face = 1;
  if( sum == 2) loc_face = 2;
  return(loc_face);

} // End GetLocalFaceTri

//****************************************************************************80
intT UnstMesh::GetLocalFaceQuad(const Array1D<intT>& face_nodes,
				const Array1D<intT>& elem_nodes)
{
  //---> Return variable
  intT loc_face = -99; //---> Initialized to - 99 to throw error if it's wrong
  Array1D<intT> loc_nodes(2);

  //---> Local Node numbers + 1 multiply to a unique value
  for(intT i = 0; i < elem_nodes.get_size(0); i++) { // elem node loop
    for(intT j = 0; j < face_nodes.get_size(0); j++) {
      if( face_nodes(j) == elem_nodes(i) ) {
        loc_nodes(j) = i;
      }
    }
  } // End elem node loop

  intT prod = 1;
  for (intT i = 0; i < 2; i++){
    prod *= (loc_nodes(i) + 1);
  }
  if(prod == 2 ) loc_face = 0;
  if(prod == 6 ) loc_face = 1;
  if(prod == 12) loc_face = 2;
  if(prod == 4 ) loc_face = 3;

  return(loc_face);

} // End GetLocalFaceQuad

//****************************************************************************80
intT UnstMesh::GetLocalFaceTet(const Array1D<intT>& face_nodes,
                                            const Array1D<intT>& elem_nodes)
{
  //---> Return variable
  intT loc_face = -99; //---> Initialized to - 99 to throw error if it's wrong
  Array1D<intT> loc_nodes(3);
  //---> Now we are in 3-D...not so easy anymore.

  for(intT i = 0; i < elem_nodes.get_size(0); i++) { // elem node loop
    for(intT j = 0; j < face_nodes.get_size(0); j++) {
      if( face_nodes(j) == elem_nodes(i) ) {
        loc_nodes(j) = i;
      }
    }
  } // End elem node loop

  intT sum = 0;
  for (intT i = 0; i < 3; i++){
    sum += loc_nodes(i);
  }
  /*---> For a tetrahedra the sum of the local node numbers for each face is
    unique.  Therefore we use this property to find the local faces numbers
    for a tet.  If the element is not a tet we're going to need to find a
    cuter trick for those.  */
  if( sum == 4 ) loc_face = 0;
  if( sum == 6 ) loc_face = 1;
  if( sum == 5 ) loc_face = 2;
  if( sum == 3 ) loc_face = 3;

  return(loc_face);
} // End GetLocalFaceTet

//****************************************************************************80
intT UnstMesh::GetLocalFaceHex(const Array1D<intT>& face_nodes,
                                            const Array1D<intT>& elem_nodes)
{
  //---> Return variable
  intT loc_face = -99; //---> Initialized to - 99 to throw error if it's wrong
  Array1D<intT> loc_nodes(4);
  //---> Now we are in 3-D...not so easy anymore.

  for(intT i = 0; i < elem_nodes.get_size(0); i++) { // elem node loop
    for(intT j = 0; j < face_nodes.get_size(0); j++) {
      if( face_nodes(j) == elem_nodes(i) ) {
        loc_nodes(j) = i;
      }
    }
  } // End elem node loop

  intT prod = 1;
  for (intT i = 0; i < 4; i++){
    prod *= (loc_nodes(i) + 1);
  }
  /*---> For a hexahedra the product of the local node numbers for each
    face is unique.  Therefore we use this property to find the local
    faces numbers for a hex. */
  if( prod == 160  ) loc_face = 0;
  if( prod == 252  ) loc_face = 1;
  if( prod == 60   ) loc_face = 2;
  if( prod == 672  ) loc_face = 3;
  if( prod == 24   ) loc_face = 4;
  if( prod == 1680 ) loc_face = 5;
  return(loc_face);
} // End GetLocalFaceHex

//****************************************************************************80
intT UnstMesh::GetLocalFacePrism(const Array1D<intT>& face_nodes,
				  const Array1D<intT>& elem_nodes)
{
  //---> Return variable
  intT loc_face = -99; //---> Initialized to - 99 to throw error if it's wrong
  Array1D<intT> loc_nodes(4);
  //---> Now we are in 3-D...not so easy anymore.

  for(intT i = 0; i < elem_nodes.get_size(0); i++) { // elem node loop
    for(intT j = 0; j < face_nodes.get_size(0); j++) {
      if( face_nodes(j) == elem_nodes(i) ) {
        loc_nodes(j) = i;
      }
    }
  } // End elem node loop

  intT prod = 1;
  for (intT i = 0; i < face_nodes.get_size(0); i++){
    prod *= (loc_nodes(i) + 1);
  }
  /*---> For a prism the product of the local node numbers for each
    face is unique.  Therefore we use this property to find the local
    faces numbers for a hex. */
  if( prod == 40  ) loc_face = 0;
  if( prod == 180 ) loc_face = 1;
  if( prod == 72  ) loc_face = 2;
  if( prod == 6   ) loc_face = 3;
  if( prod == 120 ) loc_face = 4;
  return(loc_face);
} // End GetLocalFacePrism

//****************************************************************************80
intT UnstMesh::GetLocalFacePyramid(const Array1D<intT>& face_nodes,
				    const Array1D<intT>& elem_nodes)
{
  //---> Return variable
  intT loc_face = -99; //---> Initialized to - 99 to throw error if it's wrong
  Array1D<intT> loc_nodes(4);
  //---> Now we are in 3-D...not so easy anymore.

  for(intT i = 0; i < elem_nodes.get_size(0); i++) { // elem node loop
    for(intT j = 0; j < face_nodes.get_size(0); j++) {
      if( face_nodes(j) == elem_nodes(i) ) {
        loc_nodes(j) = i;
      }
    }
  } // End elem node loop

  intT prod = 1;
  for (intT i = 0; i < face_nodes.get_size(0); i++){
    prod *= (loc_nodes(i) + 1);
  }
  /*---> For a pyramid the product of the local node numbers for each
    face is unique.  Therefore we use this property to find the local
    faces numbers for a hex. */
  if( prod == 24 ) loc_face = 0;
  if( prod == 10 ) loc_face = 1;
  if( prod == 30 ) loc_face = 2;
  if( prod == 60 ) loc_face = 3;
  if( prod == 20 ) loc_face = 4;
  return(loc_face);
} // End GetLocalFacePyramid


//****************************************************************************80
void UnstMesh::FormNode2Element()
{
  //---> Local Variables
  intT e; // The element looping index
  intT n; // The nodes on element looping index
  intT node; // The node number
  intT nne; // Number of nodes attached to element
  intT size_node2element; // Size of the array node2element

  //---> Initialize the temporary counter array
  nodeitemp.set_value(0);

  /*---> Loop over elements and count the number of elements attached to a
    node */
  for (e = 0; e < nelement; e++) { // element_loop
    //---> Get number of nodes on this element
    nne = element2node.get_ncol(e);
    //---> Loop over nodes on the element
    for(n = 0; n< nne; n++) { // node_loop
      //---> Get node number
      node = element2node(e, n);
      /*---> Add 1 to the count of number of elements attached to a node
        for node */
      nodeitemp(node + 1) += 1;
    } // End node_loop
  } // End element_loop

    /*---> Now add up all the element surrounding each node and setup
      linked list index array */
  size_node2element = 0;
  for ( n = 0; n < nnode ; n++) { // Node_loop
    //---> Size of node2element array computation
    size_node2element += nodeitemp(n + 1);
  } // End Node_loop

    //---> Allocate the memory for node2element data array
  node2element.initialize(nnode, size_node2element);
  grid_mem += node2element.get_mem();

  for ( n = 0; n < nnode ; n++) { // Node_loop
    //---> Size of node2element array computation
    node2element.set_ncol(n, nodeitemp(n + 1));
    nodeitemp(n + 1) = 0;
  } // End Node_loop

    //---> Now for the connectivity node2element
  for ( e = 0; e < nelement; e++ ) { // Element_loop
    //---> Get number of nodes on this element
    nne = element2node.get_ncol(e);
    //---> Loop over nodes on the element
    for(n = 0; n < nne; n++) { // node_loop
      //---> Get node number
      node = element2node(e, n);
      //----> Assign element e to index position of node2element
      node2element(node, nodeitemp(node)) = e;
      /*----> We've added an element to the node2element for node: node...so
        add 1 to node2elementi(node) to act as a counter for how many
        elements we've added to the node, this acts as a local index now */
      //node2elementi(node) += 1;
      nodeitemp(node) += 1;
    } // End node_loop
  } // End element_loop
  return;
}// End Function FormNode2Element

//****************************************************************************80
void UnstMesh::FormBcFaceConnectivity()
{
  //Initialize
  bc_face2elem.initialize(nbc_face);
  grid_mem += bc_face2elem.get_mem();
  bc_local_face.initialize(nbc_face);
  grid_mem += bc_local_face.get_mem();

  //---> Local Variables
  intT f; // Face looping index
  //---> For each boundary face figure out to which element it belongs
  for (f = 0; f < nbc_face; f++) {// BC face loop
    intT e = FindBcElement(f);
    bc_face2elem(f) = e;
  }// End bc_face_loop

  /*---> Now for each boundary face figure out local face number. */
  for (f = 0; f< nbc_face; f++){// BC face loop

    intT e = bc_face2elem(f); // Element attached to bc_face f
    intT etype =  element_type(e); // Element type
    intT loc_face; // Local Face number
    intT nne = element2node.get_ncol(e); // Number of nodes on the element e
    intT nnf = bc_face2node.get_ncol(f); // Number of nodes on the face f

    //---> Temp vectors
    Array1D<intT> face_nodes(nnf);
    Array1D<intT> elem_nodes(nne);

    //---> Fill elem_nodes
    for (intT i = 0; i < nne; i++){// elem nodes loop
      //---> Access element connecitivity
      elem_nodes(i) = element2node(e, i);
    }// End elem nodes loop

    //---> Fill face_nodes
    for (intT i = 0; i < nnf; i++){// face nodes loop
      //---> Acces face connecitivity
      face_nodes(i) = bc_face2node(f, i);
    }// end face nodes loop

    //---> Based on cell type figure out local face id
    switch (etype) {
    case element_types::ELEMTYPE_BAR: // 1-D Bar elements
      loc_face = GetLocalFace1D(face_nodes, elem_nodes);
      break;
    case element_types::ELEMTYPE_TRIANGLE: // 2-D Triangles
      loc_face = GetLocalFaceTri(face_nodes, elem_nodes);
      break;
    case element_types::ELEMTYPE_QUAD: // 2-D Quadrilaterals
      loc_face = GetLocalFaceQuad(face_nodes, elem_nodes);
      break;
    case element_types::ELEMTYPE_TETRAHEDRA: // 3-D Tetrahedra
      loc_face = GetLocalFaceTet(face_nodes, elem_nodes);
      break;
    case element_types::ELEMTYPE_PRISM: // 3-D Prism
      loc_face = GetLocalFacePrism(face_nodes, elem_nodes);
      break;
    case element_types::ELEMTYPE_PYRAMID: // 3-D Pyramid
      loc_face = GetLocalFacePyramid(face_nodes, elem_nodes);
      break;
    case element_types::ELEMTYPE_HEXAHEDRA: // 3-D Hex
      loc_face = GetLocalFaceHex(face_nodes, elem_nodes);
      break;
    }

    if( loc_face == -99 ) std::cout << "ERROR: Could not find local "
                                    << "bc_face index on bc-face "
                                    << f
                                    << " For element " << e
                                    << " which is of type: " << etype
                                    << std::endl;
    //---> Assign the value of bc_local_face for bc_face f
    bc_local_face(f) = loc_face;
  }// End BC face loop
  
  //---> Count the number of boundary faces for each boundary
  for(intT id = 0; id < nbc_id; id++){nface_per_bcid(id) = 0;}
  
  //---> Loop over boundary faces
  for(intT bcf = 0; bcf < nbc_face; bcf++) { // bc_face loop
    intT id = bc_face_id(bcf);
    nface_per_bcid(id) +=1;
  } // End bc_face loop
  
}// End FormBcFaceConnectivity

//****************************************************************************80
void UnstMesh::CountElementTypes()
{

  nBar = 0;
  nTri = 0;
  nQuad = 0;
  nTet = 0;
  nPrism = 0;
  nPyr = 0;
  nHex = 0;

  for (intT e = 0; e < nelement; e++) { // Count element types
    switch (element_type(e)) {
    case element_types::ELEMTYPE_BAR:
      nBar += 1;
      break;
    case element_types::ELEMTYPE_TRIANGLE:
      nTri += 1;
      break;
    case element_types::ELEMTYPE_QUAD:
      nQuad += 1;
      break;
    case element_types::ELEMTYPE_TETRAHEDRA:
      nTet += 1;
      break;
    case element_types::ELEMTYPE_PRISM:
      nPrism += 1;
      break;
    case element_types::ELEMTYPE_PYRAMID:
      nPyr += 1;
      break;
    case element_types::ELEMTYPE_HEXAHEDRA:
      nHex += 1;
      break;
    }
  }
}// End CountElementTypes

//****************************************************************************80
void UnstMesh::CountBcFaceTypes()
{
  nBcNode = 0;
  nBcEdge = 0;
  nBcTri = 0;
  nBcQuad = 0;

  for (intT f = 0; f < nbc_face; f++) { // Count bc_face types
    switch (bc_face2node.get_ncol(f)) {
    case 1:
      nBcNode += 1;
      break;
    case 2:
      nBcEdge += 1;
      break;
    case 3:
      nBcTri += 1;
      break;
    case 4:
      nBcQuad += 1;
      break;
    }
  }
}// End CoundBcFaceTypes

//****************************************************************************80
void UnstMesh::FormConnectivity()
{
  //---> Get the elements surrounding a node
  FormNode2Element();

  //---> Get boundary face connectivity
  FormBcFaceConnectivity();

  //---> Count element types and store connectivity for each type...for
  //     CGNS Output;
  CountElementTypes();

  //---> Count bc face types and store connectivity for each type...for
  //     CGNS Output;
  CountBcFaceTypes();
} // End FormConnectivity

//****************************************************************************80
void UnstMesh::ReadGridFile(std::string const & filename) 
{


}
