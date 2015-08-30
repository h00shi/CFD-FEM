#include "UnstMesh.h"

//****************************************************************************80
UnstMesh::UnstMesh(std::string const & filename, std::string const & file_type)
{
  //---> Initialize mesh to empty
  ndim_ = 0;
  nnode_ = 0;
  nelement_ = 0;
  nbc_face_ = 0;
  nbc_id_ = 0;
  ncolor_ = 0;
  nline_ = 0;
  ndof_ = 0;
  nbar_ = 0;
  ntri_ = 0;
  nquad_ = 0;
  ntet_ = 0;
  nprism_ = 0;
  npyr_ = 0;
  nhex_ = 0;
  nbc_node_ = 0;
  nbc_edge_ = 0;
  nbc_tri_ = 0;
  nbc_quad_ = 0;
  grid_mem_ = 0;
  
  if(file_type.compare("Umesh2D-Binary") == 0 ) {
    
  }
  else if(file_type.compare("Grid-NKB") == 0 ){
    ReadGridFile(filename);
  }
  else {
    std::cout << "ERROR: Invalid mesh file type specified." 
	      << "Mesh cannot be constructed." << std::endl
	      << "The supplied mesh file type is " << file_type << std::endl
	      << "Valid file types are Umesh2D-Binary, Grid-NKB" << std::endl;
  
#ifdef DEV_DEBUG
  std::cout << "Error from UnstMesh constructor at line 21 of UnstMesh.cpp" 
	    << std::endl;
#endif
  }

  //---> Now Process the Mesh
  FormConnectivity();
  

}

//****************************************************************************80
UnstMesh::~UnstMesh() {};

//****************************************************************************80
void UnstMesh::Diagnostic(std::ostream& out_stream) 
{
  out_stream << std::endl;
  out_stream << "------------------------- Mesh Diagnostics --------------------"
	    << std::endl;
  out_stream << "Mesh Data: " << std::endl;
  out_stream << "Number of Physical Dimensions: " << ndim_ << std::endl;
  out_stream << "Number of Nodes: " << nnode_ << std::endl;
  out_stream << "Number of Elements: " << nelement_ << std::endl;
  out_stream << "Number of Boundary Faces: " << nbc_face_ << std::endl;
  out_stream << "Number of Boundary IDs: " << nbc_id_ << std::endl;
  out_stream << "Number of Degrees of Freedom: " << ndof_ << std::endl;
  out_stream << "Number of Colors: " << ncolor_ << std::endl;
  out_stream << "Element Data: " << std::endl;
  out_stream << "Number of Bar Elmenets: " << nbar_ << std::endl;
  out_stream << "Number of Tri Elements: " << ntri_ << std::endl;
  out_stream << "Number of Quad Elements: " << nquad_ << std::endl;
  out_stream << "Number of Tet Elements: " << ntet_ << std::endl;
  out_stream << "Number of Prism Elemnets: " << nprism_ << std::endl;
  out_stream << "Number of Pyr Elements: " << npyr_ << std::endl;
  out_stream << "Number of Hex Elements: " << nhex_ << std::endl;
  out_stream << "Boundary Face Data: " << std::endl;
  out_stream << "Number of Boundary Node Faces: " << nbc_node_ << std::endl;
  out_stream << "Number of Boundary Edge Faces: " << nbc_edge_ << std::endl;
  out_stream << "Number of Boundary Tri Faces : " << nbc_tri_ << std::endl;
  out_stream << "Number of Boundary Quad Faces : " << nbc_quad_ << std::endl;
  out_stream << std::endl;
  out_stream << "MemoryDiagnostics: " << std::endl;
 
  out_stream << x_.MemoryDiagnostic("x");
  out_stream << element2node_.MemoryDiagnostic("element2node_"); 
  out_stream << bc_face2elem_.MemoryDiagnostic("bc_face2elem_");
  out_stream << bc_face2node_.MemoryDiagnostic("bc_face2node_");
  out_stream << element_type_.MemoryDiagnostic("element_type_");
  out_stream << bc_face_id_.MemoryDiagnostic("bc_face_id_");
  out_stream << bcid_type_.MemoryDiagnostic("bcid_type_");
  out_stream << bc_local_face_.MemoryDiagnostic("bc_local_face_");
  out_stream << nface_per_bcid_.MemoryDiagnostic("nface_per_bcid_");
  out_stream << node2element_.MemoryDiagnostic("node2element_");
  out_stream << adj_.MemoryDiagnostic("adj_");
  out_stream << element_vol_.MemoryDiagnostic("element_vol_");
  out_stream << element_sa_.MemoryDiagnostic("element_sa_");
  out_stream << edge2node_.MemoryDiagnostic("edge2node_");
  out_stream << std::endl;
  out_stream << "Total Grid Memory: " << grid_mem_ << " MB" << std::endl;
  out_stream << "---------------------- End  Mesh Diagnostics -----------------" 
	    << std::endl << std::endl;
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
  intT nbnode = bc_face2node_.get_ncol(f); //bc_face2nodei(f + 1) - bc_face2nodei(f);

  //---> Get first node on boundary face
  intT bnode1 = bc_face2node_(f,0);
  e = -1;

  //---> Loop over elements attached to node
  for( i = 0; i < node2element_.get_ncol(bnode1); i++) {
    // Element on node loop

    //---> Get the ith element attached to node1
    intT elem = node2element_(bnode1, i);
    count = 1;

    //---> Loop over the rest of the nodes on bc_face f
    for ( j = 1; j < bc_face2node_.get_ncol(f); j ++){
      // bc_node_loop
      intT node = bc_face2node_(f,j);
      //---> Now loop over all the elements of node
      for( k = 0; k < node2element_.get_ncol(node); k++){
        // node loop
        // Get the elements containing node
        intT elem2 = node2element_(node, k);
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
    for( i = 0 ; i < bc_face2node_.get_ncol(f); i++){
      std::cout << " Node " << i << ": "
                << bc_face2node_(f,i) << std::endl;
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
  std::cout << "-------------Node 2 Element Connectivity ---------------------"
	    << std::endl;
  //---> Initialize the temporary counter array
  nodeitemp_.set_value(0);

  /*---> Loop over elements and count the number of elements attached to a
    node */
  for (e = 0; e < nelement_; e++) { // element_loop
    //---> Get number of nodes on this element
    nne = element2node_.get_ncol(e);
    //---> Loop over nodes on the element
    for(n = 0; n < nne; n++) { // node_loop
      //---> Get node number
      node = element2node_(e, n);
      /*---> Add 1 to the count of number of elements attached to a node
        for node */
      nodeitemp_(node + 1) += 1;
    } // End node_loop
  } // End element_loop

    /*---> Now add up all the element surrounding each node and setup
      linked list index array */
  size_node2element = 0;
  for ( n = 0; n < nnode_ ; n++) { // Node_loop
    //---> Size of node2element array computation
    size_node2element += nodeitemp_(n + 1);
  } // End Node_loop

    //---> Allocate the memory for node2element data array
  node2element_.initialize(nnode_, size_node2element);
  grid_mem_ += node2element_.get_mem();

  for ( n = 0; n < nnode_ ; n++) { // Node_loop
    //---> Size of node2element array computation
    node2element_.set_ncol(n, nodeitemp_(n + 1));
    nodeitemp_(n + 1) = 0;
  } // End Node_loop

    //---> Now for the connectivity node2element
  for ( e = 0; e < nelement_; e++ ) { // Element_loop
    //---> Get number of nodes on this element
    nne = element2node_.get_ncol(e);
    //---> Loop over nodes on the element
    for(n = 0; n < nne; n++) { // node_loop
      //---> Get node number
      node = element2node_(e, n);
      //----> Assign element e to index position of node2element
      node2element_(node, nodeitemp_(node)) = e;
      /*----> We've added an element to the node2element for node: node...so
        add 1 to node2elementi(node) to act as a counter for how many
        elements we've added to the node, this acts as a local index now */
      //node2elementi(node) += 1;
      nodeitemp_(node) += 1;
    } // End node_loop
  } // End element_loop

  // DON'T THINK I NEED THIS 
  //---> Sort the the columns for each row so it's easier to deal with later
  // for ( n = 0; n < nnode_; n++) {
  //   intT start = 0;
  //   intT end = node2element.get_ncol(n) - 1;
  //   std::sort(node2element.get_ptr(n,start), node2element.get_ptr(n,end) + 1);
  // }
 
  std::cout << "------------- End Node 2 Element Connectivity -----------------"
	    << std::endl;
 
  return;
}// End Function FormNode2Element

//****************************************************************************80
void UnstMesh::FormBcFaceConnectivity()
{
   std::cout << "-------------- Boundary Face Connectivity --------------------"
	    << std::endl;
  //---> Local Variables
  intT f; // Face looping index
  //---> For each boundary face figure out to which element it belongs
  for (f = 0; f < nbc_face_; f++) {// BC face loop
    intT e = FindBcElement(f);
    bc_face2elem_(f) = e;
  }// End bc_face_loop

  /*---> Now for each boundary face figure out local face number. */
  for (f = 0; f < nbc_face_; f++){// BC face loop

    intT e = bc_face2elem_(f); // Element attached to bc_face f
    intT etype =  element_type_(e); // Element type
    intT loc_face; // Local Face number
    intT nne = element2node_.get_ncol(e); // Number of nodes on the element e
    intT nnf = bc_face2node_.get_ncol(f); // Number of nodes on the face f

    //---> Temp vectors
    Array1D<intT> face_nodes(nnf);
    Array1D<intT> elem_nodes(nne);

    //---> Fill elem_nodes
    for (intT i = 0; i < nne; i++){// elem nodes loop
      //---> Access element connecitivity
      elem_nodes(i) = element2node_(e, i);
    }// End elem nodes loop

    //---> Fill face_nodes
    for (intT i = 0; i < nnf; i++){// face nodes loop
      //---> Acces face connecitivity
      face_nodes(i) = bc_face2node_(f, i);
    }// end face nodes loop

    //---> Based on cell type figure out local face id
    switch (etype) {
    case element_types::ELEMTYPE_BAR: // 1-D Bar elements
      loc_face = GetLocalFace1D(face_nodes, elem_nodes);
      break;
    case element_types::ELEMTYPE_TRI: // 2-D Triangles
      loc_face = GetLocalFaceTri(face_nodes, elem_nodes);
      break;
    case element_types::ELEMTYPE_QUAD: // 2-D Quadrilaterals
      loc_face = GetLocalFaceQuad(face_nodes, elem_nodes);
      break;
    case element_types::ELEMTYPE_TET: // 3-D Tetrahedra
      loc_face = GetLocalFaceTet(face_nodes, elem_nodes);
      break;
    case element_types::ELEMTYPE_PRISM: // 3-D Prism
      loc_face = GetLocalFacePrism(face_nodes, elem_nodes);
      break;
    case element_types::ELEMTYPE_PYR: // 3-D Pyramid
      loc_face = GetLocalFacePyramid(face_nodes, elem_nodes);
      break;
    case element_types::ELEMTYPE_HEX: // 3-D Hex
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
    bc_local_face_(f) = loc_face;
  }// End BC face loop
  
  //---> Count the number of boundary faces for each boundary
  for(intT id = 0; id < nbc_id_; id++){nface_per_bcid_(id) = 0;}
  
  //---> Loop over boundary faces
  for(intT bcf = 0; bcf < nbc_face_; bcf++) { // bc_face loop
    intT id = bc_face_id_(bcf);
    nface_per_bcid_(id) +=1;
  } // End bc_face loop
  
  std::cout << "-------------- End Boundary Face Connectivity -----------------"
	    << std::endl;
}// End FormBcFaceConnectivity

//****************************************************************************80
void UnstMesh::CountElementTypes()
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
    case element_types::ELEMTYPE_BAR:
      nbar_ += 1;
      break;
    case element_types::ELEMTYPE_TRI:
      ntri_ += 1;
      break;
    case element_types::ELEMTYPE_QUAD:
      nquad_ += 1;
      break;
    case element_types::ELEMTYPE_TET:
      ntet_ += 1;
      break;
    case element_types::ELEMTYPE_PRISM:
      nprism_ += 1;
      break;
    case element_types::ELEMTYPE_PYR:
      npyr_ += 1;
      break;
    case element_types::ELEMTYPE_HEX:
      nhex_ += 1;
      break;
    }
  }
}// End CountElementTypes

//****************************************************************************80
void UnstMesh::CountBcFaceTypes()
{
  nbc_node_ = 0;
  nbc_edge_ = 0;
  nbc_tri_ = 0;
  nbc_quad_ = 0;

  for (intT f = 0; f < nbc_face_; f++) { // Count bc_face types
    switch (bc_face2node_.get_ncol(f)) {
    case 1:
      nbc_node_ += 1;
      break;
    case 2:
      nbc_edge_ += 1;
      break;
    case 3:
      nbc_tri_ += 1;
      break;
    case 4:
      nbc_quad_ += 1;
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
  
  //---> Form Adjacency
  FormAdjacency();
  
} // End FormConnectivity

//****************************************************************************80
void UnstMesh::ReadGridFile(std::string const & filename) 
{
  //---> We know we are reading a 2-D mesh file;
  ndim_ = 2;
  
  FILE * mesh_file;
 
  char char_crap[100];
  intT nnz_elem2node = 0;
  intT nnz_bc_face2node = 0;
   
  std::cout << "Opening mesh file: " << filename << std::endl;		\
  /*---> Open the mesh file.  NOTE: must cast string to char* because that 
    is the argument of mesh_file.open(); */
  
  mesh_file = fopen(filename.c_str(), "r");
  
  if( mesh_file == NULL ) { // file_check
     
    std::cout << "ERROR: Could not open mesh file: " << filename 
	      << std::endl;
    SystemModule::my_exit();
      
  } //End file_check


    //---> Get top line: nnodes nelement nbc_face nbc_id
  std::cout << "Reading Header" << std::endl;;
  fscanf(mesh_file, "%d %d %d %d\n", &nnode_, &nelement_ , &nbc_face_,
	 &nbc_id_ );
  
  //---> Memory based on stuff so far
  //----------------------------------------------------------------------------
  nodeitemp_.initialize(nnode_ + 1);
  grid_mem_ += nodeitemp_.get_mem();
  //----------------------------------------------------------------------------
  x_.initialize(nnode_, ndim_);
  grid_mem_ += x_.get_mem();
  //----------------------------------------------------------------------------
  bc_face2elem_.initialize(nbc_face_);
  grid_mem_ += bc_face2elem_.get_mem();
  //----------------------------------------------------------------------------
  element_type_.initialize(nelement_);
  grid_mem_ += element_type_.get_mem();
  //----------------------------------------------------------------------------
  bc_face_id_.initialize(nbc_face_);
  grid_mem_ += bc_face_id_.get_mem();		       
  //----------------------------------------------------------------------------
  bcid_type_.initialize(nbc_id_);
  grid_mem_ += bcid_type_.get_mem();
  //----------------------------------------------------------------------------
  bc_local_face_.initialize(nbc_face_);
  grid_mem_ += bc_local_face_.get_mem();
  //----------------------------------------------------------------------------
  nface_per_bcid_.initialize(nbc_id_);
  grid_mem_ += nface_per_bcid_.get_mem();
  //----------------------------------------------------------------------------
  //---> Get line of text: Nodes per Cell
  // Use fgets to read 100 chars of crap and move one...it's a neat trick
  fgets(char_crap, 100, mesh_file);
  std::cout << "Reading Element Type" << std::endl;
  //---> Loop over elements and get nnodes per element
  for(intT e = 0; e < nelement_; e++){
    intT n = 0; 
    fscanf(mesh_file, "%d\n", &n);
    nnz_elem2node += n;
    
    switch (n) {
    case 2:
      element_type_(e) = element_types::ELEMTYPE_BAR;
      break;
    case 3:
      element_type_(e) = element_types::ELEMTYPE_TRI;
      break;
    case 4:
      element_type_(e) = element_types::ELEMTYPE_QUAD;
      break;
    default:
      std::cout << "ERROR: Element " << e
		<< " has an invalid number of nodes " << n 
		<< ".  Not sure how this happened..."
		<< "but it needs to be fixed 2-D grids only please.  " 
		<< "Bar, Tri and Quad supported " << std::endl;
      
    } // End check number of nodes per element
  }// End loop over elements

  nnz_bc_face2node = 2*nbc_face_;
  //----------------------------------------------------------------------------
  element2node_.initialize(nelement_, nnz_elem2node);
  grid_mem_ += element2node_.get_mem();
  //----------------------------------------------------------------------------
  bc_face2node_.initialize(nbc_face_, nnz_bc_face2node);
  grid_mem_ += bc_face2node_.get_mem();
  //----------------------------------------------------------------------------
  std::cout << "Reading Nodal Coordinates" << std::endl;
  //---> Again using gets read line: Nodal Coordinates
  fgets(char_crap, 100, mesh_file);
  for(intT n = 0; n < nnode_; n++){
    fscanf(mesh_file, "%lf %lf\n", &x_(n,0), &x_(n,1));
  }
  std::cout << "Reading Element 2 Node Connecitivity" << std::endl;
  //---> Again using gets read line: Element to Node Connectivity
  fgets(char_crap, 100, mesh_file);
  //std::cout << char_crap << std::endl;
  for(intT e = 0; e < nelement_; e++){ // Element2node loop 
    intT n0, n1, n2, n3;
   
    switch (element_type_(e)) { // Choose element_type
    case element_types::ELEMTYPE_BAR:
      fscanf(mesh_file, "%d %d\n", &n0, &n1);
      
      element2node_.set_ncol(e,2);
      element2node_(e,0) = n0 - 1;
      element2node_(e,1) = n1 - 1;
      
      break;
    case element_types::ELEMTYPE_TRI:
      fscanf(mesh_file, "%d %d %d\n", &n0, &n1, &n2);
    
      element2node_.set_ncol(e,3);
      element2node_(e,0) = n0 - 1;
      element2node_(e,1) = n1 - 1;
      element2node_(e,2) = n2 - 1;
      
      break;
    case element_types::ELEMTYPE_QUAD:
      fscanf(mesh_file, "%d %d %d %d\n", &n0, &n1, &n2, &n3);
      
      element2node_.set_ncol(e,4);
      element2node_(e,0) = n0 - 1;
      element2node_(e,1) = n1 - 1;
      element2node_(e,2) = n2 - 1;
      element2node_(e,3) = n3 - 1;
      
      break;
    } // End Choose element type
  } // End Element2node loop
  std::cout << "Reading Boundary Face 2 Node Connecitivity" << std::endl;
  //---> Again using gets read line: Boundary Face Data
  fgets(char_crap, 100, mesh_file);
  for(intT f = 0; f < nbc_face_; f++){// Read Bc Faces
    realT sbc = 0.0;
    realT xbc = 0.0;
    realT ybc = 0.0;
    //---> Set 2 nodes per bc_face
    bc_face2node_.set_ncol(f,2);
   
    fscanf(mesh_file, "%d %d %d\n", &bc_face2node_(f,0), 
	   &bc_face2node_(f,1), &bc_face_id_(f));
    
    bc_face2node_(f,0) = bc_face2node_(f,0) - 1;
    bc_face2node_(f,1) = bc_face2node_(f,1) - 1;
    bc_face_id_(f) = bc_face_id_(f) - 1;
    
    for(intT i = 0; i < 11; i++){
      fscanf(mesh_file, "%lf %lf %lf\n", &sbc, &xbc, &ybc);
    }
    
  } // End Read BC faces

  fclose(mesh_file);
  std::cout << "Done Reading mesh file." << std::endl;
  
}// End ReadGridFile

//****************************************************************************80
void UnstMesh::FormAdjacency() 
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
  adj_.initialize(nnode_, nnz_adj);
  grid_mem_ +=adj_.get_mem();
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
    adj_.set_ncol(n,c);
 
    for (intT i = 0; i < c; i++){adj_(n,i) = temp[i];}
    
  } // End Node_loop 


}// End FormAdjacency
