#include "IO/UnstMeshReaderNKBGrid.h"
//****************************************************************************80
UnstMeshReaderNKBGrid::UnstMeshReaderNKBGrid(){is_open_ = false;}
//****************************************************************************80
UnstMeshReaderNKBGrid::~UnstMeshReaderNKBGrid(){is_open_ = false;}
//****************************************************************************80
void UnstMeshReaderNKBGrid::Open(const std::string& filename)
{
   //---> We know we are reading a 2-D mesh file;
  intT ndim = 2;
  
  char char_crap[100];
  intT nnode;
  intT nelement;
  intT nbc_face;
  intT nbc_id;
  intT nnz_elem2node = 0;
  intT nnz_bc_face2node = 0;
   
  SystemModule::cout << "Opening mesh file: " << filename << std::endl;
  /*---> Open the mesh file.  NOTE: must cast string to char* because that 
    is the argument of mesh_file_.open(); */
  
  mesh_file_ = fopen(filename.c_str(), "r");
  
  if( mesh_file_ == NULL ) { // file_check
     
    std::cout << "ERROR: Could not open mesh file: " << filename 
	      << std::endl;
    SystemModule::my_exit();
      
  } //End file_check


    //---> Get top line: nnodes nelement nbc_face nbc_id
  std::cout << "Reading Header" << std::endl;
  fscanf(mesh_file_, "%d %d %d %d\n", &nnode, &nelement , &nbc_face, &nbc_id);
  
  x_.initialize(nnode, ndim);
  element_type_.initialize(nelement);
  element_region_.initialize(nelement);
  bc_face_id_.initialize(nbc_face);
       
   //---> Get line of text: Nodes per Cell
  // Use fgets to read 100 chars of crap and move one...it's a neat trick
  fgets(char_crap, 100, mesh_file_);
  std::cout << "Reading Element Type" << std::endl;
  //---> Loop over elements and get nnodes per element
  for(intT e = 0; e < nelement; e++){
    intT n = 0; 
    fscanf(mesh_file_, "%d\n", &n);
    nnz_elem2node += n;
    
    switch (n) {
    case 2:
      element_type_(e) = ElementTopology::element_types::BAR;
      break;
    case 3:
      element_type_(e) = ElementTopology::element_types::TRI;
      break;
    case 4:
      element_type_(e) = ElementTopology::element_types::QUAD;
      break;
    default:
      std::cout << "ERROR: Element " << e
		<< " has an invalid number of nodes " << n 
		<< ".  Not sure how this happened..."
		<< "but it needs to be fixed 2-D grids only please.  " 
		<< "Bar, Tri and Quad supported " << std::endl;
      
    } // End check number of nodes per element
  }// End loop over elements

  nnz_bc_face2node = 2*nbc_face;
  element2node_.initialize(nelement, nnz_elem2node);
  bc_face2node_.initialize(nbc_face, nnz_bc_face2node);
 
  //----------------------------------------------------------------------------
  std::cout << "Reading Nodal Coordinates" << std::endl;
  //---> Again using gets read line: Nodal Coordinates
  fgets(char_crap, 100, mesh_file_);
  for(intT n = 0; n < nnode; n++){
    fscanf(mesh_file_, "%lf %lf\n", &x_(n,0), &x_(n,1));
  }
  std::cout << "Reading Element 2 Node Connecitivity" << std::endl;
  //---> Again using gets read line: Element to Node Connectivity
  fgets(char_crap, 100, mesh_file_);
  //std::cout << char_crap << std::endl;
  for(intT e = 0; e < nelement; e++){ // Element2node loop 
    intT n0, n1, n2, n3;
   
    switch (element_type_(e)) { // Choose element_type
    case ElementTopology::element_types::BAR:
      fscanf(mesh_file_, "%d %d\n", &n0, &n1);
      
      element2node_.set_ncol(e,2);
      element2node_(e,0) = n0 - 1;
      element2node_(e,1) = n1 - 1;
      
      break;
    case ElementTopology::element_types::TRI:
      fscanf(mesh_file_, "%d %d %d\n", &n0, &n1, &n2);
    
      element2node_.set_ncol(e,3);
      element2node_(e,0) = n0 - 1;
      element2node_(e,1) = n1 - 1;
      element2node_(e,2) = n2 - 1;
      
      break;
    case ElementTopology::element_types::QUAD:
      fscanf(mesh_file_, "%d %d %d %d\n", &n0, &n1, &n2, &n3);
      
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
  fgets(char_crap, 100, mesh_file_);
  for(intT f = 0; f < nbc_face; f++){// Read Bc Faces
    realT sbc = 0.0;
    realT xbc = 0.0;
    realT ybc = 0.0;
    //---> Set 2 nodes per bc_face
    bc_face2node_.set_ncol(f,2);
   
    fscanf(mesh_file_, "%d %d %d\n", &bc_face2node_(f,0), 
	   &bc_face2node_(f,1), &bc_face_id_(f));
    
    bc_face2node_(f,0) = bc_face2node_(f,0) - 1;
    bc_face2node_(f,1) = bc_face2node_(f,1) - 1;
    bc_face_id_(f) = bc_face_id_(f) - 1;
    
    for(intT i = 0; i < 11; i++){
      fscanf(mesh_file_, "%lf %lf %lf\n", &sbc, &xbc, &ybc);
    }
    
  } // End Read BC faces

  fclose(mesh_file_);
  std::cout << "Done Reading mesh file." << std::endl;


}// End UnstMeshReaderNKBGrid::Open
//****************************************************************************80
Array2D<realT> UnstMeshReaderNKBGrid::ReadNodes(){return x_;}
