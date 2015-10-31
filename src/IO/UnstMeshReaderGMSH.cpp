/*
 * UnstMeshReaderGMSH.cpp
 *
 *  Created on: Oct 11, 2015
 *      Author: rabbit
 */
#include "IO/UnstMeshReaderGMSH.h"
//****************************************************************************80
UnstMeshReaderGMSH::
UnstMeshReaderGMSH(const std::string& filename,
                   const std::string& imap_filename) :
                   UnstMeshReader(filename)
{
  //---> Read IDMAP
  ReadIdMap(imap_filename);

  //---> Open and insect GMSH header
  SystemModule::cout << "Opening mesh file: " << UnstMeshReader::filename_
      << std::endl;
  /*---> Open the mesh file.  NOTE: must cast string to char* because that
        is the argument of mesh_file_.open(); */

  mesh_file_ = fopen(UnstMeshReader::filename_.c_str(), "r");

  if( mesh_file_ == NULL ) { // file_check

    SystemModule::cout<< "ERROR: Could not open mesh file: "
        << UnstMeshReader::filename_
        << std::endl;
    SystemModule::my_exit();

  } //End file_check;

  //---> Read Header;
  file_type ft = ReadHeader();

  //---> Now read file based on type
  switch (ft){
    case file_type::ASCII :
      ReadFileASCII();
      break;
    case file_type::BINARY :
      ReadFileBinary();
      break;
  }

  fclose(mesh_file_);
}// End UnstMeshReaderGMSH::UnstMeshReaderGMSH
//****************************************************************************80
void UnstMeshReaderGMSH::ReadFileASCII()
{
  //---> Two lines of Crap
  SkipLine();
  SkipLine();
  ReadNodesASCII(UnstMeshReaderGMSH::read_mode::STORE);
  SkipLine();
  SkipLine();
  CountElementsASCII();

  //---> Rewind file;
  rewind(mesh_file_);
  ReadHeader();
  SkipLine();
  SkipLine();
  ReadNodesASCII(read_mode::SKIP);
  //---> Two more lines of crap
  SkipLine();
  SkipLine();
  //---> Read Elements
  ReadElementsASCII();

}// End ReadASCII

//****************************************************************************80
void UnstMeshReaderGMSH::ReadFileBinary()
{
  //---> Two lines of Crap
  SkipLine();
  SkipLine();
  ReadNodesBinary(read_mode::STORE);
  SkipLine();
  SkipLine();

  //---> Read Elements
  CountElementsBinary();

  //---> Rewind file to read in element connecitivity
  rewind(mesh_file_);
  ReadHeader();
  SkipLine();
  SkipLine();
  ReadNodesBinary(read_mode::SKIP);
  //---> Two more lines of crap
  SkipLine();
  SkipLine();
  ReadElementsBinary();

  return;
}

//****************************************************************************80
void UnstMeshReaderGMSH::ReadIdMap(const std::string& idmap_filename)
{
  SystemModule::cout << "Opening IDMAP file: " << idmap_filename
      << std::endl;
  /*---> Open the mesh file.  NOTE: must cast string to char* because that
         is the argument of mesh_file_.open(); */

  idmap_file_ = fopen(idmap_filename.c_str(), "r");

  if( idmap_file_ == NULL ) { // file_check

    SystemModule::cout<< "ERROR: Could not open IDMAP file: "
        << idmap_filename
        << std::endl;
    SystemModule::my_exit();

  } //End file_check
  char char_crap[100];
  SystemModule::cout << "Reading IDMAP File... " << std::endl;
  //---> Read Header
  fgets(char_crap, 100, idmap_file_);

  //---> Read number of GMSH ID's
  intT ngmsh_id;
  fscanf(idmap_file_, "%d \n", &ngmsh_id);
  fgets(char_crap, 100, idmap_file_);
  gmsh_id_.initialize(ngmsh_id);

  //---> Read Boundary Boundary ID Mapping
  fscanf(idmap_file_, "%d \n", &(UnstMeshReader::nbc_id_));

  fgets(char_crap, 100, idmap_file_);
  for(intT i = 0; i < UnstMeshReader::nbc_id_; i++){
    intT gmsh_id;
    intT bc_id;
    fscanf(idmap_file_,"%d %d \n", &gmsh_id, &bc_id);
    gmsh_id_(i) = gmsh_id;
    gmsh_id_2_bc_id_[gmsh_id] = bc_id;
    gmsh_id_is_bc_[gmsh_id] = true;
  }

  //---> Read Region ID Mapping
  fgets(char_crap, 100, idmap_file_);
  fscanf(idmap_file_, "%d \n", &n_region_id_);
  //region_id_2_global_id_.initialize(n_region_id_);
  fgets(char_crap, 100, idmap_file_);
  for(intT i = 0; i < n_region_id_; i++){
    intT gmsh_id;
    intT region_id;
    fscanf(idmap_file_,"%d %d \n", &gmsh_id, &region_id);
    gmsh_id_(UnstMeshReader::nbc_id_+ i) = gmsh_id;
    // region_id_2_global_id_(n2) = n1;
    gmsh_id_2_region_[gmsh_id] = region_id;
    gmsh_id_is_region_[gmsh_id] = true;
  }

  SystemModule::cout << "Boundary Mapping information " << std::endl;
  for(intT i = 0; i < UnstMeshReader::nbc_id_; i++) {
    std::cout << "Boundary : "<< i << " GMSH ID: " << gmsh_id_(i)
                            << " Boundary ID: " << gmsh_id_2_bc_id_[gmsh_id_(i)];
    if( gmsh_id_is_bc_[gmsh_id_(i)] ){
      std::cout << " This is a boundary!";
    }
    std::cout << std::endl;
  }
  for(intT i = 0; i < n_region_id_; i++) {
    intT j = UnstMeshReader::nbc_id_ + i;
    std::cout << "Region : "<< i << " GMSH ID: "
        << gmsh_id_(j)
        << " Region ID: " << gmsh_id_2_region_[gmsh_id_(j)];
    if(gmsh_id_is_region_[gmsh_id_(j)]){
      std::cout << " This is a region!";
    }
    std::cout << std::endl;
  }

  SystemModule::cout << "Done Reading IDMAP File." << std::endl;
  fclose(idmap_file_);
}// End UnstMeshReaderGMSH::ReadIdMap

//****************************************************************************80
Array2D<realT> UnstMeshReaderGMSH::ReadNodes()
{

  Array2D<realT> xtmp(x_.get_size(0), x_.get_size(1));

  memcpy(xtmp.get_ptr(0,0),
         x_.get_ptr(0,0),
         sizeof(realT)*x_.get_total_size());

  return xtmp;
}// End UnstMeshReaderGMSH::ReadNodes()

//****************************************************************************80
List2D<intT> UnstMeshReaderGMSH::ReadElement2Node()
{


  Array1D<intT> ncol_e2n(UnstMeshReader::nelement_);
  for(intT e = 0; e < UnstMeshReader::nelement_; e++){
    intT entity = elem2entity_(e);
    ncol_e2n(e) = entity2node_.get_ncol(entity);
  }
  List2D<intT> e2n_tmp(ncol_e2n);
  for(intT e = 0; e < UnstMeshReader::nelement_; e++){
    intT entity = elem2entity_(e);
    for(intT i = 0; i < e2n_tmp.get_ncol(e); i++) {
      e2n_tmp(e,i) = entity2node_(entity,i);
    }
  }

  return e2n_tmp;
}// End UnstMeshReaderGMSH::ReadElement2Node

//****************************************************************************80
Array1D<ElementTopology::element_types> UnstMeshReaderGMSH::ReadElementType()
{

  Array1D<ElementTopology::element_types> etype_tmp(UnstMeshReader::nelement_);

  for(intT e = 0; e < UnstMeshReader::nelement_; e++){
    intT entity = elem2entity_(e);
    switch(gmsh_type_map_[entity_type_(entity)]) {
      case ElementTopology::element_types::BAR :
        etype_tmp(e) = ElementTopology::element_types::BAR;
        break;
      case ElementTopology::element_types::TRI :
        etype_tmp(e) = ElementTopology::element_types::TRI;
        break;
      case ElementTopology::element_types::QUAD :
        etype_tmp(e) = ElementTopology::element_types::QUAD;
        break;
      case ElementTopology::element_types::TET :
        etype_tmp(e) = ElementTopology::element_types::TET;
        break;
      case ElementTopology::element_types::PRISM :
        etype_tmp(e) = ElementTopology::element_types::PRISM;
        break;
      case ElementTopology::element_types::PYR :
        etype_tmp(e) = ElementTopology::element_types::PYR;
        break;
      case ElementTopology::element_types::HEX :
        etype_tmp(e) = ElementTopology::element_types::HEX;
        break;
    }

  }

  return etype_tmp;

}// End UnstMeshReaderGMSH::ReadElementType
//****************************************************************************80
Array1D<intT> UnstMeshReaderGMSH::ReadElementRegion()
{

  Array1D<intT> ereg_tmp(UnstMeshReader::nelement_);

  for(intT e = 0; e < UnstMeshReader::nelement_; e++){
    intT entity = elem2entity_(e);
    ereg_tmp(e) = gmsh_id_2_region_[entity_type_(entity)];
  }

  return ereg_tmp;

}// End UnstMeshReaderGMSH::ReadElementRegion()
//****************************************************************************80
List2D<intT> UnstMeshReaderGMSH::ReadBcFace2Node()
{

  Array1D<intT> ncol(UnstMeshReader::nbc_face_);
  for(intT f = 0; f < UnstMeshReader::nbc_face_; f++){
    intT entity = bc_face2entity_(f);
    ncol(f) = entity2node_.get_ncol(entity);
  }
  List2D<intT> bcf2n_tmp(ncol);
  for(intT f = 0; f < UnstMeshReader::nbc_face_; f++){
    intT entity = bc_face2entity_(f);
    for(intT i = 0; i < bcf2n_tmp.get_ncol(f); i++){
      bcf2n_tmp(f,i) = entity2node_(entity,i);
    }
  }

  return bcf2n_tmp;

}// End UnstMeshReaderGMSH::ReadBcFace2Node
//****************************************************************************80
Array1D<intT> UnstMeshReaderGMSH::ReadBcID()
{

  Array1D<intT> bcid_tmp(UnstMeshReader::nbc_face_);
  for(intT f = 0; f < UnstMeshReader::nbc_face_; f++){
    intT entity = bc_face2entity_(f);
    bcid_tmp(f) = gmsh_id_2_bc_id_[entity_id_(entity)];
  }

  return bcid_tmp;

}// End UnstMeshReaderGMSH::ReadBcID
//****************************************************************************80
Array1D<ElementTopology::face_types> UnstMeshReaderGMSH::ReadBcFaceType()
{

  Array1D<ElementTopology::face_types> bcftype_tmp(UnstMeshReader::nbc_face_);

  for(intT f = 0; f < UnstMeshReader::nbc_face_; f++){
    intT entity = bc_face2entity_(f);
    switch(gmsh_type_map_[entity_type_(entity)]) {
      case ElementTopology::element_types::BAR :
        bcftype_tmp(f) = ElementTopology::face_types::FACE_BAR;
        break;
      case ElementTopology::element_types::TRI :
        bcftype_tmp(f) = ElementTopology::face_types::FACE_TRI;
        break;
      case ElementTopology::element_types::QUAD :
        bcftype_tmp(f) = ElementTopology::face_types::FACE_QUAD;
        break;

    }
  }

  return bcftype_tmp;

}// End UnstMeshReaderGMSH::ReadBcFaceType
//****************************************************************************80
void UnstMeshReaderGMSH::SkipLine()
{
  char char_crap[100];
  fgets(char_crap, 100, mesh_file_);
}
//****************************************************************************80
UnstMeshReaderGMSH::file_type UnstMeshReaderGMSH::ReadHeader()
{
  SystemModule::cout<< "Reading Header" << std::endl;
  file_type result;
  //---> Read crap that is the first line
  SkipLine();
  realT version;
  intT type_flag;
  intT real_size_bytes;
  intT icrap;
  fscanf(mesh_file_, "%lf %d %d\n", &version, &type_flag, &real_size_bytes);
  if(type_flag == 0){
    SystemModule::cout << "Detected GMSH ASCII File Version:  "
        << version << std::endl;
    result = file_type::ASCII;
  }
  else if(type_flag == 1){
    fread(&icrap, sizeof(intT), 1, mesh_file_);
    fscanf(mesh_file_, "\n");
    SystemModule::cout << "Detected GMSH BINARY File Version:  "
        << version << std::endl;
   result = file_type::BINARY;
  }
  else{
    SystemModule::cout<< "ERROR: GMSH File header did not specify ASCII or BINARY." << std::endl
        << " The type value read in is: " << type_flag << ", which should be 0 = ASCII or "
        << " 1 = BINARY. The file name is: " << UnstMeshReader::filename_
        << std::endl;
    SystemModule::my_exit();
  }

  return result;
}
//****************************************************************************80
void UnstMeshReaderGMSH::ReadNodesASCII(read_mode mode)
{

  intT nnode;
  switch (mode) {
    case UnstMeshReaderGMSH::STORE :
    //---> Read number of nodes
    fscanf(mesh_file_, "%d \n", &nnode);
    UnstMeshReader::nnode_= nnode;

    x_.initialize(nnode, 3);
    for(intT n = 0; n < nnode; n++){
      intT node;
      realT x;
      realT y;
      realT z;
      fscanf(mesh_file_, "%d %lf %lf %lf\n", &node, &x, &y, &z);
      x_(node - 1,0) = x;
      x_(node - 1,1) = y;
      x_(node - 1,2) = z;
    }
    break;
    case UnstMeshReaderGMSH::SKIP :
      //---> Read number of nodes;
      fscanf(mesh_file_, "%d \n", &nnode);
      for(intT n = 0; n < nnode; n++){
        intT node;
        realT x;
        realT y;
        realT z;
        fscanf(mesh_file_, "%d %lf %lf %lf\n", &node, &x, &y, &z);
      }
      break;
  }

}// End Function ReadNodesASCII
//****************************************************************************80
void UnstMeshReaderGMSH::ReadNodesBinary(read_mode mode)
{

  intT nnode;
  realT xyz[3];
  switch (mode) {
    case UnstMeshReaderGMSH::STORE :
      //---> Read number of nodes
      fscanf(mesh_file_, "%d \n", &nnode);
      UnstMeshReader::nnode_= nnode;
      x_.initialize(nnode,3);
      for(intT n = 0; n < nnode; n++){
        intT node;
        fread(&node, sizeof(intT), 1, mesh_file_);
        fread(xyz, sizeof(realT), 3, mesh_file_);
        x_(node - 1,0) = xyz[0];
        x_(node - 1,1) = xyz[1];
        x_(node - 1,2) = xyz[2];
      }

      break;
    case UnstMeshReaderGMSH::SKIP :
      //---> Read number of nodes;
      fscanf(mesh_file_, "%d \n", &nnode);
      for(intT n = 0; n < nnode; n++){
        intT node;
        fread(&node, sizeof(intT), 1, mesh_file_);
        fread(xyz, sizeof(realT), 3, mesh_file_);
      }
      break;
  }
  fscanf(mesh_file_, "\n");
}// End Function ReadNodesBindary
//****************************************************************************80
void UnstMeshReaderGMSH::CountElementsASCII()
{
  //---> Read Elements
  fscanf(mesh_file_, "%d \n", &nentity_);

  entity_type_.initialize(nentity_);
  entity_id_.initialize(nentity_);
  Array1D<intT> nnode_per_entity(nentity_);

  for(intT e = 0; e < nentity_; e++){ // Element Loop
    intT elem;
    intT type;
    intT ntag;
    intT id;
    intT icrap;
    fscanf(mesh_file_,"%d %d %d %d", &elem, &type, &ntag, &id);
    //---> Skip 1 to ntag -1 columns
    entity_type_(elem - 1) = type;
    entity_id_(elem - 1) = id;
    for(intT i = 1; i < ntag; i++){fscanf(mesh_file_,"%d",&icrap);}
    intT n;
    switch(gmsh_type_map_[type]) {
      case ElementTopology::element_types::BAR :
        n = ElementTopology::Bar::nNode;
        break;
      case ElementTopology::element_types::TRI :
        n = ElementTopology::Triangle::nNode;
        break;
      case ElementTopology::element_types::QUAD :
        n = ElementTopology::Quadrilateral::nNode;
        break;
      case ElementTopology::element_types::TET :
        n = ElementTopology::Tetrahedron::nNode;
        break;
      case ElementTopology::element_types::PRISM :
        n = ElementTopology::Prism::nNode;
        break;
      case ElementTopology::element_types::PYR :
        n = ElementTopology::Pyramid::nNode;
        break;
      case ElementTopology::element_types::HEX :
        n = ElementTopology::Hexahedron::nNode;
        break;
    }
    for(intT i = 0; i < n; i++){fscanf(mesh_file_,"%d",&icrap);}
    fscanf(mesh_file_,"\n");

    nnode_per_entity(elem-1) = n;
    if(gmsh_id_is_bc_[id]){
      UnstMeshReader::nbc_face_++;
    }
    if(gmsh_id_is_region_[id]){
      UnstMeshReader::nelement_++;
    }

  }// End Element Loop

  entity2node_.initialize(nnode_per_entity);
  bc_face2entity_.initialize(UnstMeshReader::nbc_face_);
  elem2entity_.initialize(UnstMeshReader::nelement_);
} // End  UnstMeshReaderGMSH::CountElementsASCII()
//****************************************************************************80
void UnstMeshReaderGMSH::ReadElementsASCII()
{
  //---> Read Elements
  intT icrap;
  fscanf(mesh_file_, "%d \n", &icrap);

  intT ielem = 0;
  intT ibc_face = 0;
  for(intT e = 0; e < nentity_; e++){ // Element Loop
    intT elem;
    intT type;
    intT ntag;
    intT id;
    fscanf(mesh_file_,"%d %d %d %d", &elem, &type, &ntag, &id);

    for(intT i = 1; i < ntag; i++){fscanf(mesh_file_,"%d",&icrap);}

    intT n = entity2node_.get_ncol(elem-1);

    for(intT i = 0; i < n; i++){
      fscanf(mesh_file_,"%d",&icrap);
      entity2node_(elem - 1, i) = icrap - 1;
    }

    fscanf(mesh_file_,"\n");
    if(gmsh_id_is_bc_[id]){
      bc_face2entity_(ibc_face) = elem - 1;
      ibc_face++;
    }
    if(gmsh_id_is_region_[id]){
      elem2entity_(ielem) = elem - 1;
      ielem++;
    }
  }// End Element Loop
}// End UnstMeshReaderGMSH::ReadElementsASCII()
//****************************************************************************80
void UnstMeshReaderGMSH::CountElementsBinary()
{
  //---> Read Elements
  fscanf(mesh_file_, "%d \n", &nentity_);

  entity_type_.initialize(nentity_);
  entity_id_.initialize(nentity_);
  Array1D<intT> nnode_per_entity(nentity_);

  intT ientity = 0;
  while(ientity < nentity_){
    int elem_header[3]={0,100,200};
    fread(elem_header, sizeof(intT), 3, mesh_file_);

    intT type = elem_header[0];
    ientity += elem_header[1];
    intT ntag = elem_header[2];
    intT n;
    switch(gmsh_type_map_[type]) {
      case ElementTopology::element_types::BAR :
        n = ElementTopology::Bar::nNode;
        break;
      case ElementTopology::element_types::TRI :
        n = ElementTopology::Triangle::nNode;
        break;
      case ElementTopology::element_types::QUAD :
        n = ElementTopology::Quadrilateral::nNode;
        break;
      case ElementTopology::element_types::TET :
        n = ElementTopology::Tetrahedron::nNode;
        break;
      case ElementTopology::element_types::PRISM :
        n = ElementTopology::Prism::nNode;
        break;
      case ElementTopology::element_types::PYR :
        n = ElementTopology::Pyramid::nNode;
        break;
      case ElementTopology::element_types::HEX :
        n = ElementTopology::Hexahedron::nNode;
        break;
    }

    intT ndata = 1 + ntag + n;
    intT elem_data[ndata];
    for(intT i = 0; i < elem_header[1]; i++){ // Element Loop
      fread(elem_data, sizeof(intT), ndata, mesh_file_);

      intT elem = elem_data[0];
      intT id = elem_data[1];
      //---> Skip 1 to ntag -1 columns
      entity_type_(elem - 1) = type;
      entity_id_(elem - 1) = id;
      nnode_per_entity(elem - 1) = n;
      if(gmsh_id_is_bc_[id]){

        UnstMeshReader::nbc_face_++;
      }
      if(gmsh_id_is_region_[id]){

        UnstMeshReader::nelement_++;
      }

    }// End Element Loop


  } // End While

  entity2node_.initialize(nnode_per_entity);
  bc_face2entity_.initialize(UnstMeshReader::nbc_face_);
  elem2entity_.initialize(UnstMeshReader::nelement_);
}
//****************************************************************************80
void UnstMeshReaderGMSH::ReadElementsBinary()
{
  //---> Read Elements
  fscanf(mesh_file_, "%d \n", &nentity_);

  intT ientity = 0;
  intT ibc_face = 0;
  intT ielem = 0;
  int elem_header[3];
  while(ientity < nentity_){
    fread(elem_header, sizeof(intT), 3, mesh_file_);

    intT type = elem_header[0];
    ientity += elem_header[1];
    intT ntag = elem_header[2];
    intT n;
    switch(gmsh_type_map_[type]) {
      case ElementTopology::element_types::BAR :
        n = ElementTopology::Bar::nNode;
        break;
      case ElementTopology::element_types::TRI :
        n = ElementTopology::Triangle::nNode;
        break;
      case ElementTopology::element_types::QUAD :
        n = ElementTopology::Quadrilateral::nNode;
        break;
      case ElementTopology::element_types::TET :
        n = ElementTopology::Tetrahedron::nNode;
        break;
      case ElementTopology::element_types::PRISM :
        n = ElementTopology::Prism::nNode;
        break;
      case ElementTopology::element_types::PYR :
        n = ElementTopology::Pyramid::nNode;
        break;
      case ElementTopology::element_types::HEX :
        n = ElementTopology::Hexahedron::nNode;
        break;
    }

    intT ndata = 1 + ntag + n;
    intT elem_data[ndata];
    for(intT i = 0; i < elem_header[1]; i++){ // Element Loop
      fread(elem_data, sizeof(intT), ndata, mesh_file_);

      intT elem = elem_data[0];
      intT id = elem_data[1];
      for(intT j = 0; j < n; j++){
        entity2node_(elem-1,j) = elem_data[1 + ntag + j] - 1;
      }
      if(gmsh_id_is_bc_[id]){
        bc_face2entity_(ibc_face) = elem - 1;
        ibc_face++;
      }
      if(gmsh_id_is_region_[id]){
        elem2entity_(ielem) = elem - 1;
        ielem++;
      }

    }// End Element Loop


  } // End While

}
