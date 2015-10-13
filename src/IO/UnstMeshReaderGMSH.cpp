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
  SystemModule::cout<< "Reading Header" << std::endl;
  //---> Read crap that is the first line
  char char_crap[100];
  fgets(char_crap, 12, mesh_file_);
  realT version;
  intT type_flag;
  intT real_size_bytes;
  fscanf(mesh_file_, "%lf %d %d\n", &version, &type_flag, &real_size_bytes);
  if(type_flag == 0){
    SystemModule::cout << "Detected GMSH ASCII File Version:  " << version << std::endl;
    ReadASCII();}
  else if(type_flag == 1){ReadBinary();}
  else{
    SystemModule::cout<< "ERROR: GMSH File header did not specify ASCII or BINARY." << std::endl
        << " The type value read in is: " << type_flag << ", which should be 0 = ASCII or "
        << " 1 = BINARY. The file name is: " << UnstMeshReader::filename_
        << std::endl;
    SystemModule::my_exit();
  }
  fclose(mesh_file_);
}// End UnstMeshReaderGMSH::UnstMeshReaderGMSH
//****************************************************************************80
void UnstMeshReaderGMSH::ReadASCII()
{
  //---> We know we are reading a 2-D mesh file;
  intT ndim = 3;

  char char_crap[100];
  intT nnode;

  //---> Two lines of Crap
  fgets(char_crap, 100, mesh_file_);
  fgets(char_crap, 100, mesh_file_);
  //---> Read number of nodes;
  fscanf(mesh_file_, "%d \n", &nnode);
  UnstMeshReader::nnode_= nnode;

  x_.initialize(nnode, ndim);
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
  //---> Two more lines of crap
  fgets(char_crap, 100, mesh_file_);
  fgets(char_crap, 100, mesh_file_);
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
    if(global_id_2_bc_id_(id)){
      UnstMeshReader::nbc_face_++;
    }
    if(global_id_2_region_id_(id)){
      UnstMeshReader::nelement_++;
    }

  }// End Element Loop

  entity2node_.initialize(nnode_per_entity);
  bc_face2entity_.initialize(UnstMeshReader::nbc_face_);
  elem2entity_.initialize(UnstMeshReader::nelement_);

  //std::cout << nnode_per_entity << std::endl;

  //---> Rewind file;
  rewind(mesh_file_);
  fgets(char_crap, 12, mesh_file_);
  realT rcrap;
  intT icrap;
  fscanf(mesh_file_, "%lf %d %d\n", &rcrap, &icrap, &icrap);
  //---> Two lines of Crap
  fgets(char_crap, 100, mesh_file_);
  fgets(char_crap, 100, mesh_file_);
  //---> Read number of nodes;
  fscanf(mesh_file_, "%d \n", &icrap);

  for(intT n = 0; n < nnode; n++){
    fscanf(mesh_file_, "%d %lf %lf %lf\n", &icrap, &rcrap, &rcrap, &rcrap);
  }
  //---> Two more lines of crap
  fgets(char_crap, 100, mesh_file_);
  fgets(char_crap, 100, mesh_file_);
  //---> Read Elements
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
    if(global_id_2_bc_id_(id)){
      bc_face2entity_(ibc_face) = elem - 1;
      ibc_face++;
    }
    if(global_id_2_region_id_(id)){
      elem2entity_(ielem) = elem - 1;
      ielem++;
    }
  }// End Element Loop

  SystemModule::cout << "Detected the following mesh file Parameters: " << std::endl
                     << "nElement: " << ielem << std::endl
                     << "nBcFace: " << ibc_face << std::endl;


}// End ReadASCII
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
  fscanf(idmap_file_, "%d \n", &(UnstMeshReader::nbc_id_));
  bc_id_2_global_id_.initialize(UnstMeshReader::nbc_id_);
  fgets(char_crap, 100, idmap_file_);
  for(intT i = 0; i < UnstMeshReader::nbc_id_; i++){
    intT n1;
    intT n2;
    fscanf(idmap_file_,"%d %d \n", &n1, &n2);
    bc_id_2_global_id_(n2) = n1;
  }

  fgets(char_crap, 100, idmap_file_);
  fscanf(idmap_file_, "%d \n", &n_region_id_);
  region_id_2_global_id_.initialize(n_region_id_);
  fgets(char_crap, 100, idmap_file_);
  for(intT i = 0; i < n_region_id_; i++){
    intT n1;
    intT n2;
    fscanf(idmap_file_,"%d %d \n", &n1, &n2);
    region_id_2_global_id_(n2) = n1;
  }

  global_id_2_bc_id_.initialize(UnstMeshReader::nbc_id_ + n_region_id_);
  global_id_2_region_id_.initialize(UnstMeshReader::nbc_id_ + n_region_id_);
  gmsh_id_2_my_id_.initialize(UnstMeshReader::nbc_id_ + n_region_id_);
  for(intT i = 0; i < UnstMeshReader::nbc_id_; i++){
    intT gid = bc_id_2_global_id_(i);
    gmsh_id_2_my_id_(gid) = i;
    global_id_2_bc_id_(gid) = true;
    global_id_2_region_id_(gid) = false;
  }
  for(intT i = 0; i < n_region_id_; i++){
    intT gid = region_id_2_global_id_(i);
    gmsh_id_2_my_id_(gid) = i;
    global_id_2_bc_id_(gid) = false;
    global_id_2_region_id_(gid) = true;
  }

  SystemModule::cout << "Mapping the following GMSH id numbers: " << std::endl;
  for(intT i = 0; i < UnstMeshReader::nbc_id_ + n_region_id_; i++){
   std::cout << "GMSH ID: " << i << " Maps to ID: " << gmsh_id_2_my_id_(i)
             << std::endl;
  }

  SystemModule::cout << "Done Reading IDMAP File." << std::endl;
  fclose(idmap_file_);
}// End UnstMeshReaderGMSH::ReadIdMap

//****************************************************************************80
void UnstMeshReaderGMSH::ReadBinary()
{

  return;
}
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
     ereg_tmp(e) = gmsh_id_2_my_id_(entity_type_(entity));
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
     bcid_tmp(f) = gmsh_id_2_my_id_(entity_id_(entity));
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
