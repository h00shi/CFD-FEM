/*
 * UnstMeshReaderGMSH.h
 *
 *  Created on: Oct 11, 2015
 *      Author: rabbit
 */

#ifndef UNSTMESHREADERGMSH_H_
#define UNSTMESHREADERGMSH_H_
#include "IO/UnstMeshReader.h"
#include "Mesh/ElementTopology.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/List2D.h"
#include <cstring>
#include <unordered_map>
//****************************************************************************80
//! \brief This A class to read GMSH files.
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//!
//****************************************************************************80
class UnstMeshReaderGMSH : public UnstMeshReader
{
public:

//****************************************************************************80
//! \brief Constructor
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[in] filename The name of the file
//! \param[in] type The type of file ASCII or BINARY
//****************************************************************************80
  UnstMeshReaderGMSH(const std::string& filename,
                     const std::string& imap_filename);
  Array2D<realT> ReadNodes();
  List2D<intT>  ReadElement2Node();
  Array1D<ElementTopology::element_types> ReadElementType();
  Array1D<intT> ReadElementRegion();
  List2D<intT>  ReadBcFace2Node();
  Array1D<intT> ReadBcID();
  Array1D<ElementTopology::face_types> ReadBcFaceType();
private:
  intT nentity_;
  intT n_region_id_;
  intT gmsh_type_map_[8] = {-1,
                           ElementTopology::element_types::BAR,
                           ElementTopology::element_types::TRI,
                           ElementTopology::element_types::QUAD,
                           ElementTopology::element_types::TET,
                           ElementTopology::element_types::HEX,
                           ElementTopology::element_types::PRISM,
                           ElementTopology::element_types::PYR};
  FILE* mesh_file_;
  FILE* idmap_file_;
  Array2D<realT> x_;
  List2D<intT> entity2node_;
  Array1D<intT> entity_type_;
  Array1D<intT> entity_id_;
  Array1D<intT> elem2entity_;
  Array1D<intT> bc_face2entity_;
  std::unordered_map<intT, intT> gmsh_id_2_bc_id_;
  std::unordered_map<intT, intT> gmsh_id_2_region_;
  std::unordered_map<intT, bool> gmsh_id_is_bc_;
  std::unordered_map<intT, bool> gmsh_id_is_region_;
  Array1D<intT> gmsh_id_;//!< Mapping an index to a gmsh_id


  //---> Reads file if it's ASCII
  void ReadASCII();
  //---> Reads file it it's BINARY
  void ReadBinary();
  //---> Read the IDMAP file
  void ReadIdMap(const std::string& idmap_filename);
  UnstMeshReaderGMSH() = delete;
};



#endif /* UNSTMESHREADERGMSH_H_ */
