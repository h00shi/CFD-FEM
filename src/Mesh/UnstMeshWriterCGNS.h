// -*-c++-*-

#ifndef UNSTMESHWRITERCGNS_H
#define UNSTMESHWRITERCGNS_H
//****************************************************************************80
//! \brief This is the class for writing meshes using CGNS library
//! \details Makes calls to the library APIs to write meshes to CGNS formats
//! \nick
//! \version $Rev$
//!
//****************************************************************************80
#include "my_incl.h"
#include "Array2D.h"
#include "Array1D.h"
#include "List2D.h"
#include "UnstMeshWriter.h"
#include "cgnslib.h"

class UnstMeshWriterCGNS : public UnstMeshWriter
{
private:

  Array1D<realT> x_; //!< CNGS X-coordinate vector
  Array1D<realT> y_; //!< CNGS Y-coordinate vector
  Array1D<realT> z_; //!< CNGS Z-coordinate vector
  Array2D<cgsize_t> bar_con_; //!< Connectivity of bar elements
  Array2D<cgsize_t> tri_con_; //!< Connectivity of tri elements
  Array2D<cgsize_t> quad_con_; //!< Connectivity of quad elements
  Array2D<cgsize_t> tet_con_; //!< Connectivity of tet elements
  Array2D<cgsize_t> prism_con_; //!< Connectivity of prism elements
  Array2D<cgsize_t> pyr_con_; //!< Connectivity of pyramid elements
  Array2D<cgsize_t> hex_con_; //!< Connectivity of hex elements
  Array2D<cgsize_t> bc_node_con_; //!< Connectivity of boundary nodes
  Array2D<cgsize_t> bc_edge_con_; //!< Connectivity of bc edges
  Array2D<cgsize_t> bc_tri_con_; //!< Connectivity of bc triangles
  Array2D<cgsize_t> bc_quad_con_; //!< Connectivity of bc quads
  List2D<cgsize_t> bc_id_edge_; /*!< list of boundary faces that are edge 
				  by bcid */
  List2D<cgsize_t> bc_id_tri_; //!< List of boundary faces that are tri by bcid
  List2D<cgsize_t> bc_id_quad_;//!< List of boundary faces that are quad by bcid
//****************************************************************************80
//! \brief UnstMeshWriterCGNS : Default constructor
//! \details BLOCK THIS! We don't want to instantiate a writer without a
//!                      reference to a mesh
//! \nick
//! \version $Rev$
//!
//****************************************************************************80
  UnstMeshWriterCGNS() = delete;

public:

//****************************************************************************80
//! \brief UnstMeshWriterCGNS : Default constructor
//! \details BLOCK THIS! We don't want to instantiate a writer without a
//!                      reference to a mesh
//! \nick
//! \version $Rev$
//! \param[in] grid_ref The reference to the grid you want to write
//****************************************************************************80
  UnstMeshWriterCGNS(const UnstMesh& grid_ref);

//****************************************************************************80
//! \brief ~UnstMeshWriterCGNS : Destructor
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  ~UnstMeshWriterCGNS();

//****************************************************************************80
//! \brief Write : Writes a cgns formatted file of the grid
//! \details
//! \nick
//! \version $Rev$
//! \param[in] Project String giving the project name
//****************************************************************************80
  void Write(const std::string& fbase );
  
};
#endif
