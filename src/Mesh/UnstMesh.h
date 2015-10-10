// -*-c++-*-
#ifndef UNSTMESH_H
#define UNSTMESH_H
//****************************************************************************80
//! \class UnstMesh 
//! \brief This is the header file defining the class UnstMesh
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//****************************************************************************80
#include "my_incl.h"
#include "IO/UnstMeshReader.h"
#include "SystemUtils/SystemModule.h"
#include "Mesh/ElementTopology.h"
#include "Mesh/UnstMeshGeom.h"
#include "Mesh/UnstMeshElements.h"
#include "Mesh/UnstMeshBcFaces.h"

class UnstMesh
{
  
public:
//****************************************************************************80
//! \brief  Constructor using a mesh reader
//! \details  Constructor that creates a UnstMesh object with a file input
//! \nick
//! \param[in] mesh_reader The object that reads the mesh
//****************************************************************************80
  UnstMesh(UnstMeshReader& mesh_reader);

//****************************************************************************80
//! \brief The destructor
//! \details  Destructor is pure virtual to make this class abstract
//! \nick
//! \version $Rev$
//****************************************************************************80
  virtual ~UnstMesh()=0;
  inline const UnstMeshGeom& get_MeshGeom() const {return mesh_geometry_;}
  inline const UnstMeshElements& get_MeshElements() const
    {return mesh_elements_;}
  inline const UnstMeshBcFaces& get_MeshBcFaces() const {return mesh_bc_faces_;}

//****************************************************************************80
//! \brief MemoryDiagnostic : Runs a full diagnostic of the memory for the 
//!        class.  Will print all data to standard out 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  void Diagnostic(std::ostream& out_stream);

protected: 

  //++++++++++++++++++++++++++++++ PROTECTED STUFF +++++++++++++++++++++++++++++
  UnstMeshGeom mesh_geometry_;
  UnstMeshElements mesh_elements_;
  UnstMeshBcFaces mesh_bc_faces_;

  //---> Memory diagnostics data
  realT grid_mem_; /*!< Total memory for this instance of class */

private:

//****************************************************************************80
//! \brief  Is the copy constructor
//! \details  The copy constructor is explicitly blocked.
//! \nick
//! \version
//****************************************************************************80
  UnstMesh(const UnstMesh&) = delete;

//****************************************************************************80
//! \brief operator= : Is the default assignment operator for this class.
//!                    This is blocked to prevent assignment of meshes
//! \details  The assignment operator is explicitly blocked.
//! \nick
//! \version
//****************************************************************************80
  UnstMesh& operator= (const UnstMesh&) = delete;
  
};
#endif
