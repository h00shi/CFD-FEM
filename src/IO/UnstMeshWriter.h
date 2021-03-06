// -*-c++-*-
#ifndef UNSTMESHWRITER_H
#define UNSTMESHWRITER_H

//****************************************************************************80
//! \brief This is the base class for writing unstructured meshes to file
//! \details This contains of APIs and for writing meshes.
//!   If you'd like to add a writer inherit from this base and follow
//!   the API interfaces already present.
//! \nick
//! \version $Rev$
//!
//****************************************************************************80
#include "my_incl.h"
#include "Mesh/UnstMesh.h"

class UnstMeshWriter {

public:

//****************************************************************************80
//! \brief UnstMeshWriter : Constructor taking in reference to a grid
//! \details
//! \nick
//! \version $Rev$
//! \param[in] mesh Reference to a mesh that you want to write
//****************************************************************************80
  UnstMeshWriter(const UnstMesh& mesh_ref) : mesh_(mesh_ref) {}

//****************************************************************************80
//! \brief ~UnstMeshWriter : Destructor
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  virtual ~UnstMeshWriter() { }

//****************************************************************************80
//! \brief Write : writes the mesh out to the specified file...find
//! \details
//! \nick
//! \version $Rev$
//!
//****************************************************************************80
  virtual void Write(const std::string& fbase) = 0;

protected:
  UnstMesh const & mesh_; //!< Reference to UnstMesh

private:
//****************************************************************************80
//! \brief UnstMeshWriter : Default constructor
//! \details BLOCK THIS! We don't want to instantiate a writer without a
//!                      reference to a mesh
//! \nick
//! \version $Rev$
//!
//****************************************************************************80
  UnstMeshWriter() = delete;


};


#endif
