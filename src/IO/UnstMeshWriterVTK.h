//-*-c++-*-

#ifndef UNSTGRIDVTKWRITER_H
#define UNSTGRIDVTKWRITER_H
//****************************************************************************80
//! \brief Class that implements writing a mesh to VTK Library
//! \details Makes calls to the library APIs to write meshes to VTK formats
//! \qnick
//! \version $Rev$
//!
//****************************************************************************80
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkLine.h>
#include <vtkTriangle.h>
#include <vtkQuad.h>
#include <vtkTetra.h>
#include <vtkHexahedron.h>
#include <vtkWedge.h>
#include <vtkPyramid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include "Mesh/UnstMesh.h"
#include "IO/UnstMeshWriter.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/List2D.h"
#include "ParallelComm/Communication.h"
#include <cassert>
#include <iostream>
#include <string>
#include <cstdlib>

class UnstMeshWriterVTK : public UnstMeshWriter 
{
public:

//****************************************************************************80
//! \brief UnstMeshWriterVTK : Constructor taking in reference to a grid
//! \details
//! \nick
//! \version $Rev$
//! \param[in] mesh_ref Reference to a grid that you want to write
//****************************************************************************80
  UnstMeshWriterVTK(const UnstMesh& mesh_ref);

//****************************************************************************80
//! \brief ~UnstMeshWriterVTK : Destructor
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  ~UnstMeshWriterVTK();
//****************************************************************************80
//! \brief Write : writes the mesh out to the specified file...find
//! \details
//! \nick
//! \version $Rev$
//!
//****************************************************************************80
  void Write(const std::string& fbase);

private:
//****************************************************************************80
//! \brief UnstMeshWriterVTK : Default constructor
//! \details BLOCK THIS! We don't want to instantiate a writer without a
//!                      reference to a mesh
//! \nick
//! \version $Rev$
//!
//****************************************************************************80
  UnstMeshWriterVTK() = delete;

private:
  vtkSmartPointer< vtkUnstructuredGrid > vtk_mesh_; /*!< VTK mesh datastructure
						      required to write */

};// End Class UnstMeshWriterVTK

#endif
