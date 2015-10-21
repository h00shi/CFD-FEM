/*
 * CGMesh.h
 *
 *  Created on: Oct 10, 2015
 *      Author: rabbit
 */

#ifndef CGMESH_H_
#define CGMESH_H_
#include "Mesh/UnstMesh.h"
#include "DataStructures/Graph.h"
#include "DataStructures/Array1D.h"
#include "IO/UnstMeshReader.h"
//****************************************************************************80
//! \class CGMesh
//! \brief An unstructured mesh with the proper connectivities for
//!        CG(Continuous Galerkin) methods.
//! \nick
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//****************************************************************************80
class CGMesh : public UnstMesh
{
public:
//****************************************************************************80
//! \brief  Constructor using a mesh reader
//! \details  Constructor that creates a UnstMesh object with a file input
//! \nick
//! \param[in] mesh_reader The object that reads the mesh
//****************************************************************************80
  CGMesh(UnstMeshReader& mesh_reader);
//****************************************************************************80
//! \brief  Destructor
//! \details  Constructor that creates a UnstMesh object with a file input
//! \nick
//****************************************************************************80
  ~CGMesh();
  inline const Graph& get_Graph() const {return *graph_;}
  inline const Array1D<intT>& get_ElemOrder() const {return elem_p_;}

private:
  Graph* graph_;
  Array1D<intT> elem_p_;
//****************************************************************************80
//! \brief Forms the graph for CGMeshes, which is just the node2node
//!        connectivity
//! \details  Constructor that creates a UnstMesh object with a file input
//! \nick
//****************************************************************************80
  void FormGraph();

  //---> Block Default construction
  CGMesh() = delete;
};




#endif /* CGMESH_H_ */
