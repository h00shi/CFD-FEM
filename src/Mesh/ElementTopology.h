//-*-c++-*-
#ifndef ELEMENTTOPOLOGY_H
#define ELEMENTTOPOLOGY_H
#include "my_incl.h"
#include "DataStructures/StaticList2D.h"
#include <vector>
//****************************************************************************80
//! \brief A namespace that describes element topology for each type of 
//!        supported element.
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
namespace ElementTopology
{
  /*! ELEMENT TYPES
      0 = 1-D Bar, a line segement with end points
      1 = 2-D Triangle, 3 verticies connected together
      2 = 2-D Quad, 4 verticies connected together
      3 = 3-D Tetrahedra, 4 verticies with 4 triangular faces
      4 = 3-D Prism, 6 vertices with 2 trianglular faces and
      3 quadrilateral faces
      5 = 3-D Pyramid, 5 verticies with 4 triangular faces and
      1 quadrilateral face
      6 = 3-D Hexahedra, 8 verticies with 6 quadrilateral faces
  */
  enum element_types: intT{
    BAR = 0,
    TRI = 1,
    QUAD = 2,
    TET = 3,
    PRISM = 4,
    PYR = 5,
    HEX = 6
  };

  /*! FACE TYPES
      0 = Node
      1 = Edge
      2 = Triangle
      3 = Quad
  */
  enum face_types: intT{
    FACE_NODE = 0,
    FACE_BAR = 1,
    FACE_TRI = 2,
    FACE_QUAD = 3,
  };

  //****************************************************************************80
  //! \brief A function template that tell you what the local face described
  //!        by face_nodes is of the element described by elem_nodes
  //!
  //! \details
  //! \nick
  //! \version $Rev$
  //! \date $Date$
  //! \\param[in] face_nodes Set of nodes describing a face
  //! \\param[in] elem_nodes Set of nodes describing an element
  // return local_face
  //****************************************************************************80
  template <class Topology>
  intT FindLocalFace(const std::vector<intT>& face_nodes,
                     const std::vector<intT>& elem_nodes);

  //****************************************************************************80
  template <class Topology>
  intT FindLocalFace(const std::vector<intT>& face_nodes,
                      const std::vector<intT>& elem_nodes)
  {
    intT local_face = -1;
    for(intT i = 0; i < Topology::nFace; i++){
      //---> Load nodes from element into temp vector
      std::vector<intT> eface_nodes;
      eface_nodes.reserve(Topology::Faces.get_ncol(i));
      for(intT j = 0; j < Topology::Faces.get_ncol(j); j++){
        intT node = elem_nodes[Topology::Faces(i,j)];
        eface_nodes.push_back(node);
      }

      //---> Use comparison function
      if(face_nodes == eface_nodes){local_face = i; break;}
    }
    return local_face;

  }

}// End namespace ElementTopology

#include "Mesh/ElementTopologyBar.h"
#include "Mesh/ElementTopologyTriangle.h"
#include "Mesh/ElementTopologyQuadrilateral.h"
#include "Mesh/ElementTopologyTetrahedron.h"
#include "Mesh/ElementTopologyPrism.h"
#include "Mesh/ElementTopologyPyramid.h"
#include "Mesh/ElementTopologyHexahedron.h"
#endif
