//-*-c++-*-
#ifndef ELEMENTTOPOLOGY_H
#define ELEMENTTOPOLOGY_H
#include "my_incl.h"
#include "DataStructures/StaticList2D.h"

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
    FACE_EDGE = 1,
    FACE_TRI = 2,
    FACE_QUAD = 3,
  };

//****************************************************************************80
//! \brief Bar : Describes the Topology of a bar element 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  class Bar 
  {
  public:
    static const intT Type;
    static const intT nNode;
    static const intT nEdge;
    static const intT nFace;
    static const intT Nodes[2];// = {0,1};
    static const intT Edges[1][2]; //= {0,1};
    static const StaticList2D<intT,2,2> Faces;

  }; // End Class Bar
//****************************************************************************80
//! \brief Triangle : Describes the Topology of a triangle element 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************8
  class Triangle 
  {
  public:
    static const intT Type;
    static const intT nNode;
    static const intT nEdge;
    static const intT nFace;
    static const intT Nodes[3];
    static const intT Edges[3][2];
    static const StaticList2D<intT, 3, 6> Faces;
  }; // End Class Triangle

//****************************************************************************80
//! \brief Quadrilateral : Describes the Topology of a quadrilateral element 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************8
  class Quadrilateral
  {
    static const intT Type;
    static const intT nNode;
    static const intT nEdge;
    static const intT nFace;
    static const intT Nodes[4];
    static const intT Edges[4][2];
    static const StaticList2D<intT, 4, 8> Faces;
  }; // End Class Quadrilateral

//****************************************************************************80
//! \brief Tetrahedra : Describes the Topology of a tetrahedra element 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************8
  class Tetrahedron
  {
    static const intT Type;
    static const intT nNode;
    static const intT nEdge;
    static const intT nFace;
    static const intT Nodes[4];
    static const intT Edges[6][2];
    static const StaticList2D<intT, 4, 12> Faces;
  };// End Class Tetrahedron
//****************************************************************************80
//! \brief Prism : Describes the Topology of a prism element 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************8
  class Prism
  {
    static const intT Type;
    static const intT nNode;
    static const intT nEdge;
    static const intT nFace;
    static const intT Nodes[6];
    static const intT Edges[9][2];
    static const StaticList2D<intT, 5, 18> Faces;
  };// End Class Prism
//****************************************************************************80
//! \brief Pyramid : Describes the Topology of a pyramid element 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************8
  class Pyramid
  {
    static const intT Type;
    static const intT nNode;
    static const intT nEdge;
    static const intT nFace;
    static const intT Nodes[5];
    static const intT Edges[8][2];
    static const StaticList2D<intT, 5, 16> Faces; 
  };// End Pyramid
//****************************************************************************80
//! \brief Hexahedron : Describes the Topology of a hexadedron 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************8
  class Hexahedron
  {
    static const intT Type;
    static const intT nNode;
    static const intT nEdge;
    static const intT nFace;
    static const intT Nodes[8];
    static const intT Edges[12][2];
    static const StaticList2D<intT, 6, 24> Faces;
  }; // End Hexahedron
}// End namespace ElementTopology
#endif
