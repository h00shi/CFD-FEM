//-*-c++-*-
#ifndef ELEMENT_TOPOLOGY_QUAD_H
#define ELEMENT_TOPOLOGY_QUAD_H
#include "DataStructures/StaticList2D.h"

namespace ElementTopology
{

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
}
#endif
