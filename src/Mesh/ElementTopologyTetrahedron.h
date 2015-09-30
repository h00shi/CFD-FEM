//-*-c++-*-
#ifndef ELEMENT_TOPOLOGY_TET_H
#define ELEMENT_TOPOLOGY_TET_H
#include "DataStructures/StaticList2D.h"

namespace ElementTopology
{
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
}// End Namespace ElementTopology
#endif
