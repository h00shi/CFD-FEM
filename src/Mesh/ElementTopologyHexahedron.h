//-*-c++-*-
#ifndef ELEMENT_TOPOLOGY_HEX_H
#define ELEMENT_TOPOLOGY_HEX_H
#include "DataStructures/StaticList2D.h"

namespace ElementTopology
{
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
