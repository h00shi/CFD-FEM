//-*-c++-*-
#ifndef ELEMENT_TOPOLOGY_TRI_H
#define ELEMENT_TOPOLOGY_TRI_H
#include "DataStructures/StaticList2D.h"

namespace ElementTopology
{
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
}// End Namespace ElementTopology
#endif
