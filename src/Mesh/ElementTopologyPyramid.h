//-*-c++-*-
#ifndef ELEMENT_TOPOLOGY_PYR_H
#define ELEMENT_TOPOLOGY_PYR_H
#include "DataStructures/StaticList2D.h"

namespace ElementTopology
{
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
  public:
    static const intT Type;
    static const intT nNode;
    static const intT nEdge;
    static const intT nFace;
    static const intT Nodes[5];
    static const intT Edges[8][2];
    static const StaticList2D<intT, 5, 16> Faces; 
  };// End Pyramid
}// End namespace ElementTopology
#endif
