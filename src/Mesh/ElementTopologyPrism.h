//-*-c++-*-
#ifndef ELEMENT_TOPOLOGY_PRISM_H
#define ELEMENT_TOPOLOGY_PRISM_H
#include "DataStructures/StaticList2D.h"

namespace ElementTopology
{
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
  public:
    static const intT Type;
    static const intT nNode;
    static const intT nEdge;
    static const intT nFace;
    static const intT Nodes[6];
    static const intT Edges[9][2];
    static const StaticList2D<intT, 5, 18> Faces;
  };// End Class Prism
}// End namespace ElementTopology
#endif
