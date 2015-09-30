//-*-c++-*-
#ifndef ELEMENT_TOPOLOGY_BAR_H
#define ELEMENT_TOPOLOGY_BAR_H
#include "DataStructures/StaticList2D.h"
namespace ElementTopology
{
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
}// End Namespace ElementTopology
#endif
