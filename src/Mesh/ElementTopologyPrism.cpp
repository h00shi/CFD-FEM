#include "Mesh/ElementTopology.h"
//****************************************************************************80
const intT ElementTopology::Prism::Type = ElementTopology::element_types::PRISM;
const intT ElementTopology::Prism::nNode = 6;
const intT ElementTopology::Prism::nEdge = 9;
const intT ElementTopology::Prism::nFace = 5;
const intT ElementTopology::Prism::Nodes[6] = {0,1,2,3,4,5};
const intT ElementTopology::Prism::Edges[9][2] = {{0,1},{1,2},{2,0},
                                                  {0,3},{1,4},{2,5},
                                                  {3,4},{4,5},{5,3}};
const StaticList2D<intT, 5, 18> ElementTopology::Prism::Faces({4,4,4,3,3},
                                                              {0,1,4,3,
                                                               1,2,5,3,
                                                               0,3,5,2,
                                                               0,2,1,
                                                               3,4,5});
