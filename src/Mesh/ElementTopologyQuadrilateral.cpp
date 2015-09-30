#include "Mesh/ElementTopology.h"
//****************************************************************************80
const intT ElementTopology::Quadrilateral::Type = 
  ElementTopology::element_types::QUAD;
const intT ElementTopology::Quadrilateral::nNode = 4;
const intT ElementTopology::Quadrilateral::nEdge = 4;
const intT ElementTopology::Quadrilateral::nFace = 4;
const intT ElementTopology::Quadrilateral::Nodes[4] = {0,1,2,3};
const intT ElementTopology::Quadrilateral::Edges[4][2] = {{0,1},{1,2},{2,3},
                                                          {3,0}};
const StaticList2D<intT, 4, 8> ElementTopology::Quadrilateral::Faces({2,2,2,2}, 
                                                                     {0,1,
                                                                      1,2,
                                                                      2,3,
                                                                      3,0});
