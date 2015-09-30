#include "Mesh/ElementTopology.h"
//****************************************************************************80
const intT ElementTopology::Triangle::Type = 
  ElementTopology::element_types::TRI;
const intT ElementTopology::Triangle::nNode = 3;
const intT ElementTopology::Triangle::nEdge = 3;
const intT ElementTopology::Triangle::nFace = 3;
const intT ElementTopology::Triangle::Nodes[3] = {0,1,2};
const intT ElementTopology::Triangle::Edges[3][2] = {{0,1},{1,2},{2,0}};
const StaticList2D<intT, 3, 6> ElementTopology::Triangle::Faces({2,2,2}, 
                                                                {0,1,
                                                                 1,2,
                                                                 2,0});
