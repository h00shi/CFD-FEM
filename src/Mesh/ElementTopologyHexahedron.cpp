#include "Mesh/ElementTopology.h"
//****************************************************************************80
const intT ElementTopology::Hexahedron::Type = 
  ElementTopology::element_types::HEX;
const intT ElementTopology::Hexahedron::nNode = 8;
const intT ElementTopology::Hexahedron::nEdge = 12;
const intT ElementTopology::Hexahedron::nFace = 6;
const intT ElementTopology::Hexahedron::Nodes[8] = {0,1,2,3,4,5,6,7};
const intT ElementTopology::Hexahedron::Edges[12][2] = {{0,1},{1,2},{2,3},{3,0},
                                                        {0,4},{1,5},{2,6},{3,7},
                                                        {4,5},{5,6},{6,7},
                                                        {7,4}};
const StaticList2D<intT, 6, 24> 
ElementTopology::Hexahedron::Faces({4,4,4,4,4,4}, 
                                   {3,0,4,7,
                                    1,2,6,5,
                                    0,1,5,4,
                                    2,3,7,6,
                                    3,2,1,0,
                                    4,5,6,7});
