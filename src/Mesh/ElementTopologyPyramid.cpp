#include "Mesh/ElementTopology.h"
//****************************************************************************80
const intT ElementTopology::Pyramid::Type = ElementTopology::element_types::PYR;
const intT ElementTopology::Pyramid::nNode = 5;
const intT ElementTopology::Pyramid::nEdge = 8;
const intT ElementTopology::Pyramid::nFace = 5;
const intT ElementTopology::Pyramid::Nodes[5] = {0,1,2,3,4};
const intT ElementTopology::Pyramid::Edges[8][2] = {{0,1},{1,2},{2,3},{3,0},
                                                    {0,4},{1,4},{2,4},{3,4}};
const StaticList2D<intT, 5, 16> ElementTopology::Pyramid::Faces({4,3,3,3,3}, 
                                                                {0,3,2,1,
                                                                 0,1,4,
                                                                 1,2,4,
                                                                 3,4,2,
                                                                 0,4,3});
