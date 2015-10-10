#include "Mesh/ElementTopology.h"
//****************************************************************************80
const intT ElementTopology::Tetrahedron::Type = 
  ElementTopology::element_types::TET;
const intT ElementTopology::Tetrahedron::nNode = 4;
const intT ElementTopology::Tetrahedron::nEdge = 6;
const intT ElementTopology::Tetrahedron::nFace = 4;
const intT ElementTopology::Tetrahedron::Nodes[4] = {0,1,2,3};
const intT ElementTopology::Tetrahedron::Edges[6][2] = {{0,1},{1,2},{2,0},{0,3},
                                                {1,3},{2,3}};
const StaticList2D<intT, 4, 12> ElementTopology::Tetrahedron::Faces({3,3,3,3},
                                    {0,1,3,
                                     1,2,3,
                                     2,0,3,
                                     0,2,1});
