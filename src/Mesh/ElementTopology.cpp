#include "Mesh/ElementTopology.h"
//****************************************************************************80
const intT ElementTopology::Bar::Type = ElementTopology::element_types::BAR;
const intT ElementTopology::Bar::nNode = 2;
const intT ElementTopology::Bar::nEdge = 1;
const intT ElementTopology::Bar::nFace = 2;
const intT ElementTopology::Bar::Nodes[2] = {0,1};
const intT ElementTopology::Bar::Edges[1][2] = {{0,1}};
const StaticList2D<intT,2,2> ElementTopology::Bar::Faces({1,1}, {0,1});
//****************************************************************************80
const 
