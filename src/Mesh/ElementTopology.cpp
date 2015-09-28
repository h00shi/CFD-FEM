#include "Mesh/ElementTopology.h"
//****************************************************************************80
const intT ElementTopology::Bar::Type = ElementTopology::element_types::BAR;
const intT ElementTopology::Bar::nNode = 2;
const intT ElementTopology::Bar::nEdge = 1;
const intT ElementTopology::Bar::nFace = 2;
const intT ElementTopology::Bar::Nodes[2] = {0,1};
const intT ElementTopology::Bar::Edges[1][2] = {{0,1}};
const StaticList2D<intT, 2, 2> ElementTopology::Bar::Faces({1,1}, {0,1});
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
