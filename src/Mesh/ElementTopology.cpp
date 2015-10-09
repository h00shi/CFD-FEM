/*
 * ElementTopology.cpp
 *
 *  Created on: Oct 9, 2015
 *      Author: rabbit
 */

#include "Mesh/ElementTopology.h"


//****************************************************************************80
  template
  intT ElementTopology::FindLocalFace<ElementTopology::Bar>
  (const std::vector<intT>& face_nodes,
   const std::vector<intT>& elem_nodes);
  template
  intT ElementTopology::FindLocalFace<ElementTopology::Triangle>
  (const std::vector<intT>& face_nodes,
   const std::vector<intT>& elem_nodes);
  template
  intT ElementTopology::FindLocalFace<ElementTopology::Quadrilateral>
  (const std::vector<intT>& face_nodes,
   const std::vector<intT>& elem_nodes);
  template
  intT ElementTopology::FindLocalFace<ElementTopology::Tetrahedron>
  (const std::vector<intT>& face_nodes,
   const std::vector<intT>& elem_nodes);
  template
  intT ElementTopology::FindLocalFace<ElementTopology::Prism>
  (const std::vector<intT>& face_nodes,
   const std::vector<intT>& elem_nodes);
  template
  intT ElementTopology::FindLocalFace<ElementTopology::Pyramid>
  (const std::vector<intT>& face_nodes,
   const std::vector<intT>& elem_nodes);
  template
  intT ElementTopology::FindLocalFace<ElementTopology::Hexahedron>
  (const std::vector<intT>& face_nodes,
   const std::vector<intT>& elem_nodes);
