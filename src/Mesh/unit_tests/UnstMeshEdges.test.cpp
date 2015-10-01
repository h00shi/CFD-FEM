#include "Mesh/UnstMesh.h"
#include "Mesh/UnstMeshEdges.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/List2D.h"
#include "gtest/gtest.h"
#include <set>
TEST(UnstMeshEdges, Construct)
{
  UnstMesh mesh("Square.grid", "Grid-NKB");
  UnstMeshEdges edges(mesh.get_nnode(), mesh.get_nelement(),
                      mesh.get_element2node(), mesh.get_node2element(),
                      mesh.get_element_type());
  
  const Array2D<intT>& edge2node = edges.get_edge2node();
  std::cout << edge2node << std::endl;
}

TEST(UnstMeshEdges, STLConstruct)
{
  UnstMesh mesh("Square.grid", "Grid-NKB");
  std::set<std::pair<intT>> edge2node;
  
  for(intT i = 0; i < mesh.get_nelement(); i++){
    for(intT j = 0; j < mesh.get_element2node().get_ncol(e); j++){
     
      std::pair<intT,intT> edge = 
}
