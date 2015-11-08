/*
 * NodalField.test.cpp
 *
 *  Created on: Oct 19, 2015
 *      Author: rabbit
 */
#include "gtest/gtest.h"
#include "Solution/CGElementField.h"
#include "Mesh/CGMesh.h"
#include "IO/UnstMeshReaderGMSH.h"
TEST(CGElementField, Access)
{
  UnstMeshReaderGMSH mesh_reader("Box.msh-bin", "Box.idmap");
  CGMesh mesh(mesh_reader);
  CGElementField sol_field(mesh, 1);
  intT nnode = mesh.get_MeshGeom().get_nnode();

  for(intT n = 0; n < nnode; n++){
      sol_field(n,0) = (realT)n;
  }
  //std::cout << sol_field.get_Data();
  //SystemModule::pause();

  EXPECT_EQ(nnode, sol_field.get_Data().get_lead_size());

  for(intT n = 0; n < nnode; n++){
    EXPECT_DOUBLE_EQ((realT)n, sol_field(n,0));
  }

  const List2D<intT>& elem2node = mesh.get_MeshElements().get_element2node();

  Array2D<realT> e_sol(4,1);

  for(intT e = 0; e < elem2node.get_lead_size(); e++){
    sol_field.ElementData(e, e_sol);
    for(intT i = 0; i < elem2node.get_ncol(e); i++){
        intT node = elem2node(e,i);
        EXPECT_DOUBLE_EQ((realT)node, e_sol(i,0));
    }
  }

}

TEST(CGElementField, Access2)
{
  UnstMeshReaderGMSH mesh_reader("Box.msh-bin", "Box.idmap");
  CGMesh mesh(mesh_reader);
  intT nnode = mesh.get_MeshGeom().get_nnode();
  CGElementField sol_field(mesh, 5);


  for(intT n = 0; n < nnode; n++){
    for(intT f = 0; f < sol_field.get_nVar(n); f++) {
      sol_field(n,f) = (realT)(n*f);
    }
  }

  for(intT n = 0; n < nnode; n++){
    for(intT f = 0; f < sol_field.get_nVar(n); f++) {
      EXPECT_DOUBLE_EQ((realT)(n*f), sol_field(n,f));
    }
  }

  const List2D<intT>& elem2node = mesh.get_MeshElements().get_element2node();

  Array2D<realT> e_sol(4,5);

  for(intT e = 0; e < elem2node.get_lead_size(); e++){
    sol_field.ElementData(e, e_sol);
    for(intT i = 0; i < elem2node.get_ncol(e); i++){
      intT node = elem2node(e,i);
      for(intT f = 0; f < sol_field.get_nVar(node); f++){
        EXPECT_DOUBLE_EQ((realT)(node*f), e_sol(i,f));
      }
    }
  }
 //---> Pointer check
 EXPECT_DOUBLE_EQ(sol_field.DataPtr(10)[0], 0.0);
 EXPECT_DOUBLE_EQ(sol_field.DataPtr(10)[1], 10.0);
}
