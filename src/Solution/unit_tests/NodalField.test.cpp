/*
 * NodalField.test.cpp
 *
 *  Created on: Oct 19, 2015
 *      Author: rabbit
 */
#include "gtest/gtest.h"
#include "Solution/NodalField.h"
#include "Mesh/CGMesh.h"
#include "IO/UnstMeshReaderGMSH.h"
TEST(NodalField, Access)
{
  UnstMeshReaderGMSH mesh_reader("Box.msh-bin", "Box.idmap");
  CGMesh mesh(mesh_reader);
  intT nnode = mesh.get_MeshGeom().get_nnode();

  Array1D<intT> nvar(nnode);

  nvar.set_value(1);

  NodalField sol_field(mesh, nvar);

  for(intT n = 0; n < nnode; n++){
      sol_field(n,0) = (realT)n;
  }

  EXPECT_EQ(nnode, sol_field.get_Data().get_size(0));
  EXPECT_EQ(nnode+1, sol_field.get_DataIndex().get_size(0));

  for(intT n = 0; n < nnode; n++){
    EXPECT_DOUBLE_EQ((realT)n, sol_field(n,0));
  }

  const List2D<intT>& elem2node = mesh.get_MeshElements().get_element2node();

  Array2D<realT> e_sol(4,2);

  for(intT e = 0; e < elem2node.get_lead_size(); e++){
    sol_field.ElementData(e, e_sol);
    for(intT i = 0; i < elem2node.get_ncol(e); i++){
        intT node = elem2node(e,i);
        EXPECT_DOUBLE_EQ((realT)node, e_sol(i,0));
    }
  }

}


