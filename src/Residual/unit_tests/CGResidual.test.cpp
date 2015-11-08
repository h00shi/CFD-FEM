/*
 * CGResidual.test.cpp
 *
 *  Created on: Nov 3, 2015
 *      Author: rabbit
 */
#include "gtest/gtest.h"
#include "Residual/CGResidual.h"
#include "IO/UnstMeshReaderNKBGrid.h"
#include "PDE/Equation.h"
TEST(CGResidual, Construction)
{
  UnstMeshReaderNKBGrid grid_reader("Square.grid");
  CGMesh mesh(grid_reader);
  Poisson poisson_pde;
  Array1D<intT> nvar(mesh.get_MeshGeom().get_nnode());
  nvar.set_value(static_cast<intT>(Poisson::nfld_));

  CGResidual<Poisson> cg_resid(poisson_pde, mesh, 1, 1, 2, nvar);

}
