#include "gtest/gtest.h"
#include "DataStructures/CSRMatrix.h"
#include "Mesh/UnstMesh.h"

TEST(CSRMatrix, CRTP) {
  UnstMesh mesh("Square.grid", "Grid-NKB");
  Array1D<intT> nfld(mesh.get_nnode());
  nfld.set_value(2);

  CSRMatrix<realT> csr_matrix(mesh.get_adj(), mesh.get_edge2node(), nfld);
  
  SparseMatrix< realT, CSRMatrix<realT> >& matrix = csr_matrix;
  matrix(0,0,0,0) = 1.0;

}


