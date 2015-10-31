#include "gtest/gtest.h"
#include "SparseMatrix/CSRMatrix.h"
#include "Mesh/CGMesh.h"
#include "DataStructures/Graph.h"
#include "IO/UnstMeshReaderNKBGrid.h"
TEST(CSRMatrix, OneField){
  UnstMeshReaderNKBGrid mesh_reader("Square.grid");
  CGMesh mesh(mesh_reader);
  Array1D<intT> nfld(mesh.get_MeshGeom().get_nnode());
  nfld.set_value(1);
  CSRMatrix<realT> csr_matrix(mesh.get_Graph(), nfld);

  const List2D<intT>& adj = mesh.get_Graph().get_GraphAdj();
  const List2D<intT>& adj_data_offset = csr_matrix.get_adj_data_offset();
  const Array1D<intT>& rowos = csr_matrix.get_row_offset();
  const List2D<intT>& colidx = csr_matrix.get_column_idx();
  const Array1D<intT>& nnz_node = csr_matrix.get_nnz_node();
  //---> Check rowos
  EXPECT_EQ(rowos(0), 0);
  for(intT i = 0; i < adj.get_lead_size(); i++){
    EXPECT_EQ(rowos(i+1), rowos(i) + adj.get_ncol(i));
   
    for(intT j = 0; j < adj.get_ncol(i); j++){
      EXPECT_EQ(colidx(i,j),adj(i,j));
      EXPECT_EQ(rowos(i) + j, adj_data_offset(i,j));
    }
  }

  //---> Check nnz_node
  for(intT i = 0; i < adj.get_lead_size(); i++){
    intT nnz = 0;
    for(intT j = 0; j < adj.get_ncol(i); j++){
      nnz += nfld(adj(i,j));
    }
    EXPECT_EQ(nnz, nnz_node(i));
  }
}

TEST(CSRMatrix, TwoFields) {
  UnstMeshReaderNKBGrid mesh_reader("Square.grid");
  CGMesh mesh(mesh_reader);
  Array1D<intT> nfld(mesh.get_MeshGeom().get_nnode());
  nfld.set_value(2);
  CSRMatrix<realT> csr_matrix(mesh.get_Graph(), nfld);
  const List2D<intT>& adj = mesh.get_Graph().get_GraphAdj();
  const List2D<intT>& adj_data_offset = csr_matrix.get_adj_data_offset();
  const Array1D<intT>& rowos = csr_matrix.get_row_offset();
  const List2D<intT>& colidx = csr_matrix.get_column_idx();
  const Array1D<intT>& nnz_node = csr_matrix.get_nnz_node();
  
  Array1D<intT> ncol_per_node(adj.get_lead_size());
  
  for(intT i = 0; i < adj.get_lead_size(); i++){
    for(intT j = 0; j < adj.get_ncol(i); j++){
      ncol_per_node(i) += nfld(adj(i,j));
    }
  }
  
  //---> Check Rowos
  intT row = 0;
  for(intT i = 0; i < adj.get_lead_size(); i++){
    for(intT r = 0; r < nfld(i); r++){
      EXPECT_EQ(rowos(row + r + 1), rowos(row + r) + ncol_per_node(i) );
      for(intT j = 0; j < adj.get_ncol(i); j++){
	for(intT c = 0; c < nfld(adj(i,j)); c++){
	  EXPECT_EQ(colidx(row + r, j*2 + c), 2*(adj(i,j)) + c);
	}
      }
    }
    row += nfld(i);
  }
  
  //---> Check data offset
  intT index = 0;
  for(intT i = 0; i < adj.get_lead_size(); i++){
    intT loc_index = index;
    for(intT j = 0; j < adj.get_ncol(i); j++){
      EXPECT_EQ(loc_index, adj_data_offset(i,j));  
            
      loc_index += nfld(adj(i,j));
    }
    index += ncol_per_node(i)*nfld(i) ;
  
  }
  
  //---> Check nnz_node
  for(intT i = 0; i < adj.get_lead_size(); i++){
    intT nnz = 0;
    for(intT j = 0; j < adj.get_ncol(i); j++){
      nnz += nfld(adj(i,j));
    }
    EXPECT_EQ(nnz, nnz_node(i));
  }
}

TEST(CSRMatrix, ThreeFields) {
  UnstMeshReaderNKBGrid mesh_reader("Square.grid");
  CGMesh mesh(mesh_reader);
  Array1D<intT> nfld(mesh.get_MeshGeom().get_nnode());
  nfld.set_value(3);
  CSRMatrix<realT> csr_matrix(mesh.get_Graph(), nfld);
  const List2D<intT>& adj = mesh.get_Graph().get_GraphAdj();
  const List2D<intT>& adj_data_offset = csr_matrix.get_adj_data_offset();
  const Array1D<intT>& rowos = csr_matrix.get_row_offset();
  const List2D<intT>& colidx = csr_matrix.get_column_idx();
  const Array1D<intT>& nnz_node = csr_matrix.get_nnz_node();
  
  Array1D<intT> ncol_per_node(adj.get_lead_size());
  
  for(intT i = 0; i < adj.get_lead_size(); i++){
    for(intT j = 0; j < adj.get_ncol(i); j++){
      ncol_per_node(i) += nfld(adj(i,j));
    }
  }
  
  //---> Check Rowos
  intT row = 0;
  for(intT i = 0; i < adj.get_lead_size(); i++){
    for(intT r = 0; r < nfld(i); r++){
      EXPECT_EQ(rowos(row + r + 1), rowos(row + r) + ncol_per_node(i) );
      for(intT j = 0; j < adj.get_ncol(i); j++){
	for(intT c = 0; c < nfld(adj(i,j)); c++){
	  EXPECT_EQ(colidx(row + r, j*3 + c), 3*(adj(i,j)) + c);
	}
      }
    }
    row += nfld(i);
  }
  
  //---> Check data offset
  intT index = 0;
  for(intT i = 0; i < adj.get_lead_size(); i++){
    intT loc_index = index;
    for(intT j = 0; j < adj.get_ncol(i); j++){
      EXPECT_EQ(loc_index, adj_data_offset(i,j));  
      
      loc_index += nfld(adj(i,j));
    }
    index += ncol_per_node(i)*nfld(i) ;
  
  }
  //---> Check nnz_node
  for(intT i = 0; i < adj.get_lead_size(); i++){
    intT nnz = 0;
    for(intT j = 0; j < adj.get_ncol(i); j++){
      nnz += nfld(adj(i,j));
    }
    EXPECT_EQ(nnz, nnz_node(i));
  }

  //---> Check row
  row = 0;
  for(intT i = 0; i < mesh.get_MeshGeom().get_nnode(); i++)
  {
    for(intT j = 0; j < nfld(i); j++){
      EXPECT_EQ(row, csr_matrix.Row(i,j));
      row ++;
    }
  }

  csr_matrix.Diagnostic(std::cout);
}

TEST(CSRMatrix, Tridiag){
  const intT N = 20;
  //---> Setup the adjacency and 
  List2D<intT> adj(N, (N-2)*3 + 2 + 2);
  Array1D<intT> ncol(N);
  Array2D<intT> edge2node(N-1,2);
  Array1D<intT> nfld(N);
  
  ncol(0) = 2;
  for(intT i = 1; i < N-1; i++){
    ncol(i) = 3;
  }
  ncol(N - 1) = 2;
  adj.set_ncol(ncol);
 
  adj(0,0) = 0;
  adj(0,1) = 1;
  for(intT i = 1; i < N-1; i++){
    adj(i,0) = i - 1;
    adj(i,1) = i;
    adj(i,2) = i + 1;
  }

  adj(N-1,0) = N-2;
  adj(N-1,1) = N-1;
  
  for(intT i = 0; i < N-1; i++){
    edge2node(i,0) = i;
    edge2node(i,1) = i + 1;
  }
  nfld.set_value(1);
  Graph graph(adj, edge2node);
  
  CSRMatrix<realT> matrix(graph, nfld);
  matrix(0,0,0,0) = -2.0;
  matrix(0,1,0,0) = 1.0;
  for(intT i = 1; i < N - 1; i++){
    matrix(i,0,0,0) = -1.0;
    matrix(i,1,0,0) = -2.0;
    matrix(i,2,0,0) = 1.0;
  }
  matrix(N-1,0,0,0) = -1.0;
  matrix(N-1,1,0,0) = -2.0;
 
  for(intT i = 1; i < matrix.get_nrow()-1; i++){
    EXPECT_DOUBLE_EQ(-1.0, matrix(i,0,0,0));
    EXPECT_DOUBLE_EQ(-2.0, matrix(i,1,0,0));
    EXPECT_DOUBLE_EQ(1.0, matrix(i,2,0,0));
  }

  //---> Check ro
  for(intT i = 0; i < N; i++)
  {
    EXPECT_EQ(i, matrix.Row(i,0));
  }

}
