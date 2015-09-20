#include "gtest/gtest.h"
#include "SparseMatrix/CSRMatrix.h"
#include "Mesh/UnstMesh.h"


TEST(CSRMatrix, OneField){
  UnstMesh mesh("Square.grid", "Grid-NKB");
  Array1D<intT> nfld(mesh.get_nnode());
  nfld.set_value(1);
  CSRMatrix<realT> csr_matrix(mesh.get_adj(), mesh.get_edge2node(), nfld);
  const List2D<intT>& adj = mesh.get_adj();
  const List2D<intT>& adj_data_offset = csr_matrix.get_adj_data_offset();
  const Array1D<intT>& rowos = csr_matrix.get_row_offset();
  const List2D<intT>& colidx = csr_matrix.get_column_idx();
  const Array1D<intT>& node2diag = csr_matrix.get_self_adj_index();
  std::cout << node2diag << std::endl;
  //---> Check rowos
  EXPECT_EQ(rowos(0), 0);
  for(intT i = 0; i < adj.get_lead_size(); i++){
    EXPECT_EQ(rowos(i+1), rowos(i) + adj.get_ncol(i));
   
    for(intT j = 0; j < adj.get_ncol(i); j++){
      EXPECT_EQ(colidx(i,j),adj(i,j));
      EXPECT_EQ(rowos(i) + j, adj_data_offset(i,j));
      if( adj(i,j) == i ){
	EXPECT_EQ(j, node2diag(i));
      }
    }
  }

}

TEST(CSRMatrix, TwoFields) {
  UnstMesh mesh("Square.grid", "Grid-NKB");
  Array1D<intT> nfld(mesh.get_nnode());
  nfld.set_value(2);

  CSRMatrix<realT> csr_matrix(mesh.get_adj(), mesh.get_edge2node(), nfld);
  const List2D<intT>& adj = mesh.get_adj();
  const List2D<intT>& adj_data_offset = csr_matrix.get_adj_data_offset();
  const Array1D<intT>& rowos = csr_matrix.get_row_offset();
  const List2D<intT>& colidx = csr_matrix.get_column_idx();
  const Array1D<intT>& node2diag = csr_matrix.get_self_adj_index();
  
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
      if( i == adj(i,j) ) {EXPECT_EQ(j, node2diag(i));}
      
      loc_index += nfld(adj(i,j));
    }
    index += ncol_per_node(i)*nfld(i) ;
  
  }
  
}

TEST(CSRMatrix, ThreeFields) {
  UnstMesh mesh("Square.grid", "Grid-NKB");
  Array1D<intT> nfld(mesh.get_nnode());
  nfld.set_value(3);

  CSRMatrix<realT> csr_matrix(mesh.get_adj(), mesh.get_edge2node(), nfld);
  const List2D<intT>& adj = mesh.get_adj();
  const List2D<intT>& adj_data_offset = csr_matrix.get_adj_data_offset();
  const Array1D<intT>& rowos = csr_matrix.get_row_offset();
  const List2D<intT>& colidx = csr_matrix.get_column_idx();
  const Array1D<intT>& node2diag = csr_matrix.get_self_adj_index();
  
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
      if( i == adj(i,j) ) {EXPECT_EQ(j, node2diag(i));}
      
      loc_index += nfld(adj(i,j));
    }
    index += ncol_per_node(i)*nfld(i) ;
  
  }
  std::cout << mesh.get_edge2node().get_size(0) << std::endl;
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
  
  CSRMatrix<realT> matrix(adj, edge2node, nfld);
  matrix(0,0,0,0) = -2.0;
  matrix(0,1,0,0) = 1.0;
  for(intT i = 1; i < adj.get_lead_size() - 1; i++){
    matrix(i,0,0,0) = -1.0;
    matrix(i,1,0,0) = -2.0;
    matrix(i,2,0,0) = 1.0;
  }
  matrix(N-1,0,0,0) = -1.0;
  matrix(N-1,1,0,0) = -2.0;
 
  for(intT i = 0; i < matrix.get_nrow(); i++){
    for(intT j = 0; j < matrix.get_column_idx().get_ncol(i); j++){
      std::cout << matrix(i,j,0,0) << " ";
    }
    std::cout << std::endl;
  }

}
