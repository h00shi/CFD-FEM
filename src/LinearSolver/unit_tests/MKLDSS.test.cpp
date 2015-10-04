#include "gtest/gtest.h"
#include "SparseMatrix/CSRMatrix.h"
#include "DataStructures/Graph.h"
#include "LinearSolver/MKLDSS.h"

TEST(MKLDSS, Tridiag){
  const intT N = 10;
  //---> Setup the adjacency and 
  List2D<intT> adj(N, (N-2)*3 + 2 + 2);
  Array1D<intT> ncol(N);
  Array2D<intT> edge2node(N-1,2);
  Array1D<intT> nfld(N);
  
  nfld.set_value(1);
    
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

  Graph graph(adj);
  CSRMatrix<realT> matrix(graph, nfld);
  
  List2D<realT> b(N,N);
  List2D<realT> x(N,N);
  b.set_ncol(nfld);
  x.set_ncol(nfld);

  matrix(0,0,0,0) = -2.0;
  matrix(0,1,0,0) = 1.0;
  b(0,0) = -2.0;
  for(intT i = 1; i < N - 1; i++){
    matrix(i,0,0,0) = 1.0;
    matrix(i,1,0,0) = -2.0;
    matrix(i,2,0,0) = 1.0;
    b(i) = 0.0;
  }
  matrix(N-1,0,0,0) = 1.0;
  matrix(N-1,1,0,0) = -2.0;
  b(N-1,0) = -1.0;

  MKLDSS mkl_dss_solver(matrix);
  mkl_dss_solver.Factorize();
  mkl_dss_solver.Solve(b,x);
  
  Array1D<double> matlab_ans(10);

  matlab_ans(0) = 1.909090909090909;
  matlab_ans(1) = 1.818181818181818;
  matlab_ans(2) = 1.727272727272727;
  matlab_ans(3) = 1.636363636363636;
  matlab_ans(4) = 1.545454545454546;
  matlab_ans(5) = 1.454545454545454;
  matlab_ans(6) = 1.363636363636363;
  matlab_ans(7) = 1.272727272727272;
  matlab_ans(8) = 1.181818181818182;
  matlab_ans(9) = 1.090909090909091; 

  for(int e = 0; e < N; e++){
    EXPECT_DOUBLE_EQ(matlab_ans(e),x(e,0));
  }
}
TEST(MKLDSS, Tridiag_2_Field){
  const intT N = 10;
  //---> Setup the adjacency and 
  List2D<intT> adj(N, (N-2)*3 + 2 + 2);
  Array1D<intT> ncol(N);
  Array2D<intT> edge2node(N-1,2);
  Array1D<intT> nfld(N);
  
  nfld.set_value(2);
    
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

  Graph graph(adj);
  CSRMatrix<realT> matrix(graph, nfld);
  
  List2D<realT> b(N,N*2);
  List2D<realT> x(N,N*2);
  b.set_ncol(nfld);
  x.set_ncol(nfld);

  matrix(0,0,0,0) = -2.0;
  matrix(0,0,1,1) = -4.0;
  matrix(0,1,0,0) = 1.0;
  matrix(0,1,1,1) = 2.0;

  b(0,0) = -2.0;
  b(0,1) = -10.0;
  for(intT i = 1; i < N - 1; i++){
    matrix(i,0,0,0) = 1.0;
    matrix(i,0,1,1) = 2.0;
    matrix(i,1,0,0) = -2.0;
    matrix(i,1,1,1) = -4.0;
    matrix(i,2,0,0) = 1.0;
    matrix(i,2,1,1) = 2.0;
    b(i,0) = 0.0;
    b(i,1) = 0.0;
  }
  matrix(N-1,0,0,0) = 1.0;
  matrix(N-1,0,1,1) = 2.0;
  matrix(N-1,1,0,0) = -2.0;
  matrix(N-1,1,1,1) = -4.0;

  b(N-1,0) = -1.0;
  b(N-1,1) = -2.0;

  MKLDSS mkl_dss_solver(matrix);
  mkl_dss_solver.Factorize();
  mkl_dss_solver.Solve(b,x);
  
  Array1D<double> matlab_ans(10);

  matlab_ans(0) = 1.909090909090909;
  matlab_ans(1) = 1.818181818181818;
  matlab_ans(2) = 1.727272727272727;
  matlab_ans(3) = 1.636363636363636;
  matlab_ans(4) = 1.545454545454546;
  matlab_ans(5) = 1.454545454545454;
  matlab_ans(6) = 1.363636363636363;
  matlab_ans(7) = 1.272727272727272;
  matlab_ans(8) = 1.181818181818182;
  matlab_ans(9) = 1.090909090909091; 
  
  Array1D<double> matlab_ans2(10);
  matlab_ans2(0) = 4.636363636363636;
  matlab_ans2(1) = 4.272727272727272;
  matlab_ans2(2) = 3.909090909090909;
  matlab_ans2(3) = 3.545454545454545;
  matlab_ans2(4) = 3.181818181818182;
  matlab_ans2(5) = 2.818181818181818;
  matlab_ans2(6) = 2.454545454545454;
  matlab_ans2(7) = 2.090909090909091;
  matlab_ans2(8) = 1.727272727272727;
  matlab_ans2(9) = 1.363636363636364;
  for(int e = 0; e < N; e++){
    EXPECT_DOUBLE_EQ(matlab_ans(e),x(e,0));
    EXPECT_DOUBLE_EQ(matlab_ans2(e),x(e,1));
  }
}
