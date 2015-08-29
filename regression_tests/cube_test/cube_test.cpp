#include "gtest/gtest.h"
#include "my_incl.h"
#include "UnstGrid.h" /* Definition of class UnstGrid */
#include "InputIO.h" /* Definition of class InputIO */
#include "Solution.h"
#include "Array1D.h"
#include "Equation.h"
#include "consts.h"
#include "MeshIO.h"
#include "GLSResidual.h"
//---> Namespace includes: These are header files that define namespaces
#include "system_module.h" /* Includes namespace named:system_module */
#include "setup.h" /* Includes namespace named:setup */
#include "SolGlobalVars.h"


//****************************************************************************80
//! \file cube_test.cpp
//! \brief  Takes in a cube mesh in FV-UNS format and runs several tests to 
//!         ensure that freestream in preserved.  
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! 
//****************************************************************************80
TEST(CubeTest, cond1) {
  double Fmach = .25;
  double Alpha = 0.0;
  double Beta = 0.0;

  //---> Instatiate mesh file reader (MESHIO)
  MeshIO<int,double> mesh_reader; 
  
  //---> Initalize Solver Global Vars
  double uinf = Fmach*cos(-pi/180.0*Alpha)*cos(pi/180.0*Beta);
  double vinf = Fmach*sin(pi/180.0*Beta);
  double winf = Fmach*sin(pi/180.0*Alpha)*cos(pi/180.0*Beta);
  
  //--------------------------- Instiate and Setup Grid ------------------------
  //---> Instatiate Solution and Grid
  UnstGrid<int,double> grid(3);
  mesh_reader.read_fvuns_ascii("test", grid); 
  mesh_reader.read_bc("test", grid);
  
  //---> Now setup the rest of the connectivity
  grid.init_connectivity();
  
  EulerPDESet3D<int,double> Euler(1.0, uinf, vinf, winf, 1.0/1.4);
  Array1D<double> qinf(Euler.get_nfld());
  Euler.get_qinf(qinf);
    
  Solution<int,double> soln(grid, Euler.get_nfld(), 0, qinf, 1);
  
  GLSResidual< int,double,EulerPDESet3D<int,double> > gls_resid(Euler, grid, 
								 soln);
  Array2D<double> resid(grid.nnode, Euler.get_nfld()); 
  gls_resid.GLSResid(resid);

  double l2n = gls_resid.l2norm(resid);
  EXPECT_PRED_FORMAT2(testing::DoubleLE, l2n, 1.0e-14);
  
}

TEST(CubeTest, cond2) {
  double Fmach = .5;
  double Alpha = 10.0;
  double Beta = 25.0;

  //---> Instatiate mesh file reader (MESHIO)
  MeshIO<int,double> mesh_reader; 
  
  //---> Initalize Solver Global Vars
  double uinf = Fmach*cos(-pi/180.0*Alpha)*cos(pi/180.0*Beta);
  double vinf = Fmach*sin(pi/180.0*Beta);
  double winf = Fmach*sin(pi/180.0*Alpha)*cos(pi/180.0*Beta);
  
  //--------------------------- Instiate and Setup Grid ------------------------
  //---> Instatiate Solution and Grid
  UnstGrid<int,double> grid(3);
  mesh_reader.read_fvuns_ascii("test", grid); 
  mesh_reader.read_bc("test", grid);
  
  //---> Now setup the rest of the connectivity
  grid.init_connectivity();
  
  EulerPDESet3D<int,double> Euler(1.0, uinf, vinf, winf, 1.0/1.4);
  Array1D<double> qinf(Euler.get_nfld());
  Euler.get_qinf(qinf);
    
  Solution<int,double> soln(grid, Euler.get_nfld(), 0, qinf, 1);
  
  GLSResidual< int,double,EulerPDESet3D<int,double> > gls_resid(Euler, grid, 
								 soln);
  Array2D<double> resid(grid.nnode, Euler.get_nfld()); 
  gls_resid.GLSResid(resid);

  double l2n = gls_resid.l2norm(resid);
  EXPECT_PRED_FORMAT2(testing::DoubleLE, l2n, 1.0e-14);
  
}

TEST(CubeTest, cond3) {
  double Fmach = 1.5;
  double Alpha = 100.0;
  double Beta = -25.0;

  //---> Instatiate mesh file reader (MESHIO)
  MeshIO<int,double> mesh_reader; 
  
  //---> Initalize Solver Global Vars
  double uinf = Fmach*cos(-pi/180.0*Alpha)*cos(pi/180.0*Beta);
  double vinf = Fmach*sin(pi/180.0*Beta);
  double winf = Fmach*sin(pi/180.0*Alpha)*cos(pi/180.0*Beta);
  
  //--------------------------- Instiate and Setup Grid ------------------------
  //---> Instatiate Solution and Grid
  UnstGrid<int,double> grid(3);
  mesh_reader.read_fvuns_ascii("test", grid); 
  mesh_reader.read_bc("test", grid);
  
  //---> Now setup the rest of the connectivity
  grid.init_connectivity();
  
  EulerPDESet3D<int,double> Euler(1.0, uinf, vinf, winf, 1.0/1.4);
  Array1D<double> qinf(Euler.get_nfld());
  Euler.get_qinf(qinf);
    
  Solution<int,double> soln(grid, Euler.get_nfld(), 0, qinf, 1);
  
  GLSResidual< int,double,EulerPDESet3D<int,double> > gls_resid(Euler, grid, 
								 soln);
  Array2D<double> resid(grid.nnode, Euler.get_nfld()); 
  gls_resid.GLSResid(resid);

  double l2n = gls_resid.l2norm(resid);
  EXPECT_PRED_FORMAT2(testing::DoubleLE, l2n, 1.0e-14);
  
}
