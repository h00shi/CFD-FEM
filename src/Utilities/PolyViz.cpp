//****************************************************************************80
//! \file poly_viz.cpp
//! \brief  A driver to write tecplot viz files for the basis functions
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! 
//****************************************************************************80
#include "my_incl.h"
#include "Elements/BarElement.h"
#include "Elements/TriElement.h"
#include "Elements/TetElement.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/Array3D.h"
#include <fstream>
int main(int argc, char** argv)
{
  int p = 5;
  //std::cout << "Enter the polynomial order you want to see: ";
  //std::cin >> p;
  
  int ndof1D = p + 1;
  int ndof2D = (p + 1)*(p + 2)/2;
  int ndof3D = (p + 1)*(p + 2)*(p + 3)/6;
  

  //---> Setup coordinates
  const int np = 30;
  const double dx = 2.0/((double)(np-1));

  Array2D<double> xi1D(np, 1);
  Array3D<double> xi2D(np, np, 2);
  Array4D<double> xi3D(np, np, np, 3);

  for(int i = 0; i < np; i++) {
    xi1D(i,0) = -1.0 + i*dx;
  }
  
  for(int i = 0; i < np; i++){
    for(int j = 0; j < np; j++){
      double xi = -1.0 + dx*i;
      double eta = -1.0 + dx*j;
      
      xi2D(i,j,0) = xi;
      xi2D(i,j,1) = (1.0 + eta)*(1 - xi)*half - 1.0;
    }
  }

  for(int i = 0; i < np; i++){
    for(int j = 0; j < np; j++){
      for(int k =0; k < np; k++){
  	double xi = -1.0 + dx*i;
  	double eta = -1.0 + dx*j;
  	double zeta = -1.0 + dx*k;
	
  	xi3D(i,j,k,0) = (1.0 + (1.0 + xi)*(1.0 - zeta)*half - 1.0) * 
  	  (1.0 - eta)*half - 1.0;
  	xi3D(i,j,k,1) = (1.0 + eta)*(1 - zeta)*half - 1.0;
  	xi3D(i,j,k,2) = zeta;
      }
    }
  }



  //---> Intiate elements
  BarElement bar;
  TriElement tri;
  TetElement tet;

  //---> Write 1-D data
  std::ofstream tecfile;
  tecfile.open("1_D_Basis.tec");
  tecfile << "Title = \"1-D Basis Functions\"" << std::endl;
  tecfile << "Variables = Xi Phi" << std::endl;

  for(int dof = 0 ; dof < ndof1D; dof++) {
    tecfile << "ZONE T= \"PHI-" << dof << "\"" << std::endl;
    tecfile << "ZONETYPE=ORDERED" << std::endl;
    tecfile << "I=" << np << ", J=" << 1 << ", K=" << 1 << std::endl;
    tecfile << "DATAPACKING=POINT" << std::endl;
    for(int i = 0; i < np; i++) {
      tecfile << xi1D(i,0) << " " << bar.EvalBasis(dof,xi1D(i,0))  
	      << std::endl;
    }
    
  } // End 1-D Dof loop 

  tecfile.close();

  //---> Write 2-D data
  tecfile.open("2_D_Basis.tec");
  tecfile << "Title = \"2-D Basis Functions\"" << std::endl;
  tecfile << "Variables = Xi Eta Phi" << std::endl;

  for(int dof = 0 ; dof < ndof2D; dof++) {
    tecfile << "ZONE T= \"PHI-" << dof << "\"" << std::endl;
    tecfile << "ZONETYPE=ORDERED" << std::endl;
    tecfile << "I=" << np << ", J=" << np << ", K=" << 1 << std::endl;
    tecfile << "DATAPACKING=POINT" << std::endl;
    for(int i = 0; i < np; i++) {
      for(int j =0; j < np; j++){
	tecfile << xi2D(i,j,0) << " " << xi2D(i,j,1) << " " 
		<< tri.EvalBasis(dof, xi2D(i,j,0), xi2D(i,j,1) )  
		<< std::endl;
      }
    }
  } // End 2-D Dof loop 

  tecfile.close();
  //---> Write 3-D data
  tecfile.open("3_D_Basis.tec");
  tecfile << "Title = \"3-D Basis Functions\"" << std::endl;
  tecfile << "Variables = Xi Eta Zeta Phi" << std::endl;

  for(int dof = 0 ; dof < ndof3D; dof++) {
    tecfile << "ZONE T= \"PHI-" << dof << "\"" << std::endl;
    tecfile << "ZONETYPE=ORDERED" << std::endl;
    tecfile << "I=" << np << ", J=" << np << ", K=" << np << std::endl;
    tecfile << "DATAPACKING=POINT" << std::endl;
    for(int i = 0; i < np; i++) {
      for(int j =0; j < np; j++){
	for(int k =0; k < np; k++){
	  tecfile << xi3D(i, j, k, 0) << " " << xi3D(i, j, k, 1) << " " 
		  << xi3D(i,j, k, 2) << " "  
		  << tet.EvalBasis(dof, xi3D(i,j,k,0), xi3D(i,j,k,1), 
				    xi3D(i,j,k,2) )  
		  << std::endl;
	}
      }
    }
    
  } // End 3-D Dof loop 

  tecfile.close();
  

  
}
