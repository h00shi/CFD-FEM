#include "my_incl.h"
#include "UnstGrid.h" /* Definition of class UnstGrid */
#include "InputIO.h" /* Definition of class InputIO */
#include "Solution.h"
#include "Array1D.h"
#include "Equation.h"
#include "consts.h"
#include "MeshIO.h"
#include "Solver.h"
#include "Solution.h"
//---> Namespace includes: These are header files that define namespaces
#include "system_module.h" /* Includes namespace named:system_module */
#include "setup.h" /* Includes namespace named:setup */
#include "SolGlobalVars.h"


//****************************************************************************80
//! \brief print_main : Prints the welcome screen to the user
//! \nick
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//****************************************************************************80
void print_main()
{
  
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << " ++++ ++++ ++       ++      +++  ++++ +    " << std::endl;
  std::cout << " +    +    + +     + +      +    +  + +    " << std::endl;
  std::cout << " ++++ ++++ +  +   +  + ---   ++  +  + +    " << std::endl;
  std::cout << " +    +    +   + +   +         + +  + +    " << std::endl;
  std::cout << " +    ++++ +    +    +      +++  ++++ ++++ " << std::endl;
  std::cout << std::endl;
  std::cout << "            FEM_SOL: Solver for FEM Code 2D " 
	    << std::endl;
  std::cout << "                Version 1.0 " << std::endl;
  std::cout << "                 Authored by  " << std::endl ;
  std::cout << "     Nicholas Burgess (nicholas.k.burgess.ctr@us.army.mil) " 
	    << std::endl;
  std::cout << "       Ryan Glasby (ryan.glasby.ctr@arnold.af.mil)" 
	    << std::endl;
  std::cout << std::endl;
  std::cout << "     U.S. Army Aeroflightdynamics Directorate   " << std::endl;
  std::cout << "             NASA Ames Research Center         " << std::endl;
  std::cout << "                      MS 215-1              " << std::endl;
  std::cout << "              Moffett Airfield, CA 94035     " << std::endl;
  std::cout << std::endl; 
  std::cout << std::endl;
  std::cout << std::endl;

}

//****************************************************************************80
//! \brief main : Main program of fem_sol.  
//! \nick
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//****************************************************************************80
int main(int argc, char** argv)
{
  
  //---> Print the main screen for user. 
  print_main();
  
  //---> Parse the input
  if( argc != 2 ){// check-input
    /*---> If we didn't find 2 arguments then the user has made a mistake 
      executing the code and we need to adjust him/her.  */
    
    //---> Print error message
    std::cout << "Usage Error: you specified too few or too many arguments. " 
	      << std::endl 
	      << "Code has deteced the mistake please excute as fem_sol." 
	      << std::endl
	      << " {version-number} input-file." <<  std::endl 
	      << "Example usage for version 1 with input file named: " 
	      << std::endl
	      << "pre.input..." << std::endl
	      << " Usage: fem_sol.1.0 sol.input" << std::endl << std::endl;
    
#ifdef DEV_DEBUG
    std::cout << "Developer Debug: program exited in file: preproc/main.cpp "
	      << std::endl
	      << "at line 80.  The exit was tripped because the check on " 
	      << std::endl
	      << "line 64 failed.  Wrong number of input arguments was some-" 
	      << std::endl
	      << "how specified.  " << std::endl; 
#endif
    
    //---> Exit the code
    system_module::my_exit();
  
  } //END if check-input 

  //---> Instantiate input file reader
  InputIO<int,double> input_read(argv[1]);
  
  //---> Now read the inputfile
  input_read.readyamlinput();

  
  //---> Instatiate mesh file reader (MESHIO)
  MeshIO<int,double> mesh_reader; 
  
  //---> Initalize Solver Global Vars
  SolGlobalVars::initialize(input_read.global_data_map);
  double uinf = 0;
  double vinf = 0;
  
  //--------------------------- Instiate and Setup Grid ------------------------
  //---> Instatiate Solution and Grid
  UnstGrid<int,double> grid(2);
  
  //---> Read mesh
  if(SolGlobalVars::MeshType=="ryan-2D"){
    mesh_reader.read_2Dryan_ascii(SolGlobalVars::Project, grid); 
  }
  mesh_reader.read_bc(SolGlobalVars::Project, grid);
  
  //---> Now setup the rest of the connectivity
  grid.init_connectivity();
  grid.write_tecplot(SolGlobalVars::Project);
  
  //--------------------------- Euler Equations -------------------------------
  if( SolGlobalVars::EquationType.compare("Euler") == 0) {
    uinf = SolGlobalVars::Fmach*cos(pi/180.0*SolGlobalVars::Alpha);
    vinf = SolGlobalVars::Fmach*sin(pi/180.0*SolGlobalVars::Alpha);
    
    EulerPDESet2D<int,double> Euler(1.0, uinf, vinf, 1.0/1.4);
    Array1D<double> qinf(Euler.get_nfld());
    Euler.get_qinf(qinf);
    
    Solution<int,double> soln(grid, Euler.get_nfld(), 0, qinf, 
			      SolGlobalVars::PG);
    soln.write_tecplot(SolGlobalVars::Project);
    
    std::cout << "Euler Equations Selected proceeding with solution..." 
	      << std::endl;
    
    Solver<int, double, EulerPDESet2D<int, double> > solver(Euler, grid, soln);
    solver.SteadySolve(SolGlobalVars::NTS);
     
  }
  
  //---------------------- Navier-Stokes Equations -----------------------------

  //---------------------------- Wrong Equations -------------------------------
  else {
    std::cout << "ERROR: Unsupported EquationType: "<< 
      SolGlobalVars::EquationType << std::endl;
  } // End Equation Type
     
  std::cout << "Solution Completed...Good Bye!" 
	      << std::endl;
}// END main
