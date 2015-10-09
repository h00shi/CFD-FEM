#include "my_incl.h"
#include "UnstGrid.h" /* Definition of class UnstGrid */
#include "InputIO.h" /* Definition of class InputIO */
#include "MeshIO.h" /* Definition of class MeshIO */

//---> Namespace includes: These are header files that define namespaces
#include "system_module.h" /* Includes namespace named:system_module */
#include "PreGlobalVars.h" /* Includes namespace named:PreGlobalVars */

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
  std::cout << "  ++++  ++++   ++++      ++++ ++++ ++       ++ " << std::endl;
  std::cout << "  +   + +   +  +         +    +    + +     + + " << std::endl;
  std::cout << "  ++++  ++++   ++++ ---  ++++ ++++ +  +   +  + " << std::endl;
  std::cout << "  +     +   +  +         +    +    +   + +   + " << std::endl;
  std::cout << "  +     +    + ++++      +    ++++ +    +    + " << std::endl;
  std::cout << std::endl;
  std::cout << "            PRE_FEM: Processor for FEM Code " 
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
//! \brief main : Main program of pre_fem.  
//! \nick
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//****************************************************************************80

int main(int argc, char** argv)
{
  
  MeshIO<int,double> mesh_read; // Instantiate a mesh reader class
  
  //---> Print the main screen for user. 
  print_main();
  
  //---> Parse the input
  if( argc != 2 ){// check-input
    /*---> If we didn't find 2 arguments then the user has made a mistake 
      executing the code and we need to adjust him/her.  */
    
    //---> Print error message
    std::cout << "Usage Error: you specified too few or too many arguments. " 
	      << std::endl 
	      << "Code has deteced the mistake please excute as pre_supg." 
	      << std::endl
	      << " {version-number} input-file." <<  std::endl 
	      << "Example usage for version 1 with input file named: " 
	      << std::endl
	      << "pre.input..." << std::endl
	      << " Usage: pre_supg.1.0 pre.input" << std::endl << std::endl;
    
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

  InputIO<int,double> input_read(argv[1]); // Instantiate an input reader class

  //---> Read the inputfile
  input_read.readyamlinput();

  //---> Initialize Global Variables
  PreGlobalVars::initialize(input_read.global_data_map);
  
  UnstGrid<int,double> grid(3); // Instantiate mesh using UnstGrid class  
  
  //---> Open mesh file for reading
  if( PreGlobalVars::meshtype == "fv-uns") {
    //---> Mesh is in fv-uns format open this file
    mesh_read.read_fvuns_ascii(PreGlobalVars::project,grid);
  }
  
  //---> Now setup the rest of the connectivity
  grid.init_connectivity();
  
}// END int main
