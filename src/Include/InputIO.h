// -*-c++-*-
#ifndef INPUTIO_H
#define INPUTIO_H
//****************************************************************************80
//! \class InputIO InputIO.h
//! \brief This is the header file defining the class InputIO
//! \details The functions and data present here are designed to parse an 
//!          a particular form of input file for use in this code.    
//! \nick
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//****************************************************************************80

#include "my_incl.h" // Must always include my_kinddefs.h
#include <fstream> // Need fstream for file io
#include "yaml-cpp/yaml.h"
//---> Namespace includes: These are header files that define namespace
#include "system_module.h" /* Includes namespace named:system_module */

template < typename intT, typename realT>
class InputIO 
{
private:
  //+++++++++++++++++++++++++++++++ PRIVATE STUFF +++++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++ Class Private Data +++++++++++++++++++++++++
  std::ifstream infile; /*!< File object used to access file */
  std::string line; /*!< A string to represent the line */
  char* inputfilename; /*!< The name of the input file */
//****************************************************************************80
//!
//! \brief InputIO () : Default constructor of the InputIO
//! \details Blocked...you should know your input file name before instantiating
//! \nick 
//! \version $Rev$
//! \date $Date$
//! 
//****************************************************************************80
  InputIO(){}
public: 
  //+++++++++++++++++++++++++++++++ PUBLIC STUFF ++++++++++++++++++++++++++++++
  std::map<std::string,std::string> global_data_map; /*!< A map of the global 
						       data to be decoded 
						       during global parameter 
						       initialization */

//****************************************************************************80
//!
//! \brief InputIO : Constructor, initializes memeber inputfilename
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] filename The filename to initialize infile with
//****************************************************************************80
  InputIO(char* filename) : inputfilename(filename) {}

//****************************************************************************80
//!
//! \brief readyaminput : Reads the YAML bases input file to the Code
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! 
//****************************************************************************80
  void readyamlinput()
  {
    //---> Local Variables
    YAML::Iterator di ; /*!<  "document" iterator */
    YAML::Parser yaml_parser; /*!< YAML parser object, you can think of it as
				a "reader" function */
    YAML::Node document; /*!< YAML document, YAML allows multiple for a file
			   we will keep it to 1 */
    std::string varname, value, dbname ; /*!< For reading file */
    intT i; /*!< Iteration index */
   
    
    //---> Initialize the std::map global_data_map to clear so that we know
    global_data_map.clear();
    
    //---> Tell the user what document we are reading
    std::cout << "Reading input file: " << inputfilename << std::endl;
    
    //---> Open the file for reading
    infile.open(inputfilename);
    if( infile.is_open() != true) {// input_check
      std::cout << "ERROR: Could not open input file: " << inputfilename 
		<< std::endl;
      system_module::my_exit();
    } // End input_check
    
    //---> Tell the yaml parser to read the file stored in infile
    yaml_parser.Load(infile);
    
    //---> Parse the databases of the file
    yaml_parser.GetNextDocument(document);

    //---> Loop over databases of document and get database name
    for (di = document.begin(); di != document.end(); ++di) { // Document_loop
      /*---> Using iterator first function get the database name */ 
      di.first() >> dbname ;
      
      std::cout << "Reading Database: " << dbname << std::endl;
      
      /*---> Define a reference map to the database which is also of type 
	YAML::Node.  This acts as pointer to the database so we can extract 
	data from it. */   
      const YAML::Node& db = document[dbname];
      for(YAML::Iterator it = db.begin(); it != db.end(); ++it) {
	it.first()  >> varname;
	it.second() >> value;
	
	std::cout << "  " << varname << " = " << value << std::endl;
	
	//---> Load the varname as the map key and the value into the map
	global_data_map[varname] = value;
	
      }
      std::cout << std::endl;
    } //End Document_loop
    
    
  
  } //End readyamlinput

  
};// End class InputIO
#endif
 
