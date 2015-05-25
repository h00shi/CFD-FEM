// -*-c++-*-
#ifndef MESHIO_H
#define MESHIO_H
//****************************************************************************80
//! \file MeshIO.h 
//! \class MeshIO MeshIO.h
//! \brief This is the header file defining the class MeshIO. 
//! \details The functions in this class are designed to read various mesh
//!          formats.  
//! \nick
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//****************************************************************************80
#include "my_incl.h" // Must always include my_kinddefs.h
#include <fstream> // Need fstream for file io
#include "system_module.h"
#include "UnstGrid.h" //Need unstructured grid class 

template < typename intT, typename realT>
class MeshIO {
  
private:
  //+++++++++++++++++++++++++++++++ PRIVATE STUFF ++++++++++++++++++++++++++++++
  std::ifstream mesh_file; /*!< Mesh file class, used to access the mesh 
			     file via fstream operation */
  std::ifstream bc_file; /*!< Boundary condition file object */
  //++++++++++++++++++++++++++++++ Class Private Data ++++++++++++++++++++++++++
  
public: 
  //+++++++++++++++++++++++++++++++ PUBLIC STUFF +++++++++++++++++++++++++++++++
  
//****************************************************************************80
//! \brief MeshIO : Is the constructor for this class, it's called 
//!                      automatically upon instantiation of the class
//! \details This function is just going to initialize all the data of the 
//!          class.  
//! \nick
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//****************************************************************************80
  MeshIO() 
  {
     
  }// End UnstGrid
  
//****************************************************************************80
//! \brief ~MeshIO : Is the destructor for this class, it's called 
//!                      automatically upon instantiation of the class
//! \details This function is just going to initialize all the data of the 
//!          class.  
//! \nick
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//****************************************************************************80
  ~MeshIO() 
  {
    
  }// End ~MeshIO
  
//****************************************************************************80
//!
//! \brief read_fvuns_ascii : Opens and gather initial data on a field view 
//!                            unstructured file.  Then goes back and reads the
//!                           mesh.   
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] project The project name stored in a std::string
//! \param[out] thisgrid The grid whose values are being populated
//****************************************************************************80
  void read_fvuns_ascii(const std::string project, 
			UnstGrid<int,double>& thisgrid)
  
  {
    //---> Local Variables
    std::string filename = project + std::string(".uns");
    std::string crap;
    
    std::cout << "Opening mesh file: " << filename << std::endl;	\
    /*---> Open the mesh file.  NOTE: must cast string to char* because that 
      is the argument of mesh_file.open(); */
    
    mesh_file.open( (char*)filename.c_str());
    
    if( mesh_file.is_open() != true ) { // file_check
     
      std::cout << "ERROR: Could not open mesh file: " << filename 
		<< std::endl;
      system_module::my_exit();
      
    } //End file_check
    
    //-------------------------------------------------------------------------
    //---> Read 9 lines of stupid crap
    for (intT i = 0; i < 9; i++) { //crap_loop 
      std::getline(mesh_file, crap);
    }  // End crap_loop 
   
    //--------------------------------------------------------------------------
    //---> Read number of boundary id's
    std::getline(mesh_file,crap);
    sscanf(crap.c_str(), "%d", &thisgrid.nbc_id);
       
    for ( intT i = 0; i < thisgrid.nbc_id; i++ ) { // bc_id loop
      std::getline(mesh_file, crap);
          
    }// End bc_id loop 
    
    //--------------------------------------------------------------------------
    //---> Read 4 lines of stupid crap
    for (intT i = 0; i < 4; i++) { //crap_loop 
      std::getline(mesh_file, crap);
    }  // End crap_loop 
    
    //--------------------------------------------------------------------------
    //---> Read the number of nodes
    std::getline(mesh_file,crap);
    sscanf(crap.c_str(), "%d", &thisgrid.nnode);
     
    //--------------------------------------------------------------------------
    //---> Read nnode+1 lines of stupid crap
    for (intT i = 0; i < thisgrid.nnode+1; i++) { //crap_loop 
      std::getline(mesh_file, crap); 
    }  // End crap_loop 
    
    //--------------------------------------------------------------------------
    //---> Read the number of boundary faces
    std::getline(mesh_file,crap);
    sscanf(crap.c_str(), "%d", &thisgrid.nbc_face);
       
    //---> Read through boundary faces to get the size of bc_face2node array
    thisgrid.bc_f2n_size = 0;
    //---> Loop over boundary faces
    for (intT i = 0; i < thisgrid.nbc_face; i++) {
      //---> Read bc_id into crap for now
      std::getline(mesh_file,crap);
      intT icrap;
      intT nnf;
      
      sscanf(crap.c_str(), "%d %d", &icrap, &nnf);
      //std::cout << icrap <<  "  " << nnf << std::endl;
      thisgrid.bc_f2n_size += nnf;
    }
    
    //---> Read a line of crap
    std::getline(mesh_file, crap);
    crap = "NULL";
    
    //--------------------------------------------------------------------------
    //---> Read the elements
    thisgrid.nelement = 0;
    thisgrid.e2n_size = 0;
    while( crap != "Variables") {
      
      intT nne;
      intT etype;
      intT icrap;
      
      //---> Read First bit of the line
      getline(mesh_file,crap);
      
      if( crap != "Variables") {
	//---> Increment the number elements
	thisgrid.nelement +=1;
	
	// Parse line knowing second entry is the element type
	sscanf(crap.c_str(), "%d %d", &icrap, &etype);
	
	//---> Based on etype set number of nodes per element
	switch (etype){
	case 1 : 
	  nne = 4;
	  break;
	case 2 :
	  nne = 5;
	  break;
	case 3 : 
	  nne = 6;
	  break;
	case 4 :
	  nne = 8;
	  break;
	default :
	  nne = 0;
	  break;
	}
      
	thisgrid.e2n_size +=nne;

      }
    
    }
    
    std::cout << "Mesh Information: " << std::endl;
    std::cout << "  NODES: " << thisgrid.nnode << std::endl;
    std::cout << "  NELEMENT: " << thisgrid.nelement << std::endl;
    std::cout << "  NBC_FACE: " << thisgrid.nbc_face << std::endl;
    std::cout << "  NBC_ID: " << thisgrid.nbc_id << std::endl;
#ifdef DEV_DEBUG
    std::cout << "  E2N_SIZE: " << thisgrid.e2n_size << std::endl;
    std::cout << "  BC_F2N_SIZE: " << thisgrid.bc_f2n_size << std::endl;
#endif
    std::cout << std::endl;
    
    //---> Allocate the mesh variables for reading
    thisgrid.allocate();
    
    /*---> These 2 operations rewind the file to the beginning, you know 
      F90 can do this with a single cmd??? */
    mesh_file.clear();
    mesh_file.seekg(mesh_file.beg);
    
   //-------------------------------------------------------------------------
    //---> Read 9 lines of stupid crap
    for (intT i = 0; i < 9; i++) { //crap_loop 
      std::getline(mesh_file, crap);
    }  // End crap_loop 
   
    //--------------------------------------------------------------------------
    //---> Read number of boundary id's
    std::getline(mesh_file,crap);
    for ( intT i = 0; i < thisgrid.nbc_id; i++ ) { // bc_id loop
      //---> Read boundary line out of file into string called crap
      std::getline(mesh_file, crap);
      //---> Now set the tag for the bc_id number specified by i
      thisgrid.set_bcid_tag(i, crap);
      
    }// End bc_id loop 
    
    //--------------------------------------------------------------------------
    //---> Read 5 lines of stupid crap
    for (intT i = 0; i < 5; i++) { //crap_loop 
      std::getline(mesh_file, crap);
    }  // End crap_loop 
    
    //--------------------------------------------------------------------------
    //---> Read Nodal coordinates
    for (intT i = 0; i < thisgrid.nnode; i++) {// read_node_coords_loop
      realT x, y, z;
      //---> Read line
      std::getline(mesh_file, crap);
      //---> Parse line into x y z
      sscanf(crap.c_str(), "%lf %lf %lf", &x, &y, &z);
      
      thisgrid.set_node_coord(i, 0, x);
      thisgrid.set_node_coord(i, 1, y);
      thisgrid.set_node_coord(i, 2, z);
    }// End read_node_coords_loop
   
    //--------------------------------------------------------------------------
    //---> Read 2 lines of crap
    for (intT i = 0; i < 2; i++) { //crap_loop 
      std::getline(mesh_file, crap);
    }  // End crap_loop 
    
    //-------------------------------------------------------------------------
    for (intT i = 0; i < thisgrid.nbc_face; i++) {// bc_face
      intT id, nnf;
      intT node;
      //---> Get the line 
      std::getline(mesh_file, crap);
      
      //---> Use the line as the constructor of string stream
      std::istringstream ss(crap);
      
      //---> Read id from ss
      ss >> id;
      
      //---> Read number of nodes per face from ss
      ss >> nnf;
      
      //---> Set the number of nodes for this boundary face
      thisgrid.set_nnode_on_bc_face(i, nnf);
      thisgrid.set_bc_face_id(i, id - 1);
      //---> Read nodes on face into grid class
      for (intT n = 0; n < nnf; n++) { // nodes_on_bc_face
	//---> Read the node id's out of string stream
	ss >> node;
	
	thisgrid.add_node_to_bc_face(i,n,node - 1);
	
      } // End nodes_on_bc_face
      
           
    }// End bc_face
    //--------------------------------------------------------------------------
    //---> Read 1 line of crap
    for (intT i = 0; i < 1; i++) { //crap_loop 
      std::getline(mesh_file, crap);
    }  // End crap_loop 
    
    //-------------------------------------------------------------------------
    for (intT i = 0; i < thisgrid.nelement; i++) {// element
      intT etype, nne, icrap;
      intT node;
      //---> Get the line 
      std::getline(mesh_file, crap);
      
      //---> Use the line as the constructor of string stream
      std::istringstream ss(crap);
      
      //---> Read id from ss
      ss >> icrap;
      
      //---> Read number of nodes per face from ss
      ss >> etype;
      
      //---> Based on etype set number of nodes per element
      switch (etype){
      case 1 : 
	nne = 4;
	//---> Set element type tetrahedra
	thisgrid.set_element_type(i, 3);
	break;
      case 2 :
	nne = 5;
	//---> Set element type pyramid
	thisgrid.set_element_type(i, 5);
	break;
      case 3 : 
	nne = 6;
	//---> Set element type prism
	thisgrid.set_element_type(i, 4);
	break;
      case 4 :
	nne = 8;
	//---> Set element type hex
	thisgrid.set_element_type(i, 6);
	break;
      default :
	nne = 0;
	break;
      }

      //---> Set the number of nodes for this boundary face
      thisgrid.set_nnode_on_element(i, nne);
      
      //---> Read nodes on element into grid class;
      for (intT n = 0; n < nne; n++) { // nodes_on_element
	//---> Read the node id's out of string stream
	ss >> node;
	
	thisgrid.add_node_to_element(i, n, node - 1);
	
      } // End nodes_on_element
      
      /*---> Fieldview does not define elements the way that I would define 
	them.  In the case of tets they are right handed just not in the way, 
	I do it.  In the case of hexes it's just plain screwy what they do, 
	it looks like it's a throw back to structured grids...ughhh.  
	Regardless, here we get things in order so we have the correct node 
	ordering for our code.  */
      switch(etype){
      case 1 : 
	//---> Store node 1 in temp for swapping
	intT temp0 = thisgrid.get_node_on_element(i, 0);
       
	//---> Write node 3 into node 0's location
	thisgrid.add_node_to_element(i, 0, thisgrid.get_node_on_element(i, 3) );
	
	//---> Write node 0 into node'3 location
	thisgrid.add_node_to_element(i, 3, temp0);

	break;
      }
    
    } //End element

    std::cout << "DONE Reading mesh file: " << filename << std::endl
	      << std::endl;
    mesh_file.close();
  } //End open_fvuns_ascii

//****************************************************************************80
//!
//! \brief read_2Dryan_ascii : Opens and gather initial data on a  2D Ryan
//!                            unstructured file.  Then goes back and reads the
//!                           mesh.   
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] project The project name stored in a std::string
//! \param[out] thisgrid The grid whose values are being populated
//****************************************************************************80
  void read_2Dryan_ascii(const std::string project, 
			UnstGrid<int,double>& thisgrid)
  
  {
    //---> Local Variables
    std::string filename = project + std::string(".mesh");
    std::string crap;
    
    std::cout << "Opening mesh file: " << filename << std::endl;	
    /*---> Open the mesh file.  NOTE: must cast string to char* because that 
      is the argument of mesh_file.open(); */
    
    mesh_file.open( (char*)filename.c_str());
    
    if( mesh_file.is_open() != true ) { // file_check
     
      std::cout << "ERROR: Could not open mesh file: " << filename 
		<< std::endl;
      system_module::my_exit();
 
    } //End file_check

    //---> Read 1 line of stupid crap
    std::getline(mesh_file, crap);
    intT wjunk, xjunk, yjunk, zjunk;
    
    //---> Read line
    std::getline(mesh_file, crap);
    
    //---> Parse line into w x y z
    sscanf(crap.c_str(), "%d %d %d %d", &wjunk, &xjunk, &yjunk, &zjunk);
    thisgrid.nnode = wjunk;
    thisgrid.nelement = xjunk;
    thisgrid.nbc_id = yjunk;
    thisgrid.nbc_face = zjunk;
      
    //---> Read 1 line of stupid crap
    std::getline(mesh_file, crap);
    
    //---> Read 1 line of stupid crap --> reread #nodes and do nothing
    std::getline(mesh_file, crap);
      
    thisgrid.bc_f2n_size = thisgrid.nbc_face*2;
    thisgrid.e2n_size = thisgrid.nelement*3;
    thisgrid.allocate();
    std::cout << "Mesh Information: " << std::endl;
    std::cout << "  NODES: " << thisgrid.nnode << std::endl;
    std::cout << "  NELEMENT: " << thisgrid.nelement << std::endl;
    std::cout << "  NBC_FACE: " << thisgrid.nbc_face << std::endl;
    std::cout << "  NBC_ID: " << thisgrid.nbc_id << std::endl;
#ifdef DEV_DEBUG
    std::cout << "  E2N_SIZE: " << thisgrid.e2n_size << std::endl;
    std::cout << "  BC_F2N_SIZE: " << thisgrid.bc_f2n_size << std::endl;
#endif
    std::cout << std::endl;
    //---> Read Nodal coordinates
    for (intT i = 0; i < thisgrid.nnode; i++) {// read_node_coords_loop
      realT x, y;
      //---> Read line
      std::getline(mesh_file, crap);
    
      //---> Parse line into x y 
      sscanf(crap.c_str(), "%lf %lf ", &x, &y);
      
      thisgrid.set_node_coord(i, 0, x);
      thisgrid.set_node_coord(i, 1, y);
    }// End read_node_coords_loop

    //---> Read 1 line of stupid crap
    std::getline(mesh_file, crap);
   
    //---> Read 1 line of stupid crap --> reread #triangles and do nothing
    std::getline(mesh_file, crap);
    
    for (intT i = 0; i < thisgrid.nelement; i++) {// element
      //---> Get the line 
      std::getline(mesh_file, crap);

      sscanf(crap.c_str(), "%d %d %d", &wjunk, &xjunk, &yjunk);
    
      //---> assumed all triangle grid
      thisgrid.set_element_type(i, 1); 

      //---> Set the number of nodes for this element
      thisgrid.set_nnode_on_element(i, 3);
    
      // decrease the node #'s by 1 for c formatting 
      thisgrid.add_node_to_element(i, 0, wjunk - 1);
      thisgrid.add_node_to_element(i, 1, xjunk - 1);
      thisgrid.add_node_to_element(i, 2, yjunk - 1);
    } //End element
   
    //---> Read 1 line of stupid crap
    std::getline(mesh_file, crap);
    
    //---> Read 1 line of stupid crap
    std::getline(mesh_file, crap);
    
    //---> Read 1 line of stupid crap
    std::getline(mesh_file, crap);
   
    for (intT i = 0; i < thisgrid.nbc_face; i++) {// boundary loop
      //---> Get the line 
      std::getline(mesh_file, crap);

      sscanf(crap.c_str(), "%d %d %d", &wjunk, &xjunk, &yjunk);

      thisgrid.set_bc_face_id(i, wjunk - 1);
      
      
      //---> Set the number of nodes for this boundary face
      //---> assume 2d mesh with bar elements on boundary
      thisgrid.set_nnode_on_bc_face(i, 2);
      
      thisgrid.add_node_to_bc_face(i,0,xjunk - 1);
      thisgrid.add_node_to_bc_face(i,1,yjunk - 1);
     } // end boundary loop

    std::cout << "DONE Reading mesh file: " << filename << std::endl
	      << std::endl;
    mesh_file.close();
    
  } //End open_2Dryan_ascii

//****************************************************************************80
//!
//! \brief read_bc : Reads the boundary condition "project".bc
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] project The project character string
//! \param[out] thisgrid A reference to the grid instance you want to read into
//****************************************************************************80
  void read_bc(const std::string project, UnstGrid<intT,realT>& thisgrid)
  {
    //---> Local Variables
    std::string filename = project + std::string(".bc");
    std::string crap;
    
    std::cout << "Opening bc file: " << filename << std::endl;	
    /*---> Open the mesh file.  NOTE: must cast string to char* because that 
      is the argument of mesh_file.open(); */
    
    bc_file.open( (char*)filename.c_str());
    
    if( bc_file.is_open() != true ) { // file_check
      
      std::cout << "ERROR: Could not open boundary condition file: " 
		<< filename 
		<< std::endl;
      system_module::my_exit();
 
    } //End file_check
    
    //---> Read 1 line of stupid crap
    for(intT i = 0; i < thisgrid.nbc_id; i++){// ID loop 
      //---> Read line from file into crap
      std::getline(bc_file, crap);
      //---> Some helpful variables
      intT id, bctype;
      char tag[100];
      //---> Use sscanf to parse string crap
      sscanf(crap.c_str(), "%d %d %s", &id, &bctype, tag);
      /*---> Instantiate a variable of type string called stag with tag as 
	it's value */
      std::string stag(tag);
     
      //---> Populate grid structure
      thisgrid.set_bcid_tag(id, stag);
      thisgrid.set_bcid_type(id,bctype);
      std::cout << "Boundary " << id << " " << stag << " BC-TYPE: " 
		<< bctype << std::endl; 
    }// End ID loop 
    
    
    std::cout << "Done reading bc file." << std::endl;
    
    bc_file.close();
  }// End read_bc

#endif
};
