// -*-c++-*-
#ifndef UNSTGRID_H
#define UNSTGRID_H
//****************************************************************************80
//! \class UnstGrid 
//! \brief This is the header file defining the class UnstGrid
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \tparam intT Template argument meant to mimic integer
//! \tparam realT Template argument meant to mimic real numbers
//****************************************************************************80
#include "my_incl.h"
#include "system_module.h"
#include <map> /*!< I need the std::map type */
#include<vector> 
#include<algorithm>
#include<map> 
#include<numeric>
#include<fstream>
#include "TriElement.h"
#include "TetElement.h"
#include "Array2D.h"

template < typename intT, typename realT >  
class UnstGrid {
  
private: 
  //+++++++++++++++++++++++++++++++ PRIVATE STUFF ++++++++++++++++++++++++++++++
  
  //++++++++++++++++++++++++++++++ Class Private Data ++++++++++++++++++++++++++
  //---> UnstGrid Description Pointers 
  intT* element2nodei; //!< Indexing array for nodes of element 
  intT* element2node; /*!< Element to node array, gives indicies of nodes on an
		     element */
  intT* bc_face2element; /*!< Boundary face to element array, for a boundary 
			   face this gives the element number attached 
			   to it. */  
  intT* face2nodei; //!< Indexing array for nodes of a face
  intT* face2node; //!< Gives nodes on a face for all interior faces
  intT* bc_face2nodei; //!< Indexing array for nodes of a boundary face
  intT* bc_face2node; /*!< Gives nodes on a boundary face for all boundary 
			faces */
  intT* element_type; /*!< Flag for the type of element
		     \verbatim
		     0 = 1-D Bar, a line segement with end points
		     1 = 2-D Triangle, 3 verticies connected together
		     2 = 2-D Quad, 4 verticies connected together
		     3 = 3-D Tetrahedra, 4 verticies with 4 triangular faces
		     4 = 3-D Prism, 6 vertices with 2 trianglular faces and
		     3 quadrilateral faces
		     5 = 3-D Pyramid, 5 verticies with 4 triangular faces and 
		     1 quadrilateral face
		     6 = 3-D Hexahedra, 8 verticies with 6 quadrilateral faces 
		     \endverbatim
		   */
  intT* element_p; //!< Flag for the element interpolation order 
  intT* bc_face_id; //!< For each boundary face this stores the id 
  intT* bcid_type; /*!< For each boundary id (or patch if you prefer) gives the
		     boundary condition type */
  intT* bc_local_face; /*!< Local face index for boundary faces */
  intT* nface_per_bcid; /*!< Gives number of faces for per boundary id # */
  intT* node2elementi; /*!< Index array for node 2 element linked list */
  intT* node2element; /*!< Data array for node 2 element linked list */
  intT* node_adji; /*!< Node adjacency linked list index array */
  intT* node_adj; /*!< Node adjacency linked list data array */
  
  //---> UnstGrid real values 
  realT* x; /*!< Coordinates, stores all the coordinates depending on 
	      number of dimensions.  
	      \verbatim
	      1-D : (x) 
	      2-D : (x,y)
	      3-D : (x,y,z)
	      \endverbatim 
	    */
  
  realT* element_vol; /*!< Notional volume of elements, 
		     different units for different numbers of dimensions */

  //---> Memory diagnostics data
  realT grid_mem ; /*!< Total memory for this instance of class */ 
  //---> Misc
  std::map<intT,std::string> bcid_tag; /*!< Boundary id string tag */
  
  //-------------------------- Local Face patterns -----------------------------
  //---> Tetrahedra 
  // std::vector<intT> tet_face1o1(intT, 3,-1);
  // std::vector<intT> tet_face1o2(3,-1);
  // std::vector<intT> tet_face1o3(3,-1);

  // std::vector<intT> tet_face2o1(3,-1);
  // std::vector<intT> tet_face2o2(3,-1);
  // std::vector<intT> tet_face2o3(3,-1);

  // std::vector<intT> tet_face3o1(3,-1);
  // std::vector<intT> tet_face3o2(3,-1);
  // std::vector<intT> tet_face3o3(3,-1);

  // std::vector<intT> tet_face4o1(3,-1);
  // std::vector<intT> tet_face4o2(3,-1);
  // std::vector<intT> tet_face4o3(3,-1);
  
  
//****************************************************************************80
//! \brief find_bc_elem : Given a boundary face figure out which element it 
//!         belongs to. 
/*! \verbatim
   Table of nodes defining a face for supported element types
   Tetrahedra:
    Face 0 has Nodes: 0 1 3 
    Face 1 has Nodes: 1 2 3
    Face 2 has Nodes: 2 0 3
    Face 3 has Nodes: 0 2 1
   Prism:
    Face 0 has Nodes: 0 2 1 
    Face 1 has Nodes: 0 3 4 2
    Face 2 has Nodes: 0 1 4 3
    Face 3 has Nodes: 1 2 4 5 
    Face 4 has Nodes: 3 4 5 
   Pyramid
    Face 0 has Nodes: 0 3 2 1
    Face 1 has Nodes: 0 1 4
    Face 2 has Nodes: 0 4 3
    Face 3 has Nodes: 3 4 2
    Face 4 has Nodes: 1 2 4 
   Hexahedron
    Face 0 has Nodes: 0 3 2 1 
    Face 1 has Nodes: 0 4 7 3 
    Face 2 has Nodes: 0 1 5 4
    Face 3 has Nodes: 1 2 6 5
    Face 4 has Nodes: 3 7 7 2
    Face 5 has Nodes: 4 5 6 7  
    \endverbatim */
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] f The boundary face
//! \return e The element to which boundary face f belongs
//****************************************************************************80
  intT find_bc_elem( const intT & f)
  {
    //---> Local Variables
    intT i; // Looping index
    intT j; // Looping index
    intT k; // Looping index
    intT count = 0; // Counter variable counting
    intT e; // Return value
    
    //---> Counter number of boundary nodes
    intT nbnode = bc_face2nodei[f + 1] - bc_face2nodei[f];
         
    //---> Get first node on boundary face
    intT bnode1 = bc_face2node[bc_face2nodei[f]];
    e = -1;
    
    //---> Loop over elements attached to node
    for( i = node2elementi[bnode1]; i < node2elementi[bnode1 + 1]; i++) { // Element on node loop
      
      //---> Get the ith element attached to node1
      intT elem = node2element[i];
      count = 1;
      
      //---> Loop over the rest of the nodes on bc_face f
      for ( j = bc_face2nodei[f] + 1; j < bc_face2nodei[f + 1]; j ++){ // bc_node_loop 
	intT node = bc_face2node[j];
	//---> Now loop over all the elements of node 
	for( k = node2elementi[node]; k < node2elementi[node + 1]; k++){ // node loop 
	  // Get the elements containing node 
	  intT elem2 = node2element[k];
	  /*---> Check to see if the element containing node 2 is the same 
	    as elem */
	  if(elem2 == elem) { //end if elem check
	    /*---> Since elem2 is the same as elem...then elem also contains *
	      node node.  Thefefore elem contains count + 1 nodes of the 
	      face f */	    
	    count += 1; 
	    /*---> Leave loop because we don't need to check anymore elements 
	      containg the node denoted node */
	    break;
	  } // End if elem check 
	  
	} // End node loop 
	
      }// End bc node loop 
               
      if( count == nbnode) { // check count
	e = elem;
	break;
      }// End count check
      
    } // End element on node loop 
	
    //---> Some error reporting
    if( e == -1 ) { // End if error
      std::cout << "ERROR: Could not find an element that contained all " 
		<<"nodes of the boundary face " << f 
		<< ".  There ie probably something wrong with the " 
		<< "connectivity of the mesh.  " << std::endl;
#ifdef DEV_DEBUG
      std::cout << " The face "  << f << " contains the following nodes:" 
		<< std::endl;
      for( i = bc_face2nodei[f]; i < bc_face2nodei[f + 1]; i++){
	std::cout << " Node " << i - bc_face2nodei[f] << ": " 
		  << bc_face2node[i] << std::endl;
      }
#endif	
      
    } // End if error
    
    //---> Return the element containing face f
    return(e);
    
  } // Function get_bc_elem

//****************************************************************************80
//! \brief get_loc_face_1D : Finds the local face id number for 1-D elements
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] face_nodes The nodes the face
//! \param[in] elem_nodes The node on the cell in question
//! \return loc_face The local face number
//****************************************************************************80
  intT get_loc_face_1D(const std::vector<intT>& face_nodes, 
		       const std::vector<intT>& elem_nodes)
  {
    //---> Return variable
    intT loc_face = -99; //---> Initialized to - 1 to throw error if it's wrong
    
    //---> Since it's a 1-D element it's very easy, face nodes is only 1 node
    if( face_nodes[0] == elem_nodes[0] ) loc_face = 0;
    if( face_nodes[0] == elem_nodes[1] ) loc_face = 1;

    return(loc_face);
    
  } // End get_loc_face_1D

//****************************************************************************80
//! \brief get_loc_face_tri : Finds the local face id number for 2-D triangles
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] face_nodes The nodes the face
//! \param[in] elem_nodes The node on the cell in question
//! \return loc_face The local face number
//****************************************************************************80
  intT get_loc_face_tri(const std::vector<intT>& face_nodes, 
		       const std::vector<intT>& elem_nodes)
  {
    //---> Return variable
    intT loc_face = -99; //---> Initialized to -99 to throw error if it's wrong
    
    //---> Since it's a 2-D triangle pretty easy...
   
    //----------------------------- Face 1 ------------------------------------
    if( face_nodes[0] == elem_nodes[0] &&
	face_nodes[1] == elem_nodes[1]) { // Face 1 check
      //---> It's face 1 in right handed order
      loc_face = 0;
    }
    //---> Otherwise it could be left handed
    else if(face_nodes[0] == elem_nodes[1] &&
	    face_nodes[1] == elem_nodes[0]) {
      //---> Mark as negative, which means it's this face but wrong order
      loc_face = -1;
    } // End Face 1 check
    
    
    //----------------------------- Face 2 ------------------------------------
    if( face_nodes[0] == elem_nodes[1] &&
	face_nodes[1] == elem_nodes[2]) { // Face 2 check
      //---> It's face 2 in right handed order
      loc_face = 1;
    }
    //---> Otherwise it could be left handed
    else if(face_nodes[0] == elem_nodes[2] &&
	    face_nodes[1] == elem_nodes[1]) {
      //---> Mark as negative, which means it's this face but wrong order
      loc_face = -2;
    } // End Face 2 check
    
    //----------------------------- Face 3 ------------------------------------
    if( face_nodes[0] == elem_nodes[2] &&
	face_nodes[1] == elem_nodes[0]) { // Face 3 check
      //---> It's face 3 in right handed order
      loc_face = 2;
    }
    //---> Otherwise it could be left handed
    else if(face_nodes[0] == elem_nodes[0] &&
	    face_nodes[1] == elem_nodes[2]) {
      //---> Mark as negative, which means it's this face but wrong order
      loc_face = -3;
    } // End Face 3 check
   
    return(loc_face);
    
  } // End get_loc_face_tri

//****************************************************************************80
//! \brief get_loc_face_quad : Finds the local face id number for 2-D quads
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] face_nodes The nodes the face
//! \param[in] elem_nodes The node on the cell in question
//! \return loc_face The local face number
//****************************************************************************80
  intT get_loc_face_quad(const std::vector<intT>& face_nodes, 
			 const std::vector<intT>& elem_nodes)
  {
    //---> Return variable
    intT loc_face = -99; //---> Initialized to - 99 to throw error if it's wrong
    
    //---> Since it's a 2-D quad pretty easy...
    
    //----------------------------- Face 1 ------------------------------------
    if( face_nodes[0] == elem_nodes[0] &&
	face_nodes[1] == elem_nodes[1]) { // Face 1 check
      //---> It's face 1 in right handed order
      loc_face = 0;
    }
    //---> Otherwise it could be left handed
    else if(face_nodes[0] == elem_nodes[1] &&
	    face_nodes[1] == elem_nodes[0]) {
      //---> Mark as negative, which means it's this face but wrong order
      loc_face = -1;
    } // End Face 1 check
        
    //----------------------------- Face 2 ------------------------------------
    if( face_nodes[0] == elem_nodes[1] &&
	face_nodes[1] == elem_nodes[2]) { // Face 2 check
      //---> It's face 2 in right handed order
      loc_face = 1;
    }
    //---> Otherwise it could be left handed
    else if(face_nodes[0] == elem_nodes[2] &&
	    face_nodes[1] == elem_nodes[1]) {
      //---> Mark as negative, which means it's this face but wrong order
      loc_face = -2;
    } // End Face 2 check
    
    //----------------------------- Face 3 ------------------------------------
    if( face_nodes[0] == elem_nodes[2] &&
	face_nodes[1] == elem_nodes[3]) { // Face 3 check
      //---> It's face 3 in right handed order
      loc_face = 2;
    }
    //---> Otherwise it could be left handed
    else if(face_nodes[0] == elem_nodes[3] &&
	    face_nodes[1] == elem_nodes[2]) {
      //---> Mark as negative, which means it's this face but wrong order
      loc_face = -3;
    } // End Face 2 check
    
    //----------------------------- Face 4 ------------------------------------
    if( face_nodes[0] == elem_nodes[3] &&
	face_nodes[1] == elem_nodes[0]) { // Face 4 check
      //---> It's face 4 in right handed order
      loc_face = 3;
    }
    //---> Otherwise it could be left handed
    else if(face_nodes[0] == elem_nodes[0] &&
	    face_nodes[1] == elem_nodes[3]) {
      //---> Mark as negative, which means it's this face but wrong order
      loc_face = -4;
    } // End Face 4 check
    
    return(loc_face);
    
  } // End get_loc_face_quad

//****************************************************************************80
//! \brief get_loc_face_tet : Finds the local face id number for 3-D Tetrahedra
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] face_nodes The nodes the face
//! \param[in] elem_nodes The node on the cell in question
//! \return loc_face The local face number
//****************************************************************************80
  intT get_loc_face_tet(const std::vector<intT>& face_nodes, 
			const std::vector<intT>& elem_nodes)
  {
    //---> Return variable
    intT loc_face = -99; //---> Initialized to - 99 to throw error if it's wrong
    std::vector<intT> loc_nodes(3,-1);
    //---> Now we are in 3-D...not so easy anymore.  
    
    for(intT i = 0; i < elem_nodes.size(); i++) { // elem node loop 
      for(intT j = 0; j < face_nodes.size(); j++) {
	if( face_nodes[j] == elem_nodes[i] ) {
	  loc_nodes[j] = i;
	}				 
      }
    } // End elem node loop 
 
    intT sum = 0;
    for (intT i = 0; i < 3; i++){
      sum += loc_nodes[i];
    }
    /*---> For a tetrahedra the sum of the local node numbers for each face is
      unique.  Therefore we use this property to find the local faces numbers
      for a tet.  If the element is not a tet we're going to need to find a 
      cuter trick for those.  */
    if( sum == 4 ) loc_face = 0;
    if( sum == 6 ) loc_face = 1;
    if( sum == 5 ) loc_face = 2;
    if( sum == 3 ) loc_face = 3;
      
    return(loc_face);
  } // End get_loc_face_tet


public:
  
  //+++++++++++++++++++++++++++++++ PUBLIC STUFF +++++++++++++++++++++++++++++++
  //---> Integer data that helps us count and keep track of the sizes of things
  intT ndim; /*!< Number of physical dimensions of the mesh, currently only 
	       value s of 1, 2, or 3 are supported. */
  intT nelement; /*!< Number of elements in the mesh, elements are cannonical 
		shapes such as triangles (2-D) and tetrhera (3-D) */
  intT nnode; //!< Number of nodes in the mesh 
  intT nface; /*!< The number of internal faces for this mesh, internal faces
		have elements of the mesh on both sides of it. */
  intT nbc_face; /*!< Number of boundary faces for this mesh, boundary faces
		   have only 1 mesh element attached to it.  
		   The collection of boundary faces defines the boundary of 
		   the computational domain */
  intT nbc_id; /*!< Number of unique boundary id numbers for this mesh.  
		 Boundary faces are given an id tag in order to group them 
		 into logical sub-collections.  For example all boundary faces 
		 on the surface of a wing will be given the same bc_id number 
		 because they share boundary conditions and in an engineer 
		 sense define a component.*/
  intT ncolor; /*!< Number of colors used to color mesh for parallel solver */
  intT e2n_size; /*!< The size of elem2node */
  intT bc_f2n_size; /*!< The size of bc_face2node */
  intT f2n_size; /*!< The size of face2node */
  
//****************************************************************************80
//! \brief UnstGrid : Is the constructor for this class, it's called 
//!                      automatically upon instantiation of the class
//! \details This function is just going to initialize all the data of the 
//!          class.  
//! \author Nick Burgess
//!  NASA Ames Research Center
//!  Moffett Field, CA
//!  \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] n The number of dimensions for the instantiated class
//****************************************************************************80
  UnstGrid(intT n)
  {
    /*---> Scalar data initialization.  NOTE: If all data to UnstGrid.h
      the order below should follow the definition order of UnstGrid.h.  
      If you don't do this I will delete your source (NKB). */
    ndim = n;
    nelement = 0;
    nnode = 0;
    nface = 0;
    nbc_face = 0;
    nbc_id = 0;
    ncolor = 0;
    
    /*---> Initialize all the pointers of the class to NULL...it's just good 
      practice.  NOTE: If all data to UnstGrid.h
      the order below should follow the definition order of UnstGrid.h.  
      If you don't do this I will delete your source (NKB). */    
    
    //---> Integer data
    element2nodei   = NULL;
    element2node    = NULL;
    bc_face2element = NULL;
    face2nodei      = NULL;
    face2node       = NULL;
    bc_face2nodei   = NULL;
    bc_face2node    = NULL;
    element_type    = NULL;
    element_p       = NULL;
    bc_face_id      = NULL;
    bcid_type       = NULL;
    bc_local_face   = NULL;
    nface_per_bcid  = NULL;
    node2elementi   = NULL;
    node2element    = NULL; 
    node_adji       = NULL;
    node_adj        = NULL;
    
    //---> REAL data 
    x           = NULL;
    element_vol = NULL;

    //---> Misc.
    grid_mem = 0.0;
    
  } // END UnstGrid
  
//****************************************************************************80
//! \brief ~UnstGrid : Is the destructor for this class, it's called 
//!                      automatically upon scope exit of the class
//! \details This function is just going to make sure everything is clean
//!          by setting all scalars to zero and deleting all pointers.   
//! \author Nick Burgess
//!  NASA Ames Research Center
//!  Moffett Field, CA
//!  \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//****************************************************************************80
  ~UnstGrid()
  {
    
    /*---> Scalar data deletion.  NOTE: If all data to UnstGrid.h
      the order below should follow the definition order of UnstGrid.h.  
      If you don't do this I will delete your source (NKB).  */
    
    nelement = 0;
    nnode = 0;
    nface = 0;
    nbc_face = 0;
    nbc_id = 0;
    ncolor = 0;
    
    /*---> Delete all the pointers of the class ...it's just good 
      practice.  NOTE: If all data to UnstGrid.h
      the order below should follow the definition order of UnstGrid.h.  
      If you don't do this I will delete your source (NKB). */ 
    
    /*--->  Further stylistically we check for if( pointer != NULL) for
      extra clarity we take no short cuts...if I find shortcuts in your
      contributions they  will be deleted and sent back to you (NKB).  
      It is NOT acceptable to check for non-NULL pointer as any of 
      if(pointer). */  
    
    //---> Integer data
    if( element2nodei != NULL ) delete [] element2nodei;
    if( element2node != NULL ) delete [] element2node;
    if( bc_face2element != NULL ) delete [] bc_face2element;
    if( bc_face2nodei != NULL ) delete [] bc_face2nodei;
    if( bc_face2node != NULL ) delete [] bc_face2node;
    if( element_type != NULL ) delete [] element_type;
    if( element_p != NULL ) delete [] element_p;
    if( bc_face_id != NULL ) delete [] bc_face_id;
    if( bcid_type != NULL ) delete [] bcid_type;  
    if( bc_local_face != NULL) delete [] bc_local_face;
    if( nface_per_bcid != NULL) delete [] nface_per_bcid;
    if( node2elementi != NULL ) delete [] node2elementi;
    if( node2element != NULL ) delete [] node2element;
    if( node_adji != NULL) delete [] node_adji;
    if( node_adj  != NULL) delete [] node_adj;
    
    //---> REAL data 
    if ( x !=NULL ) delete [] x;
    if ( element_vol !=NULL ) delete [] element_vol;
    
    //---> Set all pointers to NULL after deleltion
    //---> Integer data
    element2nodei   = NULL;
    element2node    = NULL;
    bc_face2element = NULL;
    bc_face2nodei   = NULL;
    bc_face2node    = NULL;
    element_type    = NULL;
    element_p       = NULL;
    bc_face_id      = NULL;
    bcid_type       = NULL;
    bc_local_face   = NULL;
    nface_per_bcid  = NULL;
    node2elementi   = NULL;
    node2element    = NULL; 
    node_adji       = NULL;
    node_adj        = NULL;
    
    //---> Real data
    x           = NULL;
    element_vol = NULL;

    //---> Misc.
    grid_mem = 0.0;
  } //END ~UnstGrid

//****************************************************************************80
//! \brief allocate : Member function to allocate memory for the class members
//! \details This will allocate the memory for the class instance.  
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//****************************************************************************80
  void allocate()
  {
    realT mem; // Memory for a given member
    
    //--------------------------------------------------------------------------
    mem = system_module::alloc_mem<intT, intT, realT>(element2nodei, 
						      nelement + 1);
    grid_mem += mem;
    //--------------------------------------------------------------------------
    mem = system_module::alloc_mem<intT, intT, realT>(element2node, e2n_size);
    grid_mem += mem;
    //--------------------------------------------------------------------------
    mem = system_module::alloc_mem<intT, intT, realT>(bc_face2element, 
						      nbc_face);
    grid_mem += mem;
    //--------------------------------------------------------------------------
    mem = system_module::alloc_mem<intT, intT, realT>(bc_face2nodei, 
						      nbc_face + 1);
    grid_mem += mem;
    //--------------------------------------------------------------------------
    mem = system_module::alloc_mem<intT, intT, realT>(bc_face2node, 
						      bc_f2n_size);
    grid_mem += mem;
    //--------------------------------------------------------------------------
    mem = system_module::alloc_mem<intT, intT, realT>(element_type, 
						      nelement);
    grid_mem += mem;
    //--------------------------------------------------------------------------
    mem = system_module::alloc_mem<intT, intT, realT>(element_p, 
						      nelement);
    grid_mem += mem;
    //--------------------------------------------------------------------------
    mem = system_module::alloc_mem<intT, intT, realT>(bc_face_id, 
						      nbc_face);
    grid_mem += mem;
    //--------------------------------------------------------------------------
    mem = system_module::alloc_mem<intT, intT, realT>(bcid_type, 
						      nbc_id);
    grid_mem += mem;
    //--------------------------------------------------------------------------
    mem = system_module::alloc_mem<intT, intT, realT>(bc_local_face, 
						      nbc_face);
    grid_mem += mem;
    //--------------------------------------------------------------------------
    mem = system_module::alloc_mem<intT, intT, realT>(nface_per_bcid, 
						      nbc_id);
    grid_mem += mem;
    //--------------------------------------------------------------------------
    mem = system_module::alloc_mem<intT, intT, realT>(node2elementi, 
						      nnode + 1);
    grid_mem += mem;
    //--------------------------------------------------------------------------
    mem = system_module::alloc_mem<intT, intT, realT>(node_adji, 
						      nnode + 1);
    grid_mem +=mem;
    //--------------------------------------------------------------------------
    mem = system_module::alloc_mem<realT, intT, realT>(x, ndim*nnode);
    grid_mem += mem;
    //--------------------------------------------------------------------------
    mem = system_module::alloc_mem<realT, intT, realT>(element_vol, nelement);
    grid_mem += mem;
    //--------------------------------------------------------------------------
    
    //Initialize pointer indexing arrays
    element2nodei[0] = 0;
    bc_face2nodei[0] = 0;

    std::cout << "Grid is using:  " << grid_mem << " MB of memory " <<std::endl;
    std::cout << std::endl;
    
  }// End Function allocate

//****************************************************************************80
//! \brief set_bcid_tag : Puts the boundary tag string into the bcid_tag 
//!  std::map 
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] id The boundary id number
//! \param[in] tag The name to give it 
//****************************************************************************80
  void set_bcid_tag(const intT& id, const std::string& tag)
  {
    //---> Set the tag of boundary id number = tag
    bcid_tag[id] = tag;
  
  }// End Function set_bcid_tag
  
//****************************************************************************80
//! \brief get_bcid_tag : Returns the bcid_tag of the specified boundary number
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] id The boundary id number
//! \return tag The tag being returned to the user 
//****************************************************************************80
  std::string get_bcid_tag(const intT& id)
  {
    //---> Return variables
    std::string tag(bcid_tag[id]);
    
    return(tag);
  
  }// End Function get_bcid_tag

//****************************************************************************80
//! \brief set_node_coord : Sets the value of a node coordinate
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] node The node who's coordinate is being set indexed from 0
//! \param[in] dim The domain or physical dimension your are setting indexed 
//!            from 0
//! \param[in] coord Value of coord
//****************************************************************************80
  void set_node_coord(const intT& node, const intT& dim, const realT& coord)
  {
    //---> Local Variables
    intT index = (node)*ndim + dim; /*!< Flat array index for x */
    x[index] = coord;
  } // End Function set_node_coord;
  
//****************************************************************************80
//! \brief get_node_coord : Gets the value of a node coordinate
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] node The node who's coordinate is being set indexed from 0
//! \param[in] dim The domain or physical dimension your are setting indexed 
//!            from 0
//! \return coord The coordinate to be returned to the user
//****************************************************************************80
  realT get_node_coord(const intT& node, const intT& dim)
  {
    //---> Local Variables
    intT index = (node)*ndim + dim; /*!< Flat array index for x */
    realT coord = x[index];
    return(coord);
  } // End Function get_node_coord

//****************************************************************************80
//! \brief set_nnode_on_element : Tells the class how many nodes make up an
//!        specified element
//! \details The number of nodes is set via the following formula.  
//!   For an element e, \formula{ n = element2nodei[e+1] - element2nodei[e] }. 
//!   Thefore knowing n means we are setting the value of element2nodei[e+1].
//!   Effectively this function establishs how many nodes there are per 
//!   element and the element2nodei linked list indexing array at the 
//!   same time.      
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] e The element number
//! \param[in] n The number of nodes for this element
//****************************************************************************80
  void set_nnode_on_element(const intT& e, const intT& n)
  {
    //---> Add n nodes to element2nodei[e+1]
    element2nodei[e+1] = element2nodei[e] + n;
    
  }// End Function set_nnode_on_element

//****************************************************************************80
//! \brief get_nnode_on_element : Returns the number of nodes that make up the
//!        element.  
//! \details This function will query the element requested by the user for the
//!          number of nodes that define it.  
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] e The element number (indexed from 0)
//! \return nne The number of nodes defining element e 
//****************************************************************************80
  intT get_nnode_on_element(const intT& e)
  {
    //---> Local Variables
    intT is; /*!< Starting index of linked list access */
    intT ie; /*!< Ending index of linked list access */
    intT nne; /*!< Number of nodes per element */
    
    is = element2nodei[e];
    ie = element2nodei[e + 1] - 1;
    
    nne = ie - is + 1;
    return(nne);
    
  } //End Function get_nnode_on_element

//****************************************************************************80
//! \brief add_node_to_element : Inserts a node into an element's connectivity
//!                              list.
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] e The element to which you are adding the node
//! \param[in] n The local node index i.e. you are adding node n of element e
//! \param[in] node The global node index which makes of node n of element e
//****************************************************************************80
  void add_node_to_element(const intT& e, const intT& n, const intT& node)
  {
    //---> Local Variables
    intT i;
    
    //---> First get the data array index from the index array
    i = element2nodei[e] + n;
    element2node[i] = node;
      
  }//End Function add_node_to_element

//****************************************************************************80
//! \brief get_node_on_element : Returns the nth node of the eth element.   
//! \details This routine will return the request global node number of an 
//!          element.  The user will supply an element and local node number 
//!          both indexing from 0 and the function returns the  
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] e The element index
//! \param[in] n The local node number of the element
//! \return node The grid node number requested by user
//****************************************************************************80
  intT get_node_on_element(const intT& e, const intT& n)
  {
    //---> Local Variables
    intT i;
    intT node;
    
    //---> First get the data array index from the index array
    i = element2nodei[e] + n;
    node = element2node[i];
    return(node);
    
  } // End Function get_node_on_element
  
//****************************************************************************80
//! \brief set_nnode_on_bc_face : Tells the class how many nodes make up an
//!        specified bc_face
//! \details The number of nodes is set via the following formula.  
//!   For an bc_face f, \formula{ n = bc_face2nodei[f+1] - bc_face2nodei[f] }. 
//!   Thefore knowing n means we are setting the value of bc_face2nodei[f+1].
//!   Effectively this function establishs how many nodes there are per 
//!   bc_face and the bc_face2nodei linked list indexing array at the 
//!   same time.      
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] f The bc_face number
//! \param[in] n The number of nodes for this bc_face
//****************************************************************************80
  void set_nnode_on_bc_face(const intT& f, const intT& n)
  {
    //---> Add n nodes to bc_face2nodei[f+1]
    bc_face2nodei[f+1] = bc_face2nodei[f] + n;
    
  }// End Function set_nnode_on_bc_face

//****************************************************************************80
//! \brief get_nnode_on_bc_face : Returns the number of nodes that make up the
//!        bc_face.  
//! \details This function will query the bc_face requested by the user for the
//!          number of nodes that define it.  
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] f The bc_face number (indexed from 0)
//! \return nnf The number of nodes defining bc_face f 
//****************************************************************************80
  intT get_nnode_on_bc_face(const intT& f)
  {
    //---> Local Variables
    intT is; /*!< Starting index of linked list access */
    intT ie; /*!< Ending index of linked list access */
    intT nnf; /*!< Number of nodes per bc_face */
    
    is = bc_face2nodei[f];
    ie = bc_face2nodei[f + 1] - 1;
    
    nnf = ie - is + 1;
    return(nnf);
    
  } //End Function get_nnode_on_bc_face

//****************************************************************************80
//! \brief add_node_to_bc_face : Inserts a node into an bc_face's connectivity
//!                              list.
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] f The bc_face to which you are adding the node
//! \param[in] n The local node index i.e. you are adding node n of bc_face f
//! \param[in] node The global node index which makes of node n of bc_face f
//****************************************************************************80
  void add_node_to_bc_face(const intT& f, const intT& n, const intT& node)
  {
    //---> Local Variables
    intT i;
    
    //---> First get the data array index from the index array
    i = bc_face2nodei[f] + n;
    bc_face2node[i] = node;
      
  }//End Function add_node_to_bc_face

//****************************************************************************80
//! \brief get_node_on_bc_face : Returns the nth node of the fth bc_face.   
//! \details This routine will return the request global node number of an 
//!          bc_face.  The user will supply an bc_face and local node number 
//!          both indexing from 0 and the function returns the  
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] f The bc_face index
//! \param[in] n The local node number of the bc_face
//****************************************************************************80
  intT get_node_on_bc_face(const intT& f, const intT& n)
  {
    //---> Local Variables
    intT i;
    
    //---> First get the data array index from the index array
    i = bc_face2nodei[f] + n;
    return(bc_face2node[i]);
    
  } //End Function get_node_on_bc_face 

//****************************************************************************80
//! \brief set_bc_face_id : Sets the id number for a specified boundary face
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] f The boundary face
//! \param[in] id The id number to assign to the face f
//****************************************************************************80
  void set_bc_face_id(const intT& f, const intT& id)
  {
    //---> Set the bc_face_id
    bc_face_id[f] = id;
    
  } // End Function set_bc_face_id

//****************************************************************************80
//! \brief get_bc_face_id : Gets the id number for a specified boundary face
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] f The boundary face
//! \return id
//****************************************************************************80
  intT get_bc_face_id(const intT& f)
  {
    //---> Set the bc_face_id
    intT id = bc_face_id[f] ;
    return(id);
    
  } // End Function get_bc_face_id

//****************************************************************************80
//!
//! \brief set_bcid_type : Sets the boundary condition type for an id'd patch
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] id The boundary id number
//! \param[in] bctype  The boundary condition type number to specified
//****************************************************************************80
  void set_bcid_type(const intT& id, const intT& bctype)
  {
    //---> Set the boundary condition type for this id #
    bcid_type[id] = bctype;
  }// End set_bcid_type
  
//****************************************************************************80
//!
//! \brief get_bcid_type : Gets the boundary condition type for an id'd patch
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] id The id number you want the boundary condition type of 
//****************************************************************************80
  intT get_bcid_type(const intT& id)
  {
    //---> Get the boundary condition for id
    return(bcid_type[id]);
  } // End get_bcid_type 

//****************************************************************************80
//! \brief set_element_type : Sets the element type for each element
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] e The element whose type you want to set
//! \param[in] etype The type of element you want set element e to
//****************************************************************************80
  void set_element_type(const intT& e, const intT& etype)
  {
    
    //---> Set the element type
    element_type[e] = etype;

  } // End Function set_element_type

//****************************************************************************80
//! \brief get_element_type : Gets the element type for each element
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] e The element whose type you want to get
//! \return etype 
//****************************************************************************80
  intT get_element_type(const intT& e)
  {
    
    //---> Get the element type
    return(element_type[e]);
  

  } // End Function get_element_type

//****************************************************************************80
//! \brief elements_surr_node : Gets the elements surrounding a node for all 
//!        nodes.  
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//****************************************************************************80
  void elements_surr_node()
  {
    //---> Local Variables
    intT e; // The element looping index 
    intT n; // The nodes on element looping index  
    intT node; // The node number 
    intT nne; // Number of nodes attached to element    
    intT size_node2element; // Size of the array node2element
    intT index; // An index temp array
  
    //---> Initialize the array 
    for (n = 0; n <= nnode; n++) node2elementi[n] = 0;

    /*---> Loop over elements and count the number of elements attached to a
      node */
    for (e = 0; e < nelement; e++) { // element_loop 
      
      //---> Get number of nodes on this element
      nne = get_nnode_on_element(e);
     
      //---> Loop over nodes on the element
      for(n = 0; n< nne; n++) { // node_loop
	
	//---> Get node number 
	node = get_node_on_element(e, n);

	/*---> Add 1 to the count of number of elements attached to a node
	  for node */
	node2elementi[node + 1] += 1;
	
      } // End node_loop 
      
    } // End element_loop
    
    /*---> Now add up all the element surrounding each node and setup 
      linked list index array */ 
    size_node2element = 0;
    for ( n = 0; n < nnode ; n++) { // Node_loop 
      //---> Size of node2element array computation
      size_node2element += node2elementi[n + 1];
      //---> Setup of node2element linked list index array 
      node2elementi[n + 1] += node2elementi[n];
     
    } // End Node_loop 
   
    //---> Allocate the memory for node2element data array 
    realT mem =  system_module::alloc_mem<intT, intT, realT>(node2element, 
    						      size_node2element);
    grid_mem += mem;
     
    //---> Now for the connectivity node2element
    for ( e = 0; e < nelement; e++ ) { // Element_loop 
       
      //---> Get number of nodes on this element
      nne = get_nnode_on_element(e);
     
      //---> Loop over nodes on the element
      for(n = 0; n < nne; n++) { // node_loop
	
    	//---> Get node number 
    	node = get_node_on_element(e, n);
	
    	//---> Get index out of node2element index array
    	index = node2elementi[node];
    	//----> Assign element e to index position of node2element 
    	node2element[index] = e;
    	/*----> We've added an element to the node2element for node: node...so
    	 add 1 to node2elementi[node] to act as a counter for how many 
    	 elements we've added to the node, this acts as a local index now */
    	node2elementi[node] += 1;
	
      } // End node_loop 
    
    } // End element_loop 
    
    //---> Now fix node2element index array 
    for (n = nnode; n > 0; n--){ //node_loop
       node2elementi[n] = node2elementi[n-1];
      
     } // End node_loop 
     node2elementi[0] = 0;
     
     return;

  }// End Function elements_surr_node
  
//****************************************************************************80
//! \brief get_nelements_on_node : For a speciefied node returns the number of
//!        elements that contain that node
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] n The node number
//! \return ne Number of elements containing node n
//****************************************************************************80
  intT get_nelement_on_node(const intT& n)
  {
    intT ne;
    /*---> Compute number of elements attached a node from node2element 
      linked list index array. */ 
    ne = node2elementi[n + 1] - node2elementi[n];
    //---> Return value to user
    return(ne);
  } // End Function get_nelement_on_node
  
//****************************************************************************80
//! \brief get_element_on_node : For a specified node returns the specified 
//!        element that contains that node.  I.e. the ith element containing 
//!        the node n. 
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] n The node numder
//! \param[in] i The index of the element containing node n we want. 
//! \return e The element number, ith element containing node n
//****************************************************************************80
  intT get_element_on_node(const intT& n, const intT& i)
  {
    //---> Local Variables
    intT e; // Element index
    intT index; // Index of data array 
    
    //---> Get index of ith element containing node n
    index = node2elementi[n] + i;
    
    //---> The ith element containing node n
    e = node2element[index];
    
    //---> Return the ith element containing node n to user
    return(e);
    
  } // End Function get_element_on_node
  
//****************************************************************************80
//! \brief get_node_adj : Computes an adjacency array of nodes NOT ELEMENTS
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//****************************************************************************80
  void get_node_adj()
  {
    intT n; // Node index
    intT ne; // Number of elements containing a node
    intT i; // Looping index
    intT j; // Looping index
    intT nn_tot;// Total number of possible nodes adjacent to node n
    realT mem;
    std::vector<intT> node_adj_temp;
      
    //---> Initialize size of list to zero
    node_adji[0] = 0;
    //---> Loop over the nodes
    for (n = 0; n < nnode; n++) { // Node_loop
      
      //---> Get number of elements surrounding a node
      ne = get_nelement_on_node(n);
    
      //---> Initialize total possible nodes adjacent to node n
      nn_tot = 0;
      
      /*---> Loop over elements attached node n and count how many nodes
	they contain */
      for(i = 0; i < ne; i++ ) {// Element_loop 
	//---> Get element number for ith element surrounding node n
	 intT e = get_element_on_node(n, i);
	 
	//---> Number of nodes on element e
	intT nn = get_nnode_on_element(e);
	
	//---> Update total
	nn_tot += nn;
      }
      
      //---> Temporary array for nodes attached to node n
      std::vector<int> temp (nn_tot,-1);
      
      //---> Initialize counter
      intT counter = 0;
      
      //---> Now loop over the elements containing that node
      for( i = 0; i < ne; i++ ) {// element_loop
	//---> Get element containing node n
	intT e = get_element_on_node(n, i);

	//---> Get number of nodes of element e
  	intT nn = get_nnode_on_element(e);
	
	//---> Loop over nodes attached to elem
  	for ( j = 0; j < nn; j++){ // Node_on_element
  	
	  //---> Get node index of jth node attached to elem
	  intT node = get_node_on_element(e, j);	 
	
	  //---> Store node in temp 
  	  temp[counter] = node;
	  
	  //---> Increment counter
	  counter += 1;
	
  	} // Node_on_element

      } // End element_loop 
      
      //---> Sort all possible nodes adjacent to node n in lexigraphical order
      std::sort (temp.begin(), temp.end());

      //---> Define variables to use unique function 
      std::vector<int>::iterator itemp;
    
      /*---> Eleminate all duplicate nodes from list temp NOTE: Must call 
	sort first for std::unique to work as expected.  */
      itemp = std::unique(temp.begin(), temp.end());
  
      //---> Resize temp to remove undefined entries
      temp.resize(std::distance( temp.begin(),itemp ) );

      //---> Return size of vector temp which is number of unique adjacent nodes
      intT c = temp.size();
 
      for( i = 0; i < c; i++){// Fill_loop 
	//---> Fill adjacency array
	node_adj_temp.push_back(temp[i]);
      }// End Fill_loop 
     
      //---> Update the node adjacency linked list index array to end of node n
      node_adji[n+1] = node_adji[n] + c;
     
    } // End Node_loop 
 
    //--------------------------------------------------------------------------
    mem = system_module::alloc_mem<intT, intT, realT>(node_adj, 
						      node_adj_temp.size());
    grid_mem +=mem;
    //--------------------------------------------------------------------------
    //---> Now copy the contents from node_adj_temp to node adj
    std::copy(node_adj_temp.begin(), node_adj_temp.end(), node_adj);

  }// End get_node_adj

//****************************************************************************80
//! \brief get_nnode_adj : Get the number of nodes in adjacency structure of
//! specified node n . 
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] n The node you want number of adjacent nodes of. 
//! \return nn The number of adjacent nodes to node 
//****************************************************************************80
  intT get_nnode_adj(const intT& n )
  {
    //---> Compute number of adjacent nodes to node
    intT nn = node_adji[n + 1] - node_adji[n];
    return(nn);
    
  } // End get_nnode_adj

//****************************************************************************80
//! \brief get_node_adj : Get the ith adjacent to node to specified node n
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] n The node you want the adjacency structure of.
//! \param[in] i Local index specifying which node adjacency to n you want.
//! \return node The ith node adjacent to node n
//****************************************************************************80
  intT get_node_adj(const intT& n, const intT& i)
  {
    //---> Get the specified node 
    intT index = node_adji[n] + i;
    intT node = node_adj[index];
    
    //---> Return value of node to user
    return(node);
    
  } //End get_node_adj

//****************************************************************************80
//! \brief initialize_bc_faces : This function does operations to initialize 
//!        boundary face connectivity data.  
//! \details This function brings together several aspects of bc faces.
//!     \verbatim
//!        1.  Finds element to which bc faces belong
//!        2.  Determines which face of the element a bc face is
//!        3.  Determines if mesh supplied node ordering of boundary face
//!            is correct.  Correct is defined as normal pointing out of 
//!            element.
//!      \endverbatim
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//****************************************************************************80
  void initialize_bc_faces() 
  {
    //---> Local Variables
    intT f; // Face looping index
    //---> For each boundary face figure out to which element it belongs
    for (f = 0; f < nbc_face; f++) {// BC face loop 
      intT e = find_bc_elem(f);
      bc_face2element[f] = e;
    }// End bc_face_loop 
    
    /*---> Now for each boundary face figure out local face number. */
    for (f = 0; f< nbc_face; f++){// BC face loop 
      
      intT e = bc_face2element[f]; // Element attached to bc_face f
      intT etype =  element_type[e]; // Element type
      intT loc_face; // Local Face number
      intT nne = get_nnode_on_element(e); // Number of nodes on the element e
      intT nnf = get_nnode_on_bc_face(f); // Number of nodes on the face f
           
      //---> Temp vectors 
      std::vector<intT> face_nodes(nnf,-1);
      std::vector<intT> elem_nodes(nne,-1);
         
      //---> Fill elem_nodes
      for (intT i = 0; i < nne; i++){// elem nodes loop
	//---> Access element connecitivity
	elem_nodes[i] = get_node_on_element(e, i);
      }// End elem nodes loop 
      
      //---> Fill face_nodes
      for (intT i = 0; i < nnf; i++){// face nodes loop 
	//---> Acces face connecitivity
	face_nodes[i] = get_node_on_bc_face(f, i);
      }// end face nodes loop 
          
      //---> Based on cell type figure out local face id
      switch (etype) {
      case 0 : // 1-D Bar elements
	loc_face = get_loc_face_1D(face_nodes, elem_nodes);
	break;
      case 1 : // 2-D Triangles
	loc_face = get_loc_face_tri(face_nodes, elem_nodes);
	break;
      case 2 : // 2-D Quadrilaterals
	loc_face = get_loc_face_quad(face_nodes, elem_nodes);
	break;
      case 3 : // 3-D Tetrahedra
	loc_face = get_loc_face_tet(face_nodes, elem_nodes);
	break;
      case 4 : // 3-D Prism 
	break;
      case 5 : // 3-D Pyramid
	break;
      case 6 : // 3-D Hex
	break;
      }
  
      if( loc_face == -99 ) std::cout << "ERROR: Could not find local face index on of face " << f << " For element " << e << std::endl;
      //---> Assign the value of bc_local_face for bc_face f
      bc_local_face[f] = loc_face;
     
      
    }// End BC face loop 

  }// End get_bc_face_to_element

//****************************************************************************80
//! \brief get_element_on_bc_face : Returns the elemnet attached to a
//!   boundary face
//! \details 
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] f The boundary face 
//! \return e The element attached to boundary face f
//****************************************************************************80
  intT get_element_on_bc_face(const intT& f) 
  {
    //---> Get element from already formed list
   intT e = bc_face2element[f];
   //---> Return to user
   return(e);
  }// End get_element_on_bc_80


//****************************************************************************80
//!
//! \brief get_bc_local_face : Gets the local face number for specified bc_face.
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] f The bc_face you want the local face number for 
//! \return bc_loc_face The local face number of bc_face f
//****************************************************************************80
  intT get_bc_local_face(const intT& f)
  {
    //---> Return bc local face of bc_face f to user
    return(bc_local_face[f]);
  } // End get_bc_local_face

//****************************************************************************80
//!
//! \brief get_element_vol : Returns the element volume
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] e The element who's volume you want to get
//! \return element_vol(e) The volume of element e
//****************************************************************************80
  realT get_element_vol(const intT& e)
  {
    return(element_vol(e));
  } // End get_element_vol

//****************************************************************************80
//!
//! \brief set_element_vol : Sets the element volume
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] e The element who's volume you want to set
//! \param[in] vol The volume you want to assign to element e
//****************************************************************************80
  void set_element_vol(const intT& e, const realT& vol)
  {
    element_vol[e] = vol;
  } // End set_element_vol

//****************************************************************************80
//! \brief init_connectivity : Initializes grid connectivity
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//****************************************************************************80
  void init_connectivity()
  {
    //---> Get the elements surrounding a node
    elements_surr_node();
    //---> Get the node adjacency 
    get_node_adj();
    //---> Get boundary face connectivity
    initialize_bc_faces();
    
    //---> Count number of boundary faces for each bc-id number
    count_faces_per_bcid();

  } // End init_connectivity
  
//****************************************************************************80
//!
//! \brief  count_faces_per_bcid : Counts the number faces for a bc id number
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! 
//****************************************************************************80
  void count_faces_per_bcid()
  {
    //---> Initalize
    for(intT id = 0; id < nbc_id; id++){nface_per_bcid[id] = 0;}
    
    //---> Loop over boundary faces

    for(intT bcf = 0; bcf < nbc_face; bcf++) { // bc_face loop  
      intT id = bc_face_id[bcf];
      nface_per_bcid[id] +=1;
    } // End bc_face loop
    
  }// End count_faces_per_bcid


//****************************************************************************80
//!
//! \brief write_tecplot : Writes the mesh to tecplot format
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! 
//****************************************************************************80
  void write_tecplot(std::string& Project)
  {
    std::ofstream tecfile;
    
    tecfile.open( (char*)(Project + std::string(".tec")).c_str() );
    tecfile << "Title = Computational Mesh" << std::endl;
    tecfile << "FILETYPE=GRID" << std::endl;
    switch (ndim){
    case 1:
      tecfile << "Variables = X ID" << std::endl;
      break;
    case 2:
      tecfile << "Variables = X, Y, ID" << std::endl;
      break;
    case 3:
      tecfile << "Variables = X, Y, Z, ID" << std::endl;
      break;
    }
    tecfile << "ZONE T = Mesh" << std::endl;
    switch (ndim) {
    case 1:
      tecfile << "ZONETYPE=FELINESEG" << std::endl;
      break;
    case 2:
      tecfile << "ZONETYPE=FETRIANGLE" << std::endl;
      break;
    case 3:
      tecfile << "ZONETYPE=FETETRAHEDRON" << std::endl;
      break;
    }
    
    tecfile << "NODES = " << nnode << std::endl;
    tecfile << "ELEMENTS=" << nelement << std::endl;
    tecfile << "DATAPACKING=BLOCK" << std::endl;
    
    switch (ndim) {
    case 1:
      tecfile << "VARLOCATION=([1]=NODAL,[2]=CELLCENTERED)" << std::endl;
      break;
    case 2:
      tecfile << "VARLOCATION=([1-2]=NODAL,[3]=CELLCENTERED)" << std::endl;
      break;
    case 3:
      tecfile << "VARLOCATION=([1-3]=NODAL,[4]=CELLCENTERED)" << std::endl;
      break;
    }
    for(intT d = 0; d < ndim; d++){ // Dimension loop 
      for(intT n = 0; n < nnode; n++){// Node loop 
     	tecfile << get_node_coord(n, d) << " ";
      } // End node loop 
      tecfile << std::endl;
    } // End Dimension loop 

    for(intT e = 0; e < nelement; e++){
      tecfile << e << std::endl;
    } // Element ID ID 
    
    for(intT e = 0; e < nelement; e++){// Element loop 
      for(intT n = 0; n < get_nnode_on_element(e); n++) { // Node loop 
	tecfile << get_node_on_element(e,n)+1 << " ";
      } // End Node loop 
      tecfile << std::endl;
    } // End element loop 
    
    //--------------------------- Write Boundary ZONES -------------------------
    for(intT id = 0; id < nbc_id; id++){// id loop 
      std::string idstr;
      std::ostringstream converter;
      converter << id;
      idstr = converter.str();
      tecfile << "ZONE T = \"Boundary-" << idstr << " " 
    	      << get_bcid_tag(id) << "\"" << std::endl;
      
      switch (ndim) {
      case 2: 
    	tecfile << "ZONETYPE=FELINESEG" << std::endl;
    	break;
      case 3:
    	tecfile << "ZONETYPE=FETRIANGLE" << std::endl;
    	break;
      }
      tecfile << "NODES = " << nnode << std::endl;
      tecfile << "ELEMENTS=" << nface_per_bcid[id] << std::endl;
      tecfile << "DATAPACKING=BLOCK" << std::endl;
      
      switch (ndim) {
      case 2:
    	tecfile << "VARLOCATION=([1-2]=NODAL,[3]=CELLCENTERED)" << std::endl;
    	tecfile << "VARSHARELIST=([1-2]=1)" << std::endl;
    	break;
      case 3:
    	tecfile << "VARLOCATION=([1-3]=NODAL,[4]=CELLCENTERED)" << std::endl;
    	tecfile << "VARSHARELIST=([1-3]=1)" << std::endl;
    	break;
      } 
      for( intT bcf = 0; bcf < nbc_face; bcf++ ){ // BC face loop 
    	if( bc_face_id[bcf] == id ){ // Check id
	  tecfile << bcf << " ";
	  tecfile << std::endl;
    	} // End Check id
      } // End BC face loop 

      for( intT bcf = 0; bcf < nbc_face; bcf++ ){ // BC face loop 
    	if( bc_face_id[bcf] == id ){ // Check id
    	  for( int n = 0; n < get_nnode_on_bc_face(bcf); n++ ) {// Node loop 
    	    tecfile << get_node_on_bc_face(bcf, n) + 1 << " ";
    	  } // End Node loop 
	  
    	  tecfile << std::endl;
    	} // End Check id
      } // End BC face loop 
      


    }// End id loop 


    tecfile.close();
        
  }// End write_tecplot

}; // End class UnstGrid
#endif
