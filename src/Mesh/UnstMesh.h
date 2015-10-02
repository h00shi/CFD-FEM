// -*-c++-*-
#ifndef UNSTMESH_H
#define UNSTMESH_H
//****************************************************************************80
//! \class UnstMesh 
//! \brief This is the header file defining the class UnstMesh
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \tparam intT Template argument meant to mimic integer
//! \tparam realT Template argument meant to mimic real numbers
//****************************************************************************80
#include "my_incl.h"
#include "SystemUtils/SystemModule.h"
#include <vector> 
#include <algorithm>
#include <fstream>
#include <map>
#include "DataStructures/Array2D.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/List2D.h"
#include <stdio.h>
#include "Mesh/UnstMeshElements.h"
#include "Mesh/ElementTopology.h"
class UnstMesh {
  
public:
  //+++++++++++++++++++++++++++++++ PUBLIC STUFF +++++++++++++++++++++++++++++++
//****************************************************************************80
//! \brief The constructor from file string
//! \details  Constructor that creates a UnstMesh object with a file input
//! \nick
//! \param[in] filename - name of file (with extension)
//! \param[in] file_type - The type of file we are reading
//****************************************************************************80
  UnstMesh(std::string const & filename, std::string const & file_type);

//****************************************************************************80
//! \brief The constructor from input stream
//! \details  Constructor that creates a UnstMesh object with an input stream
//! \nick
//! \param[in] mesh_stream - input stream that contains the mesh input
//****************************************************************************80
  UnstMesh(std::istream & mesh_stream);

//****************************************************************************80
//! \brief The destructor
//! \details  Destructor is pure virtual to make this class abstract
//! \nick
//! \version $Rev$
//****************************************************************************80
  ~UnstMesh();
 
//****************************************************************************80
//! \brief get_ndim() : Returns the number of dimensions
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
  inline intT get_ndim()  const {return(ndim_);}

//****************************************************************************80
//! \brief get_ndof() : Returns the number of degrees of freedom
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
 inline intT get_ndof()  const {return(ndof_);}

//****************************************************************************80
//! \brief get_nbc_face() : Returns the number of boundary faces in the mesh
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
  inline intT get_nbc_face()  const {return(nbc_face_);}

//****************************************************************************80
//! \brief get_nbc_face() : Returns the number of boundary IDs in the mesh
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
  inline intT get_nbc_id() const {return(nbc_id_);}

//****************************************************************************80
//! \brief get_nnode() : Returns the number of nodes in the mesh
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
  inline intT get_nnode()  const {return(nnode_);}

//****************************************************************************80
//! \brief get_nbar : Returns the number of bar elements in the mesh
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
  inline intT get_nbar() const {return(nbar_);}

//****************************************************************************80
//! \brief get_ntri : Returns the number of tri elements in the mesh
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
 inline intT get_ntri() const {return(ntri_);}

//****************************************************************************80
//! \brief get_nQuad : Returns the number of quad elements in the mesh
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
 inline intT get_nquad() const {return(nquad_);}

//****************************************************************************80
//! \brief get_ntet : Returns the number of tet elements in the mesh
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
 inline intT get_ntet() const {return(ntet_);}

//****************************************************************************80
//! \brief get_nprism : Returns the number of prism elements in the mesh
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
 inline intT get_nprism() const {return(nprism_);}

//****************************************************************************80
//! \brief get_npyr : Returns the number of pyramid elements in the mesh
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
  inline intT get_npyr() const {return(npyr_);}

//****************************************************************************80
//! \brief get_nhex : Returns the number of hex elements in the mesh
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
 inline intT get_nhex() const {return(nhex_);}

//****************************************************************************80
//! \brief get_nbc_node : Returns the number of nbc_Node elements in the mesh
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
 inline intT get_nbc_node() const {return(nbc_node_);}

//****************************************************************************80
//! \brief get_nbc_edge : Returns the number of bc edges in the mesh
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
 inline intT get_nbc_edge() const {return(nbc_edge_);}

//****************************************************************************80
//! \brief get_nbc_tri : Returns the number of bc tri faces in the mesh
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
  inline intT get_nbc_tri() const {return(nbc_tri_);}

//****************************************************************************80
//! \brief get_nbc_quad : Returns the number of bc quad faces in the mesh
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
  inline intT get_nbc_quad() const {return(nbc_quad_);}

//****************************************************************************80
//! \brief get_nelement : Returns the number of elements in the mesh
//! \details
//! \nick
//! \version $Rev$
//! \return nelement The number of elements in the mesh
//****************************************************************************80
  inline intT get_nelement() const {return(nelement_);}

//****************************************************************************80
//! \brief get_element2node : Returns reference to element to node connectivity
//! \details
//! \nick
//! \version $Rev$
//! \return element2node Element to node connectivity
//****************************************************************************80
  inline const List2D<intT>& get_element2node() const {return element2node_;}

//****************************************************************************80
//! \brief get_bcface2element : Returns reference to bc_face to element
//!  connectivity
//! \details
//! \nick
//! \version $Rev$
//! \return element2node bc_face to element connectivity
//****************************************************************************80
  inline const Array1D<intT>& get_bc_face2elem() const {return bc_face2elem_;}

//****************************************************************************80
//! \brief get_bcface2element : Returns reference to bc_face to element
//!  connectivity
//! \details
//! \nick
//! \version $Rev$
//! \return element2node bc_face to element connectivity
//****************************************************************************80
  inline const Array1D<intT>& get_bc_face2element() const 
  { 
    return bc_face2elem_;
  }

//****************************************************************************80
//! \brief get_bcface2node : Returns reference the bc_face2node connectivity
//! \details
//! \nick
//! \version $Rev$
//! \return bc_face2node bc_face to node connecitivity
//****************************************************************************80
  inline const List2D<intT>& get_bc_face2node() const {return bc_face2node_;}

//****************************************************************************80
//! \brief get_element_type : Returns reference to the element type
//! \details
//! \nick
//! \version $Rev$
//! \return element_type The type of element
//****************************************************************************80
  inline const Array1D<intT>& get_element_type() const {return element_type_;}

//****************************************************************************80
//! \brief get_bc_face_id : Returns reference to bc_face_id
//! \details
//! \nick
//! \version $Rev$
//! \return bc_face_id The patch ID number of a boundary face
//****************************************************************************80
  inline const Array1D<intT>& get_bc_face_id() const {return bc_face_id_;}

//****************************************************************************80
//! \brief get_bcid_type : Returns reference to bcid_type
//! \details
//! \nick
//! \version $Rev$
//! \return bcid_type The boundary condition type for each patch ID
//****************************************************************************80
  inline const Array1D<intT>& get_bcid_type() const {return bcid_type_;}

//****************************************************************************80
//! \brief get_bc_local_face : Returns reference to bc_local_face
//! \details
//! \nick
//! \version $Rev$
//! \return bc_local_face The local face number of each boundary face
//****************************************************************************80
  inline const Array1D<intT>& get_bc_local_face() const {return bc_local_face_;}

//****************************************************************************80
//! \brief get_nface_per_bcid : Returns reference to nface_per_bcid
//! \details
//! \nick
//! \version $Rev$
//! \return nface_per_bcid The number of faces on each bc-ID patch
//****************************************************************************80
  inline const Array1D<intT>& get_nface_per_bcid() const 
  {
    return nface_per_bcid_;
  }
  
//****************************************************************************80
//! \brief get_node2element : Returns reference to node2element
//! \details
//! \nick
//! \version $Rev$
//! \return node2element The elements that surround a node
//****************************************************************************80
  inline const List2D<intT>& get_node2element() const {return node2element_;}

//****************************************************************************80
//! \brief get_adj : Returns reference to adj
//! \details
//! \nick
//! \version $Rev$
//! \return adj The adjacency array
//****************************************************************************80
  inline const List2D<intT>& get_adj() const {return adj_;}

//****************************************************************************80
//! \brief get_x : Returns reference to x array, the grid point coordinates
//! \details
//! \nick
//! \version $Rev$
//! \return x The grid point coordinates
//****************************************************************************80
  inline const Array2D<realT>& get_x() const {return x_;}

//****************************************************************************80
//! \brief get_element_vol : Returns reference to element_vol
//! \details
//! \nick
//! \version $Rev$
//! \return element_vol The volume of each element in the mesh
//****************************************************************************80
  inline const Array1D<realT>& get_element_vol() const {return element_vol_;}

//****************************************************************************80
//! \brief get_element_sa : Returns reference to element_sa
//! \details
//! \nick
//! \version $Rev$
//! \return element_sa The surface area of element in the mesh
//****************************************************************************80
  inline const Array1D<realT>& get_element_sa() const {return element_sa_;}

//****************************************************************************80
//! \brief get_edge2node : Returns reference edge2node
//! \details
//! \nick
//! \version $Rev$
//! \return edge2node_ The nodes that make up an edge
//****************************************************************************80
  inline const Array2D<intT>& get_edge2node() const {return edge2node_;}

//****************************************************************************80
//! \brief get_VTKType : Gets the vtk element type corresponding our element
//!        type
//! \details
//! \nick
//! \version $Rev$
//! \param[in] etype Our element type
//****************************************************************************80
  inline intT get_VTKType(const intT& etype) const {return vtk_type[etype];}

//****************************************************************************80
//! \brief get_VTKFaceType : Gets the vtk face type corresponding to our face
//!        type.  
//! \details
//! \nick
//! \version $Rev$
//! \param[in] ftype Our face type
//****************************************************************************80
  inline intT get_VTKFaceType(const intT& ftype) const 
  {
    return vtk_face_type[ftype];
  }// End get_VTKFaceType

//****************************************************************************80
//! \brief MemoryDiagnostic : Runs a full diagnostic of the memory for the 
//!        class.  Will print all data to standard out 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  void Diagnostic(std::ostream& out_stream);
  /*! ELEMENT TYPES
      0 = 1-D Bar, a line segement with end points
      1 = 2-D Triangle, 3 verticies connected together
      2 = 2-D Quad, 4 verticies connected together
      3 = 3-D Tetrahedra, 4 verticies with 4 triangular faces
      4 = 3-D Prism, 6 vertices with 2 trianglular faces and
      3 quadrilateral faces
      5 = 3-D Pyramid, 5 verticies with 4 triangular faces and
      1 quadrilateral face
      6 = 3-D Hexahedra, 8 verticies with 6 quadrilateral faces
  */

protected: 

  //++++++++++++++++++++++++++++++ PROTECTED STUFF +++++++++++++++++++++++++++++
  //---> Class scalar data
  intT ndim_; /*!< Number of physical dimensions of the mesh, currently only
               value s of 1, 2, or 3 are supported. */
  intT nelement_; /*!< Number of elements in the mesh, elements are cannonical
                shapes such as triangles (2-D) and tetrhera (3-D) */
  intT nnode_; //!< Number of nodes in the mesh
  intT nbc_face_; /*!< Number of boundary faces for this mesh, boundary faces
                   have only 1 mesh element attached to it.
                   The collection of boundary faces defines the boundary of
                   the computational domain */
  intT nbc_id_; /*!< Number of unique boundary id numbers for this mesh.
                 Boundary faces are given an id tag in order to group them
                 into logical sub-collections.  For example all boundary faces
                 on the surface of a wing will be given the same bc_id number
                 because they share boundary conditions and in an engineer
                 sense define a component.*/
  intT ncolor_; /*!< Number of colors used to color mesh for parallel solver */
  intT ndof_; /*!< Number of degrees of freedom associated with the mesh */
  intT nline_; /*!< Number of graph lines in the mesh */
  intT nbar_; //!< Number of 1-D bar elements in the mesh
  intT ntri_; //!< Number of 2-D tri elements in the mesh
  intT nquad_; //!< Number of 2-D quad elements in the mesh
  intT ntet_; //!< Number of 3-D tet elements in the mesh
  intT nprism_; //!< Number of 3-D prism elements in the mesh
  intT npyr_; //!< Number of 3-D pyramid elements in the mesh
  intT nhex_; //!< Number of 3-D hex elements in the mesh
  intT nbc_node_; /*!< Number of 1-D boundary nodes in the mesh 
		    (1-D meshes only) */
  intT nbc_edge_; /*!< Number of 2-D boundary edges in the mesh 
		    (2-D meshes only) */
  intT nbc_tri_; /*!< Number of 3-D boundary tris in the mesh
		   (3-D meshes only) */
  intT nbc_quad_; /*!< Number fo 3-D boundary quads in the mesh 
		    (3-D meshes only) */

  //---> UnstMesh Description Pointers
  List2D<intT> element2node_; /*!< Element to node array, gives indicies of
                                nodes on an element */
  Array1D<intT> bc_face2elem_; /*!< Boundary face to element array, for a
                                   boundary
                                   face this gives the element number attached
                                   to it. */
  List2D<intT> bc_face2node_; /*!< Gives nodes on a boundary face for all
                                boundary faces */
  Array1D<intT> element_type_;   //!< Flag for the type of element>
  Array1D<intT> bc_face_id_; //!< For each boundary face this stores the id
  Array1D<intT> bcid_type_; /*!< For each boundary id (or patch if you prefer)
                             gives the boundary condition type */
  Array1D<intT> bc_local_face_; /*!< Local face index for boundary faces */
  Array1D<intT> nface_per_bcid_; /*!< Gives number of faces for per boundary
                                  id # */
  List2D<intT> node2element_; /*!< Data array for node 2 element linked list */
  List2D<intT> adj_; /*!< Node adjacency linked list data array - uninitialized*/
  Array2D<intT> edge2node_; /*!< Give the two nodes that make up an edge */
  
  //---> UnstMesh real values
  Array2D<realT> x_; /*!< Coordinates, stores all the coordinates depending on
              number of dimensions.
              \verbatim
              1-D : (x)
              2-D : (x,y)
              3-D : (x,y,z)
              \endverbatim
            */

  Array1D<realT> element_vol_; /*!< Notional volume of elements,
                     different units for different numbers of dimensions */
  Array1D<realT> element_sa_; /*!< Notional Surface area of elements */

  //---> Memory diagnostics data
  realT grid_mem_; /*!< Total memory for this instance of class */
  //---> Misc
  std::map<intT,std::string> bcid_tag_; /*!< Boundary id string tag */
  intT vtk_type[7] = {3, 5, 9, 10, 13, 14, 12};
  intT vtk_face_type[4] = {3, 5, 9, 10};
//****************************************************************************80
//! \brief AllocateMemory : Allocates the memory for the mesh
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] nnz_element2node Total number of non-zeros in element2node 
//! \param[in] nnz_bc_face2node Total number of non-zeros in bc_face2node
//****************************************************************************80
  void AllocateMemory(const intT& nnz_element2node, 
		      const intT& nnz_bc_face2node); 
  
//****************************************************************************80
//! \brief FindBcElement : Find the element containing the specified boundary 
//!        face
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] f The boundary face
//! \return e The element containg boundary face f
//****************************************************************************80
  intT FindBcElement(const intT& f);

//****************************************************************************80
//! \brief GetLocalFace1D : Get the loca face number of specified 1D element
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] face_nodes Nodes on the face
//! \param[in] elem_nodes Nodes on the element
//! \return loc_face The local face contains the specified face_nodes of 
//! elem_nodes
//****************************************************************************80
  intT GetLocalFace1D(const Array1D<intT>& face_nodes,
		      const Array1D<intT>& elem_nodes);

//****************************************************************************80
//! \brief GetLocalFaceTri : Get the local face number of specified triangle 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] face_nodes Nodes on the face
//! \param[in] elem_nodes Nodes on the element
//! \return loc_face The local face contains the specified face_nodes of 
//! elem_nodes
//****************************************************************************80
  intT GetLocalFaceTri(const Array1D<intT>& face_nodes,
		       const Array1D<intT>& elem_nodes);

//****************************************************************************80
//! \brief GetLocalFaceQuad : Get the local face number of specified quad
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] face_nodes Nodes on the face
//! \param[in] elem_nodes Nodes on the element
//! \return loc_face The local face contains the specified face_nodes of 
//! elem_nodes
//****************************************************************************80
  intT GetLocalFaceQuad(const Array1D<intT>& face_nodes,
			const Array1D<intT>& elem_nodes);

//****************************************************************************80
//! \brief GetLocalFaceTet : Get the local face number of specified tetrahedra 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] face_nodes Nodes on the face
//! \param[in] elem_nodes Nodes on the element
//! \return loc_face The local face contains the specified face_nodes of 
//! elem_nodes
//****************************************************************************80
  intT GetLocalFaceTet(const Array1D<intT>& face_nodes,
		       const Array1D<intT>& elem_nodes);

//****************************************************************************80
//! \brief GetLocalFaceHex : Get the local face number of specified hexahedra 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] face_nodes Nodes on the face
//! \param[in] elem_nodes Nodes on the element
//! \return loc_face The local face contains the specified face_nodes of 
//! elem_nodes
//****************************************************************************80
  intT GetLocalFaceHex(const Array1D<intT>& face_nodes,
		       const Array1D<intT>& elem_nodes);

//****************************************************************************80
//! \brief GetLocalFacePrism : Get the local face number of specified prism 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] face_nodes Nodes on the face
//! \param[in] elem_nodes Nodes on the element
//! \return loc_face The local face contains the specified face_nodes of 
//! elem_nodes
//****************************************************************************80
  intT GetLocalFacePrism(const Array1D<intT>& face_nodes,
			 const Array1D<intT>& elem_nodes);    

//****************************************************************************80
//! \brief GetLocalFacePyramid : Get the local face number of specified pyramid
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] face_nodes Nodes on the face
//! \param[in] elem_nodes Nodes on the element
//! \return loc_face The local face contains the specified face_nodes of 
//! elem_nodes
//****************************************************************************80
  intT GetLocalFacePyramid(const Array1D<intT>& face_nodes,
			   const Array1D<intT>& elem_nodes);    

//****************************************************************************80
//! \brief FormNode2Element : Forms the node2element connecitivity
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  void FormNode2Element();

//****************************************************************************80
//! \brief FormBcFaceConnectivity : Forms the boundary face connectivity
//! \details The boundary face connectivity is all connectivity for the boundary
//!          faces 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  void FormBcFaceConnectivity();

//****************************************************************************80
//! \brief CountElementTypes : Counts up the various element types
//! \details While we don't really need this in the code it's helpful for output
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  void CountElementTypes();

//****************************************************************************80
//! \brief CountBcFaceTypes : Counts up the various boundary face types 
//! \details Again, we don't need this but it's helpful for output
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  void CountBcFaceTypes();

//****************************************************************************80
//! \brief FormConnectivity : A wrapper to form all connectivity of the mesh
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  void FormConnectivity();

//****************************************************************************80
//! \brief FormAdjacency : Forms the adjacency of the mesh
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  void FormAdjacency();
  
//****************************************************************************80
//! \brief ReadGridFile : Reads a .grid formatted file.  This is an ascii file
//!        of my own design.
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] filename A string containing the filename to read
//****************************************************************************80
  void ReadGridFile(std::string const & filename);

private:

//****************************************************************************80
//! \brief  Is the copy constructor
//! \details  The copy constructor is explicitly blocked.
//! \nick
//! \version
//****************************************************************************80
  UnstMesh(const UnstMesh&) = delete;

//****************************************************************************80
//! \brief operator= : Is the default assignment operator for this class.
//!                    This is blocked to prevent assignment of meshes
//! \details  The assignment operator is explicitly blocked.
//! \nick
//! \version
//****************************************************************************80
  UnstMesh& operator= (const UnstMesh&) = delete;

  
};
#endif
