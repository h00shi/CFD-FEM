/*
 * UnstMeshBcFaces.h
 *
 *  Created on: Oct 9, 2015
 *      Author: rabbit
 */
// -*-c++-*-
#ifndef UNSTMESHBCFACES_H
#define UNSTMESHBCFACES_H
#include "DataStructures/List2D.h"
#include "DataStructures/Array1D.h"
#include "IO/UnstMeshReader.h"
#include "Mesh/UnstMeshElements.h"

//****************************************************************************80
//! \class UnstMeshBcFaces
//! \brief Container for boundary face information.
//! \nick
//! \version
//! \date
//****************************************************************************80
class UnstMeshBcFaces {
public:
//****************************************************************************80
//! \brief Constructor for UnstMeshBcFaces
//! \nick
//! \version
//! \date
//! \param[in] mesh_reader The class that reads the mesh
//****************************************************************************80
  UnstMeshBcFaces(UnstMeshReader& mesh_reader);
//****************************************************************************80
//! \brief Constructor for UnstMeshBcFaces including elements
//! \nick
//! \version
//! \date
//! \param[in] mesh_reader The class that reads the mesh
//****************************************************************************80
  UnstMeshBcFaces(UnstMeshReader& mesh_reader, UnstMeshElements& mesh_elements);

//****************************************************************************80
//! \brief Destructor for UnstMeshBcFaces
//! \nick
//! \version
//! \date
//****************************************************************************80
  ~UnstMeshBcFaces();
//****************************************************************************80
//! \brief Forms connections between boundary face and element it's attached to
//!        also finds what local face number of the attached element the
//!        boundary face is.
//! \nick
//! \version
//! \date
//! \param[in] mesh_elements Mesh elements, used to figure out element that
//!            each bc face is attached to.
//****************************************************************************80
 void FormBcFace2Element(UnstMeshElements& mesh_elements);

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
 //! \brief get_nbc_node : Returns the number of nbc_Node elements in the mesh
 //! \details
 //! \nick
 //! \version $Rev$
 //****************************************************************************80
 inline intT get_nbc_node() const {return(nbc_node_);}

 //****************************************************************************80
 //! \brief get_nbc_bar : Returns the number of bc bars in the mesh
 //! \details
 //! \nick
 //! \version $Rev$
 //****************************************************************************80
 inline intT get_nbc_bar() const {return(nbc_bar_);}

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
 //! \brief get_bcface2node : Returns reference the bc_face2node connectivity
 //! \details
 //! \nick
 //! \version $Rev$
 //! \return bc_face2node bc_face to node connectivity
 //****************************************************************************80
 inline const List2D<intT>& get_bc_face2node() const {return bc_face2node_;}

 //****************************************************************************80
 //! \brief get_bcface2element : Returns reference to bc_face to element
 //!  connectivity
 //! \details
 //! \nick
 //! \version $Rev$
 //! \return element2node bc_face to element connectivity
 //****************************************************************************80
 inline const Array1D<intT>& get_bc_face2element() const{
   return bc_face2elem_;
 }
 //****************************************************************************80
 //! \brief get_bc_face_id : Returns reference to bc_face_id
 //! \details
 //! \nick
 //! \version $Rev$
 //! \return bc_face_id The patch ID number of a boundary face
 //****************************************************************************80
   inline const Array1D<intT>& get_bc_face_id() const {return bc_face_id_;}

 //****************************************************************************80
 //! \brief Returns reference to bc_face_type_
 //! \details
 //! \nick
 //! \version $Rev$
 //! \return bcid_type The boundary condition type for each patch ID
 //****************************************************************************80
   inline const Array1D<ElementTopology::face_types>& get_bc_face_type()
       const {return bc_face_type_;}

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
 private:
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
  intT nbc_node_; /*!< Number of 1-D boundary nodes in the mesh
  		    (1-D meshes only) */
  intT nbc_bar_; /*!< Number of 2-D boundary edges in the mesh
  		    (2-D meshes only) */
  intT nbc_tri_; /*!< Number of 3-D boundary tris in the mesh
  		   (3-D meshes only) */
  intT nbc_quad_; /*!< Number fo 3-D boundary quads in the mesh
  		    (3-D meshes only) */
  List2D<intT> bc_face2node_; //!< Nodes on boundary faces
  Array1D<ElementTopology::face_types> bc_face_type_; //!< Type of boundary face
  Array1D<intT> bc_face_id_; //!< Boundary face id number
  Array1D<intT> bc_local_face_; /*!< Local face index for boundary faces */
  Array1D<intT> nface_per_bcid_; /*!< Gives number of faces for per boundary
				    id # */
  Array1D<intT> bc_face2elem_; /*!< Boundary face to element array, for a
				    boundary
				    face this gives the element number attached
				    to it. */
//****************************************************************************80
//! \brief FindBcElement : Find the element containing the specified boundary
//!        face
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[in] f The boundary face
//! \param[in] node2element The node2lement mapping
//! \return e The element containing boundary face f
//****************************************************************************80
  intT FindBcElement(const intT& f, const List2D<intT>& node2elemenet);

//****************************************************************************80
//! \brief A template function to extract the local face number for a face
//!  given the nodes on  the element it it attached to.
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[in] f The boundary face
//! \param[in] node2element The node2lement mapping
//! \return e The element containing boundary face f
//****************************************************************************80
  template<class Topology>
  inline intT ExtractLocalFace(const intT& f,
                               const List2D<intT>& elem2node)
  {
    intT e = bc_face2elem_(f);

    std::vector<intT> face_nodes;
    face_nodes.reserve(bc_face2node_.get_ncol(f));

    //---> Load face nodes into temp vector
    for(intT i = 0; i < bc_face2node_.get_ncol(f); i++){
      face_nodes.push_back(bc_face2node_(f,i));
    }

    std::vector<intT> elem_nodes;
    elem_nodes.reserve(Topology::nNode);

    for(intT i = 0; i < Topology::nNode; i++){
     elem_nodes.push_back(elem2node(e,i));
    }

    return ElementTopology::FindLocalFace<Topology>(face_nodes, elem_nodes);
  } // End UnstMeshBcFaces::ExtractLocalFace

  UnstMeshBcFaces() = delete;
};
#endif
