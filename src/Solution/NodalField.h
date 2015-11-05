/*
 * NodalField.h
 *
 *  Created on: Oct 16, 2015
 *      Author: rabbit
 */

#ifndef NODALFIELD_H_
#define NODALFIELD_H_
#include "Mesh/ElementTopology.h"
#include "Mesh/UnstMesh.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "Solution/MeshField.h"
//****************************************************************************80
//! \class NodalField
//! \brief And interpretation of a MeshField as a field over the nodes of the
//!        mesh.  This field is naturally continuous at all element boundaries
//! \nick
//! \version $Rev$
//! \date $Date$
//!
//****************************************************************************80
class NodalField : public MeshField
{
public:
//****************************************************************************80
//! \brief Constructor
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[in] mesh The mesh on which we are defining the field
//! \param[in] nvar The number of variables per node
//****************************************************************************80
  NodalField(const UnstMesh& mesh, const Array1D<intT>& nvar);
//****************************************************************************80
//! \brief Destructor
//! \nick
//! \version $Rev$
//! \date $Date$
//!
//****************************************************************************80
  ~NodalField();
//****************************************************************************80
//! \brief Provides access to the data indexing array.
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[in] n Degree of freedom for which you want data
//! \param[out] data The output as an Array1D.
//****************************************************************************80
  inline const Array1D<intT>& get_DataIndex() const {return data_index_;}
  inline const Array1D<intT>& get_Nvar() const {return nvar_;}
//****************************************************************************80
//! \brief Operator for getting data for node n and variable j
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[in] n Degree of freedom for which you want data
//! \param[out] data The output as an Array1D.
//****************************************************************************80
  inline realT& operator()(const intT& n, const intT& j)
  {
    return MeshField::data_(data_index_(n)+j);
  }
//****************************************************************************80
//! \brief Operator for getting data for node n and variable j const version
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[in] n Degree of freedom for which you want data
//! \param[out] data The output as an Array1D.
//****************************************************************************80
  inline const realT& operator()(const intT& n, const intT& j)const
  {
    return MeshField::data_(data_index_(n)+j);
  }
//****************************************************************************80
//! \brief Provides the data as an Array1D for a node of the mesh;
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[in] n Degree of freedom for which you want data
//! \param[out] data The output as an Array1D.
//****************************************************************************80
  inline void NodeData(const intT& n, Array1D<realT>& data)
  {
  for(intT i = 0; i < nvar_(n); i++){
      data(i) = MeshField::data_(data_index_(n)+i);
    }
  }

//****************************************************************************80
//! \brief Provides the beginning pointer to node data
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[in] n The node that begins the data you want
//! \return Pointer to data beginning at node n
//****************************************************************************80
  inline realT* NodeDataBegin(const intT& n)
  {
    return MeshField::data_.get_ptr(data_index_(n));
  }
//****************************************************************************80
//! \brief Provides the beginning pointer to node data const version.
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[in] n The node that begins the data you want
//! \return Pointer to data beginning at node n
//****************************************************************************80
  inline const realT* NodeDataBegin(const intT& n) const
  {
    return MeshField::data_.get_ptr(data_index_(n));
  }

//****************************************************************************80
//! \brief Provides the data as an Array1D for an element of freedom of
//!        the mesh.
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[in] n Degree of freedom for which you want data
//! \param[out] data The output as an Array1D.
//****************************************************************************80
  inline void ElementData(const intT& e, Array1D<realT>& data)
  {
    intT nnode = mesh_.get_MeshElements().get_element2node().get_ncol(e);

    intT k = 0;
    for(intT i = 0; i < nnode; i++){
      intT node = mesh_.get_MeshElements().get_element2node()(e,i);

      for(intT j = 0; j < nvar_(node); j++){
        data(k) = MeshField::data_(data_index_(node) + j);
        ++k;
      }

    }

  }// End ElementData

//****************************************************************************80
//! \brief Provides the data as an Array2D for an element of freedom of
//!        the mesh.
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[in] n Degree of freedom for which you want data
//! \param[out] data The output as an Array2D.
//****************************************************************************80
  inline void ElementData(const intT& e, Array2D<realT>& data)
  {
    intT nnode = mesh_.get_MeshElements().get_element2node().get_ncol(e);

    for(intT i = 0; i < nnode; i++){
      intT node = mesh_.get_MeshElements().get_element2node()(e,i);

      for(intT j = 0; j < nvar_(node); j++){
        data(i,j) = MeshField::data_(data_index_(node) + j);
      }

    }

  }// End ElementData

private:
  Array1D<intT> data_index_; //!< data_indexing array for NodalData
  Array1D<intT> nvar_;//!< Number of variables per node
  const UnstMesh& mesh_;//!< Reference to mesh
  NodalField() = delete;
};



#endif /* NODALFIELD_H_ */
