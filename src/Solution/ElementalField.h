
#ifndef ELEMENTALFIELD_H_
#define ELEMENTALFIELD_H_
#include "Solution/Field.h"
#include "Mesh/UnstMesh.h"
#include "Elements/Element.h"
//****************************************************************************80
//! \class ElementalField
//! \brief A mesh field that is defined on the elements and naturally
//!        discontinuous at element boundaries.
//! \nick
//! \version $Rev$
//! \date $Date$
//!
//****************************************************************************80
class ElementalField : public Field
{
public:
//****************************************************************************80
//! \brief Constructor
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[in] mesh The mesh over which we define the field
//! \param[in] p The polynomial order of the elements
//! \param[in] nvar The number of variables per degree of freedom.
//****************************************************************************80
  ElementalField(const UnstMesh& mesh, const intT& ndof, const intT& nvar,
                 const List2D<intT>& elem2dof);
//****************************************************************************80
//! \brief Destructor
//! \nick
//! \version $Rev$
//! \date $Date$
//****************************************************************************80
  virtual ~ElementalField() = 0;

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
    intT ndof = elem2dof_.get_ncol(e);
    intT k = 0;
    for(intT i = 0; i < ndof; i++){
      intT dof = elem2dof_(e,i);

      for(intT j = 0; j < nvar_(dof); j++){
        data(k) = Field::data_(dof,j);
        k++;
      }
    }

  }// End ElementData

private:
  Array1D<intT> nvar_;//!< Number of variables per node
  const UnstMesh& mesh_;//!< Reference to mesh
  ElementalField() = delete;
};
#endif /* ELEMENTALFIELD_H_ */
