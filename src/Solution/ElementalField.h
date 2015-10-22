
#ifndef ELEMENTALFIELD_H_
#define ELEMENTALFIELD_H_
#include "Solution/MeshField"
//****************************************************************************80
//! \class ElementalField
//! \brief A mesh field that is defined on the elements and naturally
//!        discontinuous at element boundaries.
//! \nick
//! \version $Rev$
//! \date $Date$
//!
//****************************************************************************80
class ElementalField : public MeshField
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
  ElementalField(const UnstMesh& mesh,
                 const Array1D<intT> p,
                 const Array1D<intT> nvar);
//****************************************************************************80
//! \brief Destructor
//! \nick
//! \version $Rev$
//! \date $Date$
//****************************************************************************80
  ~ElementalField();


private:
  List2D<realT> data_index_;
  ElementalField() = delete;
};
