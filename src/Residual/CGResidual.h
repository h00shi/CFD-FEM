//-*-c++-*-
#ifndef CGRESIDUAL_H
#define CGRESIDUAL_H

#include "my_incl.h"
#include "consts.h"
#include "PDE/Equation.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/Array3D.h"
#include "Mesh/CGMesh.h"
#include "Solution/NodalField.h"
#include "Elements/Element.h"
#include "Elements/BarElement.h"
#include "Elements/TriElement.h"
#include "Elements/TetElement.h"
#include "SystemUtils/SystemModule.h"
#include "Residual/Residual.h"
//****************************************************************************80
//!
//! \brief A class to compute Galerkin Least Squares Residuals
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \tparam intT Integter data type
//! \tparam realT Real data type
//! \tparam PDET Partial differential equation type 
//****************************************************************************80
template<class PDET>
class CGResidual : public Residual{

public: 

//****************************************************************************80
//!
//! \brief CGResidual : Constructor taking PDE, grid and soln references
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] PDE_ref The reference to the PDE class
//! \param[in] grid_ref The reference to the grid class
//! \param[in] soln_ref The reference to the soln class
//****************************************************************************80
  CGResidual(PDET & PDE_ref, CGMesh& grid_ref, intT psoln, intT pmap,
             intT pquad) :
      PDE_(PDE_ref), grid_(grid_ref)
  {

    //---> Instantiate Elements
    BarElement bar;
    bar.initialize(psoln, pmap, pquad);
    TriElement tri;
    tri.initialize(psoln, pmap, pquad);
    TetElement tet;
    tet.initialize(psoln, pmap, pquad);

  } // End CGResidual

~CGResidual(){}

//****************************************************************************80
//!
//! \brief CGResidual : ComputeResidual
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$

//****************************************************************************80
  void ComputeResidual(NodalField& soln_field, NodalField& resid_field)
  {
    Array2D<realT> qhat(3,PDET::nfld_);
    for(intT e = 0; e < grid_.get_MeshElements().get_nelement(); e++){
      ElementTopology::element_types etype;
      etype = grid_.get_MeshElements().get_element_type()(e);

      soln_field.ElementData(e,qhat);



    }



  }// End ComputeResidual

//****************************************************************************80
//! \brief Interface for computing a residual and jacobian matrix
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[out] resid Residual argument
//! \param[out] jac Jacobian matrix
//****************************************************************************80
  template<class MatrixType>
  void ComputeResidualJacobian(MatrixType& jac)
  {

  }

private:
  PDET& PDE_; /*!< Reference to PDE class already instantiated */
  CGMesh& grid_; /*!< Reference to grid class already
        instantiated. */

//****************************************************************************80
//!
//! \brief  Solver : Default constructor
//! \details This is blocked...so you can't screw up
//! \nick
//! \version $Rev$
//! \date $Date$
//!
//****************************************************************************80
  CGResidual() = delete;

};

#endif

