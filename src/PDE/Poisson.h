// -*-c++-*-
#ifndef POISSON_H
#define POISSON_H
//****************************************************************************80
//! \class Poisson Poisson.h
//! \brief This header file describes the class Poisson, which is defines the
//!        system of equations for a poisson equation. 
//! \details Contains operators for the PDE \f$\nabla \cdot (D \nabla u) \f$
//! \nick 
//! \version $Rev$
//! \date $Date$
//! 
//****************************************************************************80
#include "my_incl.h"
#include "consts.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/Array3D.h"

class Poisson {

private:
//++++++++++++++++++++++++++++++++ PRIVATE STUFF +++++++++++++++++++++++++++++++
  realT D; /*!< Diffusion coefficient */
  
public:
  
//****************************************************************************80
//!
//! \brief Poisson : The constructor for this class
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! 
//****************************************************************************80
  Poisson() {
    //---> Set Diffusion coefficient to 1 by default
    D = 1.0;
  }// End Poisson
  
//****************************************************************************80
//!
//! \brief
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! 
//****************************************************************************80
  ~Poisson() {
    //---> Nothing to destroy yet
  } // End ~Poisson 

//****************************************************************************80
//!
//! \brief cg_vol_res : Computes a continuous Galerkin volume residual for
//!                     Poisson's equation. 
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] ndim The number of physical dimensions
//! \param[in] nfld The number of Poisson equations we are solving
//! \param[in] q The solution 
//! \param[in] dq The solution gradient
//! \param[out] flux The flux flux vector
//! \param[out] src The source term vector 
//****************************************************************************80
  void cg_vol_res(const int& ndim, const int& nfld, const Array1D<realT>&,
	     const Array2D<realT>& dq, Array2D<realT>& flux,
	     Array1D<realT>& src) {

    //---> Assemble the 
    for (intT f = 0; f < nfld; f++){//PDE_loop
      for(intT d = 0; d < ndim; d++){//DIM_loop
	flux(d,f) = D*dq(d,f);
	src(f) = 2.0;
      }// End DIM_loop
    }// End PDE_loop
    
  }// End cg_vol_res 

}; // End class Poisson

#endif
