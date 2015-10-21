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

  
public:
static const intT nfld = 1;
//****************************************************************************80
//!
//! \brief Poisson : The constructor for this class
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! 
//****************************************************************************80
  Poisson()
  {

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
  ~Poisson()
  {
    //---> Nothing to destroy yet
  } // End ~Poisson 

  template<class qtype>
  void CalcAccumulation(const Array1D<qtype>& q,
                        const qtype& Vr,
                        Array1D<qtype>& accum)
  {
    for (intT f = 0; f < Poisson::nfld; f++){// Field_Loop
      accum(f) = q(f)*Vr;
    }
  }

//****************************************************************************80
//!
//! \brief CalcNormFlux : Compute flux in normal direction F(c,D,q,dq)\dotn
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] ndim The number of physical dimensions
//! \param[in] norm The normal direction
//! \param[in]       Convection Speed is not required
//! \param[in] dq The solution gradient
//! \param[out] flux The flux flux vector
//****************************************************************************80
  template<class qtype>
  void CalcNormFlux(const intT ndim,
                    const Array1D<realT>& norm,
                    const Array1D<realT>& c,
                    const Array2D<realT>& D,
                    const Array1D<qtype>& q,
                    const Array2D<qtype>& dq,
                    Array1D<realT>& flux)
  {

      flux(0) = 0.0;

      for(intT id = 0; id < ndim; id++){
        flux(0) += 0.0*q(0)*c(id)*norm(id);
      }

      for(intT id = 0; id < ndim; id++){//Out Dim Loop
        for(intT jd = 0; jd < ndim; jd++){// Inner Dim Loop
          flux(0) -= D(id,jd)*dq(jd,0)*norm(id);
        }// End DIM_loop
      }// End Outer DIMM

  }// End cg_vol_res 

}; // End class Poisson

#endif
