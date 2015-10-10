// -*-c++-*-
#ifndef EULERPDESET2D_H
#define EULERPDESET2D_H
//****************************************************************************80
//!
//! \brief A class to describe the Euler equations of compressible gas dynamics.
//! \details This class implements the euler equations according to the
//!          the standard PDE intefaces specified in the document 
//!          (put hyperlink here)
//! \nick 
//! \version $Rev$
//! \date $Date$
//! 
//****************************************************************************80
#include "my_incl.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/Array3D.h"
#include "DataStructures/SquareMatrix.h"

class EulerPDESet2D 
{

protected: 
//+++++++++++++++++++++++++++++++ PROTECED STUFF +++++++++++++++++++++++++++++++
  const realT gamma; 
  const realT gm1; 
 
private:
  realT rhoinf;
  realT uinf;
  realT vinf;
  realT pinf;

//****************************************************************************80
//!
//! \brief EulerPDESet2D : Default constructor, assumes 2-D.  Sets all const 
//!  members of class and calls PDESet constructor to ensure interface 
//!  variables are initialized to sensible numbers. 
//! \details Block class user from instantiating without specfying freestream
//! \nick 
//! \version $Rev$
//! \date $Date$
//! 
//****************************************************************************80
  EulerPDESet2D() = delete;
//****************************************************************************80
//!
//! \brief fs_bc : Sets boundary condition to freestream values.
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[out] qb The boundary state vector
//****************************************************************************80
  template<class FTtype>
  void fs_bc(Array1D<FTtype>& qb)
  {
    //---> Set qb to the freestream value
    qb(0) = rhoinf;
    qb(1) = rhoinf*uinf;
    qb(2) = rhoinf*vinf;
    qb(3) = pinf/gm1 + half*rhoinf*(uinf*uinf + vinf*vinf);
    
  } // End fs_bc;

//****************************************************************************80
//!
//! \brief slip_wall_bc : Sets boundary condition for a slip wall.
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] norm The normal vector at the surface
//! \param[in] vg The grid velocity vector
//! \param[in] q The state vector from the discretization
//! \param[out] qb The boundary state vector
//****************************************************************************80
  template<class FTtype>
  void slip_wall_bc(const Array1D<realT>& norm, const Array1D<realT>&vg, 
		    const Array1D<FTtype>& q, Array1D<FTtype>& qb) 
  {
    //---> Set the boundary condition for a slip wall 
    FTtype Nmag = sqrt(norm(0)*norm(0) + norm(1)*norm(1));
    FTtype Nx = norm(0)/Nmag;
    FTtype Ny = norm(1)/Nmag;
  
    //---> Extract Velocities
    FTtype u = q(1)/q(0);
    FTtype v = q(2)/q(0);
       
    FTtype un = u*Nx + v*Ny;
    FTtype vgn = vg(0)*Nx + vg(1)*Ny;
    
    //---> Get the boudary velocity
    FTtype ub = u - (un - vgn)*Nx;
    FTtype vb = v - (un - vgn)*Ny;

    //---> Boundary condition
    qb(0) = q(0);
    qb(1) = ub*qb(0);
    qb(2) = vb*qb(0);;
    qb(3) = q(3) - q(0)*half*(u*u + v*v) + q(0)*half*(ub*ub + vb*vb);

  } // End slip_wall_bc

//****************************************************************************80
//!
//! \brief char_bc : Sets the characteristic boundary condition for in/out flow
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] norm The normal vector at the surface
//! \param[in] vg The grid velocity vector
//! \param[in] q The state vector from the discretization
//! \param[out] qb The boundary state vector
//****************************************************************************80
  template<class FTtype>
  void char_bc(const Array1D<realT>& norm, const Array1D<realT>& vg, 
	       const Array1D<FTtype>& q, Array1D<FTtype>& qb)
  {

    realT one_over_one_minus_gamma = 1.0/(1.0 - gamma);
    
    //---> Local Variables
    realT Nmag = sqrt(norm(0)*norm(0) + norm(1)*norm(1));
    realT Nx = norm(0)/Nmag;
    realT Ny = norm(1)/Nmag;
    
    //---> Extract Velocities
    FTtype rho = q(0);
    FTtype u = q(1)/rho;
    FTtype v = q(2)/rho;
    FTtype E = q(3)/rho;
    
    FTtype p = rho*gm1*(E - half*(u*u + v*v));
    FTtype a = sqrt(gamma*p/rho);
    FTtype ainf = sqrt(gamma*pinf/rhoinf);
    
    FTtype un = u*Nx + v*Ny;
    FTtype uninf = uinf*Nx + vinf*Ny;
    realT vgn = vg(0)*Nx + vg(1)*Ny;
    
    FTtype Mn = (un - vgn)/a;
    
    FTtype Rp = un + 2.0*a/gm1;
    FTtype Rm = uninf - 2.0*ainf/gm1;
    
    FTtype unb = half*(Rp + Rm);
    FTtype ab = gm1/4.0*(Rp - Rm);
    
    FTtype sb;
    FTtype utb;
    
    //---> Supersonic inflow
    if( Mn < -1.0) { fs_bc(qb); }
    
    //---> Supersonic outflow
    else if( Mn > 1.0) { sup_out_bc(q, qb); }
  
    else { //---> Sub-sonic flow
      //---> 
      if( (un - vgn) <= 0.0) { //---> Inflow
	//---> Entropy
	sb = pinf/(pow(rhoinf,gamma));
	//---> Tangential Velocity 
	utb = -uinf*Ny + vinf*Nx;
      }
      else { //---> Outflow
	//---> Entropy
	sb = p/(pow(rho,gamma));
	//---> Tangential Velocity
	utb = -u*Ny + v*Nx;
      }
     
      //---> Boundary density
      FTtype rhob = pow( sb*gamma/(ab*ab), one_over_one_minus_gamma);
      //---> Boundary pressure
      FTtype pb = sb*pow(rhob,gamma);
      //---> Boundary u-velocity;
      FTtype ub = Nx*unb - Ny*utb;
      //---> Boundary v-velocity;
      FTtype vb = Nx*utb + Ny*unb;
     
      qb(0) = rhob;
      qb(1) = rhob*ub;
      qb(2) = rhob*vb;
      qb(3) = pb/gm1 + half*rhob*(ub*ub + vb*vb);

    } // End Check for super-sonic flow

  } // End char_bc

//****************************************************************************80
//!
//! \brief sup_out_bc : Computes the super-sonic outflow bc
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] q The state vector
//! \param[out] qb The boundary state vector
//****************************************************************************80
  template<class FTtype>
  void sup_out_bc(const Array1D<FTtype>& q, Array1D<FTtype>& qb)
  {
    qb(0) = q(0);
    qb(1) = q(1);
    qb(2) = q(2);
    qb(3) = q(3);
  } // End sup_out_bc

public:
  static const intT ndim = 2;
  static const intT nfld = 4;
//****************************************************************************80
//!
//! \brief EulerFluxDotVec : Inviscid flux dotted with a vector
//! \details \f$ \vec{F_{c}} \cdot \vec{n} \f$
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] n The vector to dot the flux with
//! \param[in] vg The grid velocity
//! \param[in] q The q-vector 
//! \param[out] flux The flux vector (passed by reference)
//****************************************************************************80
  template<class FTtype>
  void FluxDotVec(const Array1D<realT>& n, const Array1D<realT>& vg, 
		  const Array1D<FTtype>& q, Array1D<FTtype>& flux) {
    
    //---> Local Variables
    FTtype rho, p, E, u, v, un, dbyrho;
    realT vgn;
    
    //---> Extract Primatives
    rho = q(0);
    dbyrho = 1.0/rho;
    u = q(1)*dbyrho;
    v = q(2)*dbyrho;
    E = q(3)*dbyrho;
    
    //---> dot(u,n)
    un = u*n(0) + v*n(1);
    //---> dot(vg,n)
    vgn = vg(0)*n(0) + vg(1)*n(1);
    
    //---> Compute pressure according to perfect gas equation of state
    p = gm1*rho*(E - half*(u*u + v*v));
    
    //------------------------------- Flux Vector -----------------------------
    //---> Mass Equation
    flux(0) += rho*(un - vgn);
    
    //---> Momentum Equations
    flux(1) += q(1)*(un - vgn) + p*n(0);
    flux(2) += q(2)*(un - vgn) + p*n(1);

    //---> Energy Equation
    flux(3) += (q(3) + p)*un - q(3)*vgn;
    
  } // End EulerFluxDotVec

//****************************************************************************80
//!
//! \brief DerivEulerFluxDotVecTimesX : Inviscid
//!           \f$ \frac{\partial \vec{F_{c}} }{\partial q}\left\{X\right\} \f$
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] n The vector to dot the flux with
//! \param[in] vg The grid velocity
//! \param[in] q The q-vector 
//! \param[in] X The vector we are multiplying by
//! \param[out] dFdqtX The flux vector (passed by reference)
//****************************************************************************80
  template<class FTtype>
  void DerivFluxDotVecTimesX(const Array1D<realT>& n, 
			     const Array1D<realT>& vg, 
			     const Array1D<FTtype>& q,
			     const Array1D<FTtype>& X,
			     Array1D<FTtype>& dFdqtX)
  {
    //---> Local Variables
    FTtype rho, p, E, u, v, un, dbyrho;
    realT vgn;
    Array1D<FTtype> dundq(nfld), dpdq(nfld);
    
    //---> Extract Primatives
    rho = q(0);
    dbyrho = 1.0/rho;
    u = q(1)*dbyrho;
    v = q(2)*dbyrho;
    E = q(3)*dbyrho;
    
    //---> dot(u,n)
    un = u*n(0) + v*n(1);
    //---> dot(vg,n)
    vgn = vg(0)*n(0) + vg(1)*n(1);
    
    //---> Compute pressure according to perfect gas equation of state
    p = gm1*rho*(E - half*(u*u + v*v));
    
    //---> d(un)/d(q);
    dundq(0) = (-q(1)*n(0) - q(2)*n(1))*dbyrho*dbyrho;
    dundq(1) = dbyrho*n(0);
    dundq(2) = dbyrho*n(1);
    dundq(3) = 0.0;

    //---> d(p)/d(q);
    dpdq(0) = half*gm1*(u*u + v*v);
    dpdq(1) = -gm1*u; 
    dpdq(2) = -gm1*v;
    dpdq(3) = gm1;

    //------------------------------- Flux Vector -----------------------------
    //---> Mass Equation
    //flux(0) += rho*(un - vgn);
    dFdqtX(0) += ((un - vgn) + rho*dundq(0))*X(0) +
      rho*dundq(1)*X(1) + 
      rho*dundq(2)*X(2) + 
      rho*dundq(3)*X(3);
    
    //---> Momentum Equations
    //flux(1) += q(1)*(un - vgn) + p*n(0);
    dFdqtX(1) += ( q(1)*dundq(0) + dpdq(0)*n(0) )*X(0) + 
      ( (un - vgn) + q(1)*dundq(1) + dpdq(1)*n(0) )*X(1) + 
      ( q(1)*dundq(2) + dpdq(2)*n(0) )*X(2) + 
      ( q(1)*dundq(3) + dpdq(3)*n(0) )*X(3);
      
    //flux(2) += q(2)*(un - vgn) + p*n(1);
    dFdqtX(2) += ( q(2)*dundq(0) + dpdq(0)*n(1) )*X(0) + 
      ( q(2)*dundq(1) + dpdq(1)*n(1) )*X(1) + 
      ( (un - vgn) + q(2)*dundq(2) + dpdq(2)*n(1))*X(2) + 
      ( q(2)*dundq(3) + dpdq(3)*n(1) )*X(3);
    
    //---> Energy Equation
    //flux(3) += (q(3) + p)*un - q(3)*vgn;
    dFdqtX(3) += ( dpdq(0)*un + (q(3) + p)*dundq(0) )*X(0) + 
      ( dpdq(1)*un + (q(3) + p)*dundq(1) )*X(1) + 
      ( dpdq(2)*un + (q(3) + p)*dundq(2) )*X(2) + 
      ( (1.0 + dpdq(3))*un + (q(3) + p)*dundq(3) - vgn )*X(3);  
    
  } // End EulerFluxDotVec

//****************************************************************************80
//!
//! \brief EigenSystem : Returns the eigen system of the Euler Equations.
//! \details Returns The product of \f$ [T][|\Lambda|][T]^{-1} \f$ for the 
//!          \f$ \vec{F}\cdot \vec{n} \f$
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] n The vector that is dotted in with the flux
//! \param[in] vg The grid velocities, needed for eigenvalues
//! \param[in] q The solution vector for the state at which we are doing this
//! \param[out] EigSys The eigensystem matrix with \f$[|\Lambda|]\f$ 
//****************************************************************************80
  template<class FTtype>
  void EigenSystem(const Array1D<realT>& n, 
		   const Array1D<realT>& vg, 
		   const Array1D<FTtype>& q,
		   SquareMatrix<FTtype>& EigSys)
  {
    
    //---> Local Variables
    realT Nx = n(0);
    realT Ny = n(1);
    realT Nmag = sqrt(Nx*Nx + Ny*Ny);
    FTtype rho = q(0);
    FTtype dbyrho = 1.0/q(0);
    FTtype u = q(1)*dbyrho;
    FTtype v = q(2)*dbyrho;
    FTtype un = u*Nx + v*Ny;
    realT vgn = vg(0)*Nx + vg(1)*Ny;
    FTtype phi = gm1*half*(u*u + v*v);
    FTtype p = gm1*q(3) - rho*phi;
    FTtype c = sqrt(gamma*p/rho);
    //  Helper division variables
    FTtype dbyc2 = 1.0/(c*c);
    realT dbyNmag = 1.0/(Nmag);
    realT dbygm1 = 1.0/gm1;
    //  Eigen values
    FTtype lam1 = abs(un - vgn);
    FTtype lam2 = abs(un - vgn);
    FTtype lam3 = abs(un - vgn - c*Nmag);
    FTtype lam4 = abs(un - vgn + c*Nmag);
    //  Temp matrix for [|Lambda|]*T^-1
    Array2D<FTtype> LamTinv(nfld,nfld);
      
    //---------------- First put [|Lambda|]*[T]^{-1} in LamTinv ----------------
    /*---> NOTE: HEY YOU READING THE CODE!  When reading this translate dby"*"
      as 1.0/() where () is always the name of the local variable you want
      to divide by.  Division is expensive...so I'm trying to avoid it.  
      I'm also trying to make this readable */
    
    //---> Row 0
    LamTinv(0,0) = (1.0 - phi*dbyc2)*lam1; 
    LamTinv(0,1) = (gm1*u*dbyc2)*lam1;
    LamTinv(0,2) = (gm1*v*dbyc2)*lam1;
    LamTinv(0,3) = -(gm1*dbyc2)*lam1;
    
    //---> Row 1
    LamTinv(1,0) = ((Ny*u - Nx*v)*dbyNmag*dbyNmag*dbyrho)*lam2;
    LamTinv(1,1) = (-Ny*dbyNmag*dbyNmag*dbyrho)*lam2;
    LamTinv(1,2) = (Nx*dbyNmag*dbyNmag*dbyrho)*lam2;
    LamTinv(1,3) = 0.0;
    
    //---> Row 2
    LamTinv(2,0) = (half*phi + half*c*un*dbyNmag)*lam3;
    LamTinv(2,1) = (-half*c*Nx*dbyNmag - half*gm1*u)*lam3;
    LamTinv(2,2) = (-half*c*Ny*dbyNmag - half*gm1*v)*lam3;
    LamTinv(2,3) = (half*gm1)*lam3;

    //---> Row 3
    LamTinv(3,0) = (half*phi - half*c*un*dbyNmag)*lam4;
    LamTinv(3,1) = (half*c*Nx*dbyNmag - half*gm1*u)*lam4;
    LamTinv(3,2) = (half*c*Ny*dbyNmag - half*gm1*v)*lam4;
    LamTinv(3,3) = (half*gm1)*lam4;
    
    //------------------------- Now mulptiy by [T] and store -------------------
    
    //---> Row 0
    EigSys(0,0) += LamTinv(0,0) + dbyc2*LamTinv(2,0) + dbyc2*LamTinv(3,0);
    EigSys(0,1) += LamTinv(0,1) + dbyc2*LamTinv(2,1) + dbyc2*LamTinv(3,1);
    EigSys(0,2) += LamTinv(0,2) + dbyc2*LamTinv(2,2) + dbyc2*LamTinv(3,2);
    EigSys(0,3) += LamTinv(0,3) + dbyc2*LamTinv(2,3) + dbyc2*LamTinv(3,3);

    //---> Row 1
    EigSys(1,0) += u*LamTinv(0,0) - rho*Ny*LamTinv(1,0) + 
      (u*dbyc2 - Nx*dbyc2*c*dbyNmag)*LamTinv(2,0) + 
      (u*dbyc2 + Nx*dbyc2*c*dbyNmag)*LamTinv(3,0);
    
    EigSys(1,1) += u*LamTinv(0,1) - rho*Ny*LamTinv(1,1) + 
      (u*dbyc2 - Nx*dbyc2*c*dbyNmag)*LamTinv(2,1) + 
      (u*dbyc2 + Nx*dbyc2*c*dbyNmag)*LamTinv(3,1);
    
    EigSys(1,2) += u*LamTinv(0,2) - rho*Ny*LamTinv(1,2) + 
      (u*dbyc2 - Nx*dbyc2*c*dbyNmag)*LamTinv(2,2) + 
      (u*dbyc2 + Nx*dbyc2*c*dbyNmag)*LamTinv(3,2);
    
    EigSys(1,3) += u*LamTinv(0,3) - rho*Ny*LamTinv(1,3) + 
      (u*dbyc2 - Nx*dbyc2*c*dbyNmag)*LamTinv(2,3) + 
      (u*dbyc2 + Nx*dbyc2*c*dbyNmag)*LamTinv(3,3);

    //---> Row 2
    EigSys(2,0) += v*LamTinv(0,0) + rho*Nx*LamTinv(1,0) + 
      (v*dbyc2 - Ny*dbyc2*c*dbyNmag)*LamTinv(2,0) + 
      (v*dbyc2 + Ny*dbyc2*c*dbyNmag)*LamTinv(3,0);
    
    EigSys(2,1) += v*LamTinv(0,1) + rho*Nx*LamTinv(1,1) + 
      (v*dbyc2 - Ny*dbyc2*c*dbyNmag)*LamTinv(2,1) + 
      (v*dbyc2 + Ny*dbyc2*c*dbyNmag)*LamTinv(3,1);
    
    EigSys(2,2) += v*LamTinv(0,2) + rho*Nx*LamTinv(1,2) + 
      (v*dbyc2 - Ny*dbyc2*c*dbyNmag)*LamTinv(2,2) + 
      (v*dbyc2 + Ny*dbyc2*c*dbyNmag)*LamTinv(3,2);
    
    EigSys(2,3) += v*LamTinv(0,3) + rho*Nx*LamTinv(1,3) + 
      (v*dbyc2 - Ny*dbyc2*c*dbyNmag)*LamTinv(2,3) + 
      (v*dbyc2 + Ny*dbyc2*c*dbyNmag)*LamTinv(3,3);
    
    //---> Row 3
    EigSys(3,0) += phi*dbygm1*LamTinv(0,0) + rho*(Nx*v - Ny*u)*LamTinv(1,0) + 
      (dbygm1 + phi*dbyc2*dbygm1 - un*dbyc2*c*dbyNmag)*LamTinv(2,0) + 
      (dbygm1 + phi*dbyc2*dbygm1 + un*dbyc2*c*dbyNmag)*LamTinv(3,0);
    
    EigSys(3,1) += phi*dbygm1*LamTinv(0,1) + rho*(Nx*v - Ny*u)*LamTinv(1,1) + 
      (dbygm1 + phi*dbyc2*dbygm1 - un*dbyc2*c*dbyNmag)*LamTinv(2,1) + 
      (dbygm1 + phi*dbyc2*dbygm1 + un*dbyc2*c*dbyNmag)*LamTinv(3,1);
 
    EigSys(3,2) += phi*dbygm1*LamTinv(0,2) + rho*(Nx*v - Ny*u)*LamTinv(1,2) + 
      (dbygm1 + phi*dbyc2*dbygm1 - un*dbyc2*c*dbyNmag)*LamTinv(2,2) + 
      (dbygm1 + phi*dbyc2*dbygm1 + un*dbyc2*c*dbyNmag)*LamTinv(3,2);
    
    EigSys(3,3) += phi*dbygm1*LamTinv(0,3) + rho*(Nx*v - Ny*u)*LamTinv(1,3) + 
      (dbygm1 + phi*dbyc2*dbygm1 - un*dbyc2*c*dbyNmag)*LamTinv(2,3) + 
      (dbygm1 + phi*dbyc2*dbygm1 + un*dbyc2*c*dbyNmag)*LamTinv(3,3);
      
  }// End EigenSystem

//****************************************************************************80
//!
//! \brief get_bc_q : Gets the boundary condition solution vector based on 
//!                  bc_type;
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] bc_type The boundary condition type 
//! \param[in] norm The normal vector to the face
//! \param[in] tang The tangent vector to the face
//! \param[in] vg The grid velocity
//! \param[in] q The state vector
//! \param[out] qb The boundary state vector
//****************************************************************************80
  template<class FTtype>
  void get_bc_q(const intT& bc_type, const Array1D<realT>& norm, 
		const Array1D<realT>&, const Array1D<realT>& vg,
		const Array1D<FTtype>& q, Array1D<FTtype>& qb)
  {

    switch(bc_type) {// Select boundary condition
    case 1: //---> Freestream boundary condition
      fs_bc<FTtype>(qb);
      break;
    case 2: //---> Slip-wall boundary condition
      slip_wall_bc<FTtype>(norm, vg, q, qb);
      break;
    case 5: //---> Characteristic boundary condition
      char_bc<FTtype>(norm, vg, q, qb);
      break;
    case 8: //---> Supersonic boundary condition
      sup_out_bc<FTtype>(q, qb);
      break;
    } // End Select boundary condition
    
  }// End get_bc_q;



//****************************************************************************80
//!
//! \brief EulerPDESet : Default constructor, assumes 3-D. 
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] rho The free stream density
//! \param[in] u The free stream u-velocity
//! \param[in] v The free stream v-velocity
//! \param[in] p The free stream pressure
//****************************************************************************80
  EulerPDESet2D(const realT& rho, const realT& u, const realT& v, 
		const realT& p) : gamma(1.4), gm1(.4), rhoinf(rho), uinf(u), 
				  vinf(v), pinf(p)
    {
      //---> All the initialization work is done above
    }

//****************************************************************************80
//!
//! \brief EulerPDESet2D :  Constructor with gamma specified, assumes 2-D.  
//!  Sets all const members of class and calls PDESet constructor to ensure 
//!  interface variables are initialized to sensible numbers. 
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] rho The free stream density
//! \param[in] u The free stream u-velocity
//! \param[in] v The free stream v-velocity
//! \param[in] p The free stream pressure
//! \param[in] gam_in Input value of gamma 
//****************************************************************************80
  EulerPDESet2D(const realT& rho, const realT& u, const realT& v, 
		const realT& p, const realT& gam_in) : 
    gamma(gam_in), gm1(gam_in - 1.0),
    rhoinf(rho), uinf(u), vinf(v), pinf(p) 
  {
    //---> All work done in code proceeding : above.
  }

//****************************************************************************80
//! ~EulerPDESet : The desctructor for the EulerPDESet
//! \brief
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! 
//****************************************************************************80
  ~EulerPDESet2D() { }

//****************************************************************************80
//!
//! \brief get_eigen_values : Returns the eigen values of 
//!      \f$ \vec{F}_{c}\cdot\vec{n} \f$ in a 1-D Array structure. 
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] n The vector we have dotted flux wtih  
//! \param[in] vg The grid velocity
//! \param[in] q The solution vector 
//! \param[in] eigen_vals The eigenvalues 
//****************************************************************************80
  template<class FTtype>
  void get_eigen_values(const Array1D<realT>& n, const Array1D<realT>& vg, 
			const Array1D<FTtype>& q, Array1D<FTtype>& eigen_vals)
  {
    //---> Extract primatives
    //---> Local Variables
    realT Nx = n(0);
    realT Ny = n(1);
    realT Nmag = sqrt(Nx*Nx + Ny*Ny);
    FTtype rho = q(0);
    FTtype dbyrho = 1.0/q(0);
    FTtype u = q(1)*dbyrho;
    FTtype v = q(2)*dbyrho;
    FTtype un = u*Nx + v*Ny;
    realT vgn = vg(0)*Nx + vg(1)*Ny;
    FTtype phi = gm1*half*(u*u + v*v);
    FTtype p = gm1*q(3) - rho*phi;
    FTtype c = sqrt(gamma*p/rho);
  
    //--->Eigen values
    eigen_vals(0) = un - vgn;
    eigen_vals(1) = un - vgn;
    eigen_vals(2) = un - vgn - c*Nmag;
    eigen_vals(3) = un - vgn + c*Nmag;

  } // End get_eigen_values

//****************************************************************************80
//!
//! \brief get_qinf : Return the infinity or farfield condition q-vector
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[out] qinf q-vector for infinity a.k.a farfield condition
//****************************************************************************80
  void get_qinf(Array1D<realT>& qinf)
  {
    //---> Set values of qinf to freestream boundary condition
    fs_bc(qinf);
  }
//****************************************************************************80
//!
//! \brief get_nfld : A public accessor method to obtain variable nfld.
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \return nfld The number of fields or PDEs
//****************************************************************************80
  inline static const intT& get_nfld () 
  { 
    return( EulerPDESet2D::nfld);
  }

//****************************************************************************80
//!
//! \brief get_ndim : A public accessor method to obtain variable ndim.
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \return nfld The number of fields or PDEs
//****************************************************************************80
  inline static const intT& get_ndim() 
  {
    return(EulerPDESet2D::ndim);
  }

}; //End Class EulerPDESet

#endif
