// -*-c++-*-
#ifndef EULERPDESET3D_H
#define EULERPDESET3D_H
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
class EulerPDESet3D 
{

protected: 
//+++++++++++++++++++++++++++++++ PROTECED STUFF +++++++++++++++++++++++++++++++
  const realT gamma;
  const realT gm1;
private:
  realT rhoinf;
  realT uinf;
  realT vinf;
  realT winf;
  realT pinf;
//****************************************************************************80
//!
//! \brief EulerPDESet3D : Default constructor, assumes 3-D.  Sets all const 
//!  members of class and calls PDESet constructor to ensure interface 
//!  variables are initialized to sensible numbers. 
//! \details Block class user from instantiating without specfying freestream
//! \nick 
//! \version $Rev$
//! \date $Date$
//! 
//****************************************************************************80
  EulerPDESet3D() = delete;

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
    qb(3) = rhoinf*winf;
    qb(4) = pinf/gm1 + half*rhoinf*(uinf*uinf + vinf*vinf + winf*winf);
    
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
    realT Nmag = sqrt(norm(0)*norm(0) + norm(1)*norm(1) + norm(2)*norm(2));
    realT Nx = norm(0)/Nmag;
    realT Ny = norm(1)/Nmag;
    realT Nz = norm(2)/Nmag;
  
    //---> Extract Velocities
    FTtype u = q(1)/q(0);
    FTtype v = q(2)/q(0);
    FTtype w = q(3)/q(0);
       
    FTtype un = u*Nx + v*Ny + w*Nz;
    realT vgn = vg(0)*Nx + vg(1)*Ny + vg(2)*Nz;
    
    //---> Get the boudary velocity
    FTtype ub = u - (un - vgn)*Nx;
    FTtype vb = v - (un - vgn)*Ny;
    FTtype wb = w - (un - vgn)*Nz;
    
    //---> Boundary condition
    qb(0) = q(0);
    qb(1) = ub*qb(0);
    qb(2) = vb*qb(0);
    qb(3) = wb*qb(0);
    qb(4) = q(4) - q(0)*half*(u*u + v*v + w*w) + 
                   q(0)*half*(ub*ub + vb*vb + wb*wb);

  } // End slip_wall_bc

//****************************************************************************80
//!
//! \brief char_bc : Sets the characteristic boundary condition for in/out flow
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] norm The normal vector at the surface
//! \param[in] tang The tangent vector at the surface
//! \param[in] vg The grid velocity vector
//! \param[in] q The state vector from the discretization
//! \param[out] qb The boundary state vector
//****************************************************************************80
  template<class FTtype>
  void char_bc(const Array1D<realT>& norm, const Array1D<realT>& tang, 
	       const Array1D<realT>& vg, 
	       const Array1D<FTtype>& q, 
	       Array1D<FTtype>& qb)
  {

    realT one_over_one_minus_gamma = 1.0/(1.0 - gamma);
    
    //---> Local Variables
    realT Nmag = sqrt(norm(0)*norm(0) + norm(1)*norm(1) + norm(2)*norm(2));
    realT Nx = norm(0)/Nmag;
    realT Ny = norm(1)/Nmag;
    realT Nz = norm(2)/Nmag;
    
    realT Tmag = sqrt(tang(0)*tang(0) + tang(1)*tang(1) + tang(2)*tang(2));
    realT Tx = tang(0)/Tmag;
    realT Ty = tang(1)/Tmag;
    realT Tz = tang(2)/Tmag;
    
    //---> Take cross product cross(norm,tang);

    realT Bx = (Ny*Tz - Nz*Ty); 
    realT By = -(Nx*Tz - Tx*Nz);
    realT Bz = (Nx*Ty - Tx*Ny);
    realT Bmag = sqrt(Bx*Bx + By*By + Bz*Bz);
    Bx /= Bmag;
    By /= Bmag;
    Bz /= Bmag;

    //---> Extract Velocities
    FTtype rho = q(0);
    FTtype u = q(1)/rho;
    FTtype v = q(2)/rho;
    FTtype w = q(3)/rho;
    FTtype E = q(4)/rho;
    
    FTtype p = rho*gm1*(E - half*(u*u + v*v + w*w));
    FTtype a = sqrt(gamma*p/rho);
    FTtype ainf = sqrt(gamma*pinf/rhoinf);
    
    FTtype un = u*Nx + v*Ny + w*Nz;
    realT uninf = uinf*Nx + vinf*Ny + winf*Nz;
    realT vgn = vg(0)*Nx + vg(1)*Ny + vg(2)*Nz;
    
    FTtype Mn = (un - vgn)/a;
    
    FTtype Rp = un + 2.0*a/gm1;
    FTtype Rm = uninf - 2.0*ainf/gm1;
    
    FTtype unb = half*(Rp + Rm);
    FTtype ab = gm1/4.0*(Rp - Rm);
    
    FTtype sb;
    FTtype utb;
    FTtype ubnb;
    
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
	utb = uinf*Tx + vinf*Ty + winf*Tz;
	ubnb = uinf*Bx + vinf*By + winf*Bz;
      }
      else { //---> Outflow
	//---> Entropy
	sb = p/(pow(rho,gamma));
	//---> Tangential Velocity
	utb = u*Tx + v*Ty + w*Tz;
	ubnb = u*Bx + v*By + w*Bz;
      }
      
      //---> Boundary density
      FTtype rhob = pow( sb*gamma/(ab*ab), one_over_one_minus_gamma);
      //---> Boundary pressure
      FTtype pb = sb*pow(rhob,gamma);
      
      //----> Set these by inverting the matrix known as M in the notes
      realT DetM = -Nx*By*Tz + Nx*Ty*Bz - Ty*Bx*Nz - Tx*Ny*Bz + 
	By*Tx*Nz + Bx*Ny*Tz;
            
      FTtype ub = ( (Ty*Bz - By*Tz)*unb + 
		    (By*Nz - Ny*Bz)*utb + 
		    (Ny*Tz - Nz*Ty)*ubnb)/DetM;
     
      FTtype vb = ( (Bx*Tz - Tx*Bz)*unb + 
		    (Nx*Bz - Bx*Nz)*utb + 
		    (Tx*Nz - Nx*Tz)*ubnb)/DetM;
           
      FTtype wb = ( (Tx*By - Ty*Bx)*unb + 
		    (Bx*Ny - Nx*By)*utb + 
		    (Nx*Ty - Tx*Ny)*ubnb)/DetM;
      
      qb(0) = rhob;
      qb(1) = rhob*ub;
      qb(2) = rhob*vb;
      qb(3) = rhob*wb;
      qb(4) = pb/gm1 + half*rhob*(ub*ub + vb*vb + wb*wb);

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
    qb(4) = q(4);

  } // End sup_out_bc

public:
  static const intT ndim = 3;
  static const intT nfld = 5;
//****************************************************************************80
//!
//! \brief FluxDotVec : Inviscid flux dotted with a vector
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
    FTtype rho, p, E, u, v, w, un, dbyrho;
    realT vgn;
    
    //---> Extract Primatives
    rho = q(0);
    dbyrho = 1.0/rho;
    u = q(1)*dbyrho;
    v = q(2)*dbyrho;
    w = q(3)*dbyrho;
    E = q(4)*dbyrho;
    
    //---> dot(u,n)
    un = u*n(0) + v*n(1) + w*n(2);
    //---> dot(vg,n)
    vgn = vg(0)*n(0) + vg(1)*n(1) + vg(2)*n(2);
    
    //---> Compute pressure according to perfect gas equation of state
    p = gm1*rho*(E - half*(u*u + v*v + w*w));
    
    //------------------------------- Flux Vector -----------------------------
    //---> Mass Equation
    flux(0) += rho*(un - vgn);
    
    //---> Momentum Equations
    flux(1) += q(1)*(un - vgn) + p*n(0);
    flux(2) += q(2)*(un - vgn) + p*n(1);
    flux(3) += q(3)*(un - vgn) + p*n(2);
    
    //---> Energy Equation
    flux(4) += (q(4) + p)*un - q(4)*vgn;
    
  } // End FluxDotVec

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
    FTtype rho, p, E, u, v, w, un, dbyrho;
    realT vgn;
    Array1D<FTtype> dundq(nfld), dpdq(nfld);
    
    //---> Extract Primatives
    rho = q(0);
    dbyrho = 1.0/rho;
    u = q(1)*dbyrho;
    v = q(2)*dbyrho;
    w = q(3)*dbyrho;
    E = q(4)*dbyrho;
    
    //---> dot(u,n)
    un = u*n(0) + v*n(1) + w*n(2);
    //---> dot(vg,n)
    vgn = vg(0)*n(0) + vg(1)*n(1) + vg(2)*n(2);
    
    //---> Compute pressure according to perfect gas equation of state
    p = gm1*rho*(E - half*(u*u + v*v + w*w));
    
    //---> d(un)/d(q);
    dundq(0) = (-q(1)*n(0) - q(2)*n(1) - q(3)*n(2))*dbyrho*dbyrho;
    dundq(1) = dbyrho*n(0);
    dundq(2) = dbyrho*n(1);
    dundq(3) = dbyrho*n(2);
    dundq(4) = 0.0;

    //---> d(p)/d(q);
    dpdq(0) = half*gm1*(u*u + v*v + w*w);
    dpdq(1) = -gm1*u; 
    dpdq(2) = -gm1*v;
    dpdq(3) = -gm1*w;
    dpdq(4) = gm1;

    //------------------------------- Flux Vector -----------------------------
    //---> Mass Equation
    //flux(0) += rho*(un - vgn);
    dFdqtX(0) = ( (un - vgn) + rho*dundq(0) )*X(0) +
      rho*dundq(1)*X(1) + 
      rho*dundq(2)*X(2) + 
      rho*dundq(3)*X(3) + 
      rho*dundq(4)*X(4);
    
    //---> Momentum Equations
    //flux(1) += q(1)*(un - vgn) + p*n(0);
    dFdqtX(1) += 
      ( q(1)*dundq(0) + dpdq(0)*n(0) )*X(0) + 
      ( (un - vgn) + q(1)*dundq(1) + dpdq(1)*n(0) )*X(1) + 
      ( q(1)*dundq(2) + dpdq(2)*n(0) )*X(2) + 
      ( q(1)*dundq(3) + dpdq(3)*n(0) )*X(3) +
      ( q(1)*dundq(4) + dpdq(4)*n(0) )*X(4);
    
    //flux(2) += q(2)*(un - vgn) + p*n(1);
    dFdqtX(2) += 
      ( q(2)*dundq(0) + dpdq(0)*n(1) )*X(0) + 
      ( q(2)*dundq(1) + dpdq(1)*n(1) )*X(1) + 
      ( (un - vgn) + q(2)*dundq(2) + dpdq(2)*n(1))*X(2) + 
      ( q(2)*dundq(3) + dpdq(3)*n(1) )*X(3) + 
      ( q(2)*dundq(4) + dpdq(4)*n(1) )*X(4);

    //flux(3) += q(3)*(un - vgn) + p*n(2);
    dFdqtX(3) += 
      ( q(3)*dundq(0) + dpdq(0)*n(2) )*X(0) + 
      ( q(3)*dundq(1) + dpdq(1)*n(2) )*X(1) + 
      ( q(3)*dundq(2) + dpdq(2)*n(2) )*X(2) + 
      ( (un - vgn) + q(3)*dundq(3) + dpdq(3)*n(2) )*X(3) +
      ( q(3)*dundq(4) + dpdq(4)*n(2) )*X(4);
    
     
    //---> Energy Equation
    //flux(4) += (q(4) + p)*un - q(4)*vgn;
    dFdqtX(4) += 
      ( dpdq(0)*un + (q(4) + p)*dundq(0) )*X(0) + 
      ( dpdq(1)*un + (q(4) + p)*dundq(1) )*X(1) + 
      ( dpdq(2)*un + (q(4) + p)*dundq(2) )*X(2) + 
      ( dpdq(3)*un + (q(4) + p)*dundq(3) )*X(3) + 
      ( (1.0 + dpdq(4))*un + (q(4) + p)*dundq(4) - vgn )*X(4);  
    
  } // End 

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
		   SquareMatrix< FTtype>& EigSys)
  {
    
    //---> Local Variables
    realT Nx = n(0);
    realT Ny = n(1);
    realT Nz = n(2);
    realT Nmag = sqrt(Nx*Nx + Ny*Ny + Nz*Nz);
    FTtype rho = q(0);
    FTtype dbyrho = 1.0/q(0);
    FTtype u = q(1)*dbyrho;
    FTtype v = q(2)*dbyrho;
    FTtype w = q(3)*dbyrho;
    FTtype un = u*Nx + v*Ny + w*Nz;
    realT vgn = vg(0)*Nx + vg(1)*Ny + vg(2)*Nz;
    FTtype phi = gm1*half*(u*u + v*v + w*w);
    FTtype p = gm1*q(4) - rho*phi;
    FTtype c = sqrt(gamma*p/rho);
    //  Helper division variables
    FTtype dbyc2 = 1.0/(c*c);
    realT dbyNmag = 1.0/(Nmag);
    realT dbygm1 = 1.0/gm1;
    realT dbyNmag2 = dbyNmag*dbyNmag;
    //  Eigen values
    FTtype lam1 = abs(un - vgn);
    FTtype lam2 = abs(un - vgn);
    FTtype lam3 = abs(un - vgn);
    FTtype lam4 = abs(un - vgn - c*Nmag);
    FTtype lam5 = abs(un - vgn + c*Nmag);
    //  Temp matrix for [|Lambda|]*T^-1
    Array2D<FTtype> LamTinv(nfld,nfld);
        
    //---------------- First put [|Lambda|]*[T]^{-1} in LamTinv ----------------
    /*---> NOTE: HEY YOU READING THE CODE!  When reading this translate dby"*"
      as 1.0/() where () is always the name of the local variable you want
      to divide by.  Division is expensive...so I'm trying to avoid it.  
      I'm also trying to make this readable */

    //---> Row 0
    LamTinv(0,0) = ( (Nx*rho - v*Nz + w*Ny)*dbyrho*dbyNmag2 - 
		     phi*Nx*dbyc2*dbyNmag2 )*lam1; 
    
    LamTinv(0,1) = (gm1*u*Nx*dbyc2*dbyNmag2)*lam1;
    LamTinv(0,2) = (Nz*dbyrho*dbyNmag2 + v*Nx*gm1*dbyc2*dbyNmag2)*lam1;
    LamTinv(0,3) = (-Ny*dbyrho*dbyNmag2 + w*Nx*gm1*dbyc2*dbyNmag2)*lam1;
    LamTinv(0,4) = (-gm1*Nx*dbyc2*dbyNmag2)*lam1;
    
    //---> Row 1
    LamTinv(1,0) = ((Ny*rho + Nz*u - Nx*w)*dbyrho*dbyNmag2 - 
		    phi*Ny*dbyc2*dbyNmag2)*lam2;
    LamTinv(1,1) = (-Nz*dbyrho*dbyNmag2 + Ny*u*gm1*dbyc2*dbyNmag2)*lam2;
    LamTinv(1,2) = (gm1*Ny*v*dbyc2*dbyNmag2)*lam2;
    LamTinv(1,3) = (Nx*dbyrho*dbyNmag2 + Ny*w*gm1*dbyc2*dbyNmag2)*lam2;
    LamTinv(1,4)=  (-gm1*Ny*dbyc2*dbyNmag2)*lam2;
    
    //---> Row 2
    LamTinv(2,0) =((Nz*rho - Ny*u + Nx*v)*dbyrho*dbyNmag2 - 
		   phi*Nz*dbyc2*dbyNmag2)*lam3; 
    LamTinv(2,1) = (Ny*dbyrho*dbyNmag2 + Nz*u*gm1*dbyc2*dbyNmag2)*lam3;
    LamTinv(2,2) = (-Nx*dbyrho*dbyNmag2 + Nz*v*gm1*dbyc2*dbyNmag2)*lam3;
    LamTinv(2,3) = (gm1*Nz*w*dbyc2*dbyNmag2)*lam3;
    LamTinv(2,4) = (-gm1*Nz*dbyc2*dbyNmag2)*lam3;
    
    //---> Row 3
    LamTinv(3,0) = (half*phi + half*c*un*dbyNmag)*lam4;
    LamTinv(3,1) = (-half*c*Nx*dbyNmag - half*gm1*u)*lam4;
    LamTinv(3,2) = (-half*c*Ny*dbyNmag - half*gm1*v)*lam4;
    LamTinv(3,3) = (-half*c*Nz*dbyNmag - half*gm1*w)*lam4;
    LamTinv(3,4) = (half*gm1)*lam4;

    //---> Row 4
    LamTinv(4,0) = (half*phi - half*c*un*dbyNmag)*lam5;
    LamTinv(4,1) = (half*c*Nx*dbyNmag - half*gm1*u)*lam5;
    LamTinv(4,2) = (half*c*Ny*dbyNmag - half*gm1*v)*lam5;
    LamTinv(4,3) = (half*c*Nz*dbyNmag - half*gm1*w)*lam5;
    LamTinv(4,4) = (half*gm1)*lam5;
  
    //------------------------- Now mulptiy by [T] and store -------------------
    
    //---> Row 0
    EigSys(0,0) += LamTinv(0,0)*Nx + LamTinv(1,0)*Ny + LamTinv(2,0)*Nz + 
      dbyc2*LamTinv(3,0) + dbyc2*LamTinv(4,0);
    
    EigSys(0,1) += LamTinv(0,1)*Nx + LamTinv(1,1)*Ny + LamTinv(2,1)*Nz + 
      dbyc2*LamTinv(3,1) + dbyc2*LamTinv(4,1);
    
    EigSys(0,2) += LamTinv(0,2)*Nx + LamTinv(1,2)*Ny + LamTinv(2,2)*Nz + 
      dbyc2*LamTinv(3,2) + dbyc2*LamTinv(4,2);
    
    EigSys(0,3) += LamTinv(0,3)*Nx + LamTinv(1,3)*Ny + LamTinv(2,3)*Nz + 
      dbyc2*LamTinv(3,3) + dbyc2*LamTinv(4,3);
   
    EigSys(0,4) += LamTinv(0,4)*Nx + LamTinv(1,4)*Ny + LamTinv(2,4)*Nz + 
      dbyc2*LamTinv(3,4) + dbyc2*LamTinv(4,4);

    //---> Row 1
    EigSys(1,0) += Nx*u*LamTinv(0,0) + (-Nz*rho + Ny*u)*LamTinv(1,0) + 
      (Ny*rho + Nz*u)*LamTinv(2,0) + 
      (-Nx*c*dbyc2*dbyNmag + u*dbyc2)*LamTinv(3,0) + 
      (Nx*c*dbyc2*dbyNmag + u*dbyc2)*LamTinv(4,0); 
                       
    EigSys(1,1) += Nx*u*LamTinv(0,1) + (-Nz*rho + Ny*u)*LamTinv(1,1) + 
      (Ny*rho + Nz*u)*LamTinv(2,1) + 
      (-Nx*c*dbyc2*dbyNmag + u*dbyc2)*LamTinv(3,1) + 
      (Nx*c*dbyc2*dbyNmag + u*dbyc2)*LamTinv(4,1); 
    
    EigSys(1,2) += Nx*u*LamTinv(0,2) + (-Nz*rho + Ny*u)*LamTinv(1,2) + 
      (Ny*rho + Nz*u)*LamTinv(2,2) + 
      (-Nx*c*dbyc2*dbyNmag + u*dbyc2)*LamTinv(3,2) + 
      (Nx*c*dbyc2*dbyNmag + u*dbyc2)*LamTinv(4,2); 

    EigSys(1,3) += Nx*u*LamTinv(0,3) + (-Nz*rho + Ny*u)*LamTinv(1,3) + 
      (Ny*rho + Nz*u)*LamTinv(2,3) + 
      (-Nx*c*dbyc2*dbyNmag + u*dbyc2)*LamTinv(3,3) + 
      (Nx*c*dbyc2*dbyNmag + u*dbyc2)*LamTinv(4,3); 
   
    EigSys(1,4) += Nx*u*LamTinv(0,4) + (-Nz*rho + Ny*u)*LamTinv(1,4) + 
      (Ny*rho + Nz*u)*LamTinv(2,4) + 
      (-Nx*c*dbyc2*dbyNmag + u*dbyc2)*LamTinv(3,4) + 
      (Nx*c*dbyc2*dbyNmag + u*dbyc2)*LamTinv(4,4); 

    //--- Row 2
    EigSys(2,0) += (Nz*rho + Nx*v)*LamTinv(0,0) + Ny*v*LamTinv(1,0) + 
      (-Nx*rho + Nz*v)*LamTinv(2,0) + 
      (-Ny*c*dbyc2*dbyNmag + v*dbyc2)*LamTinv(3,0) + 
      ( Ny*c*dbyc2*dbyNmag + v*dbyc2)*LamTinv(4,0);
    
    EigSys(2,1) += (Nz*rho + Nx*v)*LamTinv(0,1) + Ny*v*LamTinv(1,1) + 
      (-Nx*rho + Nz*v)*LamTinv(2,1) + 
      (-Ny*c*dbyc2*dbyNmag + v*dbyc2)*LamTinv(3,1) + 
      ( Ny*c*dbyc2*dbyNmag + v*dbyc2)*LamTinv(4,1);

    EigSys(2,2) += (Nz*rho + Nx*v)*LamTinv(0,2) + Ny*v*LamTinv(1,2) + 
      (-Nx*rho + Nz*v)*LamTinv(2,2) + 
      (-Ny*c*dbyc2*dbyNmag + v*dbyc2)*LamTinv(3,2) + 
      ( Ny*c*dbyc2*dbyNmag + v*dbyc2)*LamTinv(4,2);

    EigSys(2,3) += (Nz*rho + Nx*v)*LamTinv(0,3) + Ny*v*LamTinv(1,3) + 
      (-Nx*rho + Nz*v)*LamTinv(2,3) + 
      (-Ny*c*dbyc2*dbyNmag + v*dbyc2)*LamTinv(3,3) + 
      ( Ny*c*dbyc2*dbyNmag + v*dbyc2)*LamTinv(4,3);

    EigSys(2,4) += (Nz*rho + Nx*v)*LamTinv(0,4) + Ny*v*LamTinv(1,4) + 
      (-Nx*rho + Nz*v)*LamTinv(2,4) + 
      (-Ny*c*dbyc2*dbyNmag + v*dbyc2)*LamTinv(3,4) + 
      ( Ny*c*dbyc2*dbyNmag + v*dbyc2)*LamTinv(4,4);
    
    //---> Row 3
    EigSys(3,0) += (-Ny*rho + Nx*w)*LamTinv(0,0) + 
      (Nx*rho + Ny*w)*LamTinv(1,0) + Nz*w*LamTinv(2,0) + 
      (-Nz*c*dbyc2*dbyNmag + w*dbyc2)*LamTinv(3,0) + 
      ( Nz*c*dbyc2*dbyNmag + w*dbyc2)*LamTinv(4,0);
    
    EigSys(3,1) += (-Ny*rho + Nx*w)*LamTinv(0,1) + 
      (Nx*rho + Ny*w)*LamTinv(1,1) + Nz*w*LamTinv(2,1) + 
      (-Nz*c*dbyc2*dbyNmag + w*dbyc2)*LamTinv(3,1) + 
      ( Nz*c*dbyc2*dbyNmag + w*dbyc2)*LamTinv(4,1);
    
    EigSys(3,2) += (-Ny*rho + Nx*w)*LamTinv(0,2) + 
      (Nx*rho + Ny*w)*LamTinv(1,2) + Nz*w*LamTinv(2,2) + 
      (-Nz*c*dbyc2*dbyNmag + w*dbyc2)*LamTinv(3,2) + 
      ( Nz*c*dbyc2*dbyNmag + w*dbyc2)*LamTinv(4,2);
    
    EigSys(3,3) += (-Ny*rho + Nx*w)*LamTinv(0,3) + 
      (Nx*rho + Ny*w)*LamTinv(1,3) + Nz*w*LamTinv(2,3) + 
      (-Nz*c*dbyc2*dbyNmag + w*dbyc2)*LamTinv(3,3) + 
      ( Nz*c*dbyc2*dbyNmag + w*dbyc2)*LamTinv(4,3);

    EigSys(3,4) += (-Ny*rho + Nx*w)*LamTinv(0,4) + 
      (Nx*rho + Ny*w)*LamTinv(1,4) + Nz*w*LamTinv(2,4) + 
      (-Nz*c*dbyc2*dbyNmag + w*dbyc2)*LamTinv(3,4) + 
      ( Nz*c*dbyc2*dbyNmag + w*dbyc2)*LamTinv(4,4);
    
    //---> Row 4
    EigSys(4,0) += (rho*(Nz*v - Ny*w) + Nx*phi*dbygm1)*LamTinv(0,0) + 
      (rho*(Nx*w - Nz*u) + Ny*phi*dbygm1)*LamTinv(1,0) + 
      (rho*(Ny*u - Nx*v) + Nz*phi*dbygm1)*LamTinv(2,0) + 
      (dbygm1 + phi*dbyc2*dbygm1 - un*c*dbyc2*dbyNmag)*LamTinv(3,0) + 
      (dbygm1 + phi*dbyc2*dbygm1 + un*c*dbyc2*dbyNmag)*LamTinv(4,0);
    
    EigSys(4,1) += (rho*(Nz*v - Ny*w) + Nx*phi*dbygm1)*LamTinv(0,1) + 
      (rho*(Nx*w - Nz*u) + Ny*phi*dbygm1)*LamTinv(1,1) + 
      (rho*(Ny*u - Nx*v) + Nz*phi*dbygm1)*LamTinv(2,1) + 
      (dbygm1 + phi*dbyc2*dbygm1 - un*c*dbyc2*dbyNmag)*LamTinv(3,1) + 
      (dbygm1 + phi*dbyc2*dbygm1 + un*c*dbyc2*dbyNmag)*LamTinv(4,1);

    EigSys(4,2) += (rho*(Nz*v - Ny*w) + Nx*phi*dbygm1)*LamTinv(0,2) + 
      (rho*(Nx*w - Nz*u) + Ny*phi*dbygm1)*LamTinv(1,2) + 
      (rho*(Ny*u - Nx*v) + Nz*phi*dbygm1)*LamTinv(2,2) + 
      (dbygm1 + phi*dbyc2*dbygm1 - un*c*dbyc2*dbyNmag)*LamTinv(3,2) + 
      (dbygm1 + phi*dbyc2*dbygm1 + un*c*dbyc2*dbyNmag)*LamTinv(4,2);

    EigSys(4,3) += (rho*(Nz*v - Ny*w) + Nx*phi*dbygm1)*LamTinv(0,3) + 
      (rho*(Nx*w - Nz*u) + Ny*phi*dbygm1)*LamTinv(1,3) + 
      (rho*(Ny*u - Nx*v) + Nz*phi*dbygm1)*LamTinv(2,3) + 
      (dbygm1 + phi*dbyc2*dbygm1 - un*c*dbyc2*dbyNmag)*LamTinv(3,3) + 
      (dbygm1 + phi*dbyc2*dbygm1 + un*c*dbyc2*dbyNmag)*LamTinv(4,3);

    EigSys(4,4) += (rho*(Nz*v - Ny*w) + Nx*phi*dbygm1)*LamTinv(0,4) + 
      (rho*(Nx*w - Nz*u) + Ny*phi*dbygm1)*LamTinv(1,4) + 
      (rho*(Ny*u - Nx*v) + Nz*phi*dbygm1)*LamTinv(2,4) + 
      (dbygm1 + phi*dbyc2*dbygm1 - un*c*dbyc2*dbyNmag)*LamTinv(3,4) + 
      (dbygm1 + phi*dbyc2*dbygm1 + un*c*dbyc2*dbyNmag)*LamTinv(4,4);
    
  }

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
  void get_bc_q(const intT& bc_type, const Array1D<FTtype>& norm, 
		const Array1D<FTtype>& tang, const Array1D<FTtype>& vg, 
		const Array1D<FTtype>& q, Array1D<FTtype>& qb)
  {

    switch(bc_type) {// Select boundary condition
    case 1: //---> Freestream boundary condition
      fs_bc(qb);
      break;
    case 2: //---> Slip-wall boundary condition
      slip_wall_bc(norm, vg, q, qb);
      break;
    case 5: //---> Characteristic boundary condition
      char_bc(norm, tang, vg, q, qb);
      break;
    case 8: //---> Supersonic boundary condition
      sup_out_bc(q, qb);
      break;
    } // End Select boundary condition
    
  }// End get_bc_q;

//****************************************************************************80
//!
//! \brief EulerPDESet3D : Default constructor. 
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] rho The free stream density
//! \param[in] u The free stream u-velocity
//! \param[in] v The free stream v-velocity
//! \param[in] w The free stream w-velocity
//! \param[in] p The free stream pressure
//****************************************************************************80
  EulerPDESet3D(const realT& rho, const realT& u, const realT& v, 
		const realT& w, const realT& p) : 
    gamma(1.4), gm1(.4), rhoinf(rho), uinf(u), vinf(v), winf(w), pinf(p)
  {
    //---> All work done in code proceeding : above. 
  }

//****************************************************************************80
//!
//! \brief EulerPDESet3D : Constructor with gamma specificiation 
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] rho The free stream density
//! \param[in] u The free stream u-velocity
//! \param[in] v The free stream v-velocity
//! \param[in] w The free stream w-velocity
//! \param[in] p The free stream pressure
//! \param[in] gam_in Input value of gamma
//****************************************************************************80
  EulerPDESet3D(const realT& rho, const realT& u, const realT& v, 
		const realT& w, const realT& p, const realT& gam_in) : 
    gamma(gam_in), gm1(gam_in - 1.0), rhoinf(rho), uinf(u), vinf(v), winf(w), 
    pinf(p)
  {
    //---> All work done in code proceeding : above. 
  }

//****************************************************************************80
//! ~EulerPDESet3D : The desctructor for the EulerPDESet3D
//! \brief
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! 
//****************************************************************************80
  ~EulerPDESet3D() { }

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
    return( EulerPDESet3D::nfld);
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
    return(EulerPDESet3D::ndim);
  }

}; //End Class EulerPDESet

#endif
