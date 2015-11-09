// -*-c++-*-
#ifndef ELEMENT_H
#define ELEMENT_H
//****************************************************************************80
//! 
//! \brief This is the header file of the element base class Element.    
//! \details The idea here is that functions and data memebers that are 
//!  common to all elements live in the base class.  For example all elements 
//!  integrate and they all have basis functions.   
//! \nick 
//!  $Rev: 6 $
//!  $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//****************************************************************************80
#include "my_incl.h"
#include "consts.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/Array3D.h"
#include "DataStructures/Array4D.h"
#include "Surreal/Forward/Surreal.h"
#include "DataStructures/SquareMatrix.h"
#include "PDE/Equation.h"

class Element {
private:

protected:
//+++++++++++++++++++++++++++++++ PROTECTED STUFF ++++++++++++++++++++++++++++++
  /*---> The following variables are intended to be accessed from the derived 
    classes. */
  
  intT p_; /*!< The degree of polynomials to store, 
	    specifies how many polynomials to store */
  intT pmap_; /*!< The degree of the mapping that transforms from standard 
	       to physical space */
  intT deg_; /*!< The degree of polynomial to integrate, 
	      specifies how many points are required for integration  */
  intT ndof_; /*!< The number of degrees of freedom, this standard element has*/
  intT ndof_map_; /*!< The number of degrees of freedom, for mapping */
  intT nqp_; /*!< The number of quadrature points */
  intT nqp_face_; /*< The number of face quadrature points */
  intT ndim_; /*!< The number of dimensions physical dimensions */
  intT ndof_face_; /*!< The number of degrees of freedom on a face */
  Array2D<realT> xiq_; /*!< The quadrature points \f$ \xi_{q} \f$ */
  Array1D<realT> wq_; /*!< The quadrature weights \f$ w_{q} \f$*/
  Array2D<realT> xiq_face_; /*!< Face quadrature points */
  Array1D<realT> wq_face_; /*!< Face quadrature weights */
  Array2D<realT> phi_; /*!< The basis functions evaluated at the quadrature 
			points: 
			\f$ \phi_{i}\left(\xi_{q}\right) \f$  */
  Array3D<realT> dphi_dxi_; /*!< The spatial derivative of the basis functions 
			     evaluated at the quadrature points:
			     \f$ \frac{d \phi}{d \xi}\vert_{\xi_{q}} \f$ */
  Array2D<realT> phi_map_; /*!< The mapping basis functions evaluated at the 
			    quadrature points */
  Array3D<realT> dphi_map_; /*!< The spatial derivative of the mapping basis 
			     functions at the quadrature points */
  
  Array3D<realT> phi_face_; /*!< The basis functions evaluated at the 
			      quadrature 
			      points: 
			      \f$ \phi_{i}\left(\xi_{q}\right) \f$, on faces of 
			      standard element. */
  Array4D<realT> dphi_dxi_face_; /*!< The spatial derivative of the basis 
				  functions evaluated at the quadrature points:
				  \f$ \frac{d \phi}{d \xi}\vert_{\xi_{q}} \f$,
				on faces of standard element. */
  Array3D<realT> phi_map_face_; /*!< The mapping basis functions evaluated 
				 at the quadrature points, on faces of 
			       standard element. */
  Array4D<realT> dphi_map_face_; /*!< The spatial derivative of the mapping 
				  basis functions at the quadrature points on
				faces of standar element. */
  Array3D<realT> dphi_ds_; /*!< Derivative of the dim - 1 parameterization*/ 

  Array2D<realT> face_dof_map_;/*!< Mapping face local dofs to element
				global dofs */
  
public:

//****************************************************************************80
//!
//! \brief Element: The constructor, it's empty ...but we declare anyway
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! 
//****************************************************************************80
  Element(){}

//****************************************************************************80
//!
//! \brief ~Element : The destructor but we declare it anyway
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! 
//****************************************************************************80
  virtual ~Element() {} 
  
//****************************************************************************80
//!
//! \brief ProjectToQP : Takes in coefficients and does a projection to the
//!                  quadrature points.  
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] nfld The number of fields 
//! \param[in] qp The quadrature point to project to 
//! \param[in] qcoeff The coefficients of q
//! \param[out] q The projected value 
//!          \f$ q = \sum_{i} \phi_{i}\left(\xi_{qp}\right) \hat{coeff}_{i} \f$
//****************************************************************************80
  template<class qtype>
  void ProjectToQP(const intT& nfld, const int& qp, 
                   const Array2D<realT>& qcoeff, Array1D<qtype>& q)
  {

    //---> Initialize q array to zero
    q.set_value(0.0);

    for (intT f = 0; f < nfld; f++) { // pde_loop
      for(intT dof = 0; dof < ndof_; dof++) { // mode_loop 
        //---> Do projection as a sum
        q(f) += qcoeff(f, dof)*phi_(qp, dof);

      }// End pde_loop
    }// End mode_loop 

  }// ProjectToQP

  //****************************************************************************80
  //!
  //! \brief ProjectGradToQP : Takes in coefficients and does a projection of the
  //!                          gradients to the quadrature point specified
  //! \details
  //! \nick
  //! \version $Rev: 6 $
  //! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
  //! \param[in] nfld The number of fields
  //! \param[in] qp The quadrature point index
  //! \param[in] qcoeff The coefficients of
  //! \param[out] dq The resulting gradient
  //****************************************************************************80
  template<class qtype>
  void ProjectGradToQP(const intT& nfld, const int& qp, 
                       const Array2D<realT>& qcoeff, Array2D<qtype>& dq)
  {
    //---> Initialize dq array to zero
    dq.set_value(0.0);
    for (intT f = 0; f < nfld; f++) { //pde_loop 
      for(intT d = 0; d < ndim_; d++) {//dimension_loop
        for(intT i = 0; i < ndof_; i++) { //mode_loop
          dq(f,d) += qcoeff(f,i)*dphi_dxi_(qp,d,i);
        } // End dimension_loop
      } // End pde_loop
    } // mode_loop

  }// End ProjectGradToQP

//****************************************************************************80
//!
//! \brief ProjectToFaceQP: Takes in coefficients and does a projection to the
//!                  quadrature points on a face.  
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] side The "side" number of the face, side is same as local face
//! \param[in] nfld The number of fields 
//! \param[in] qp The quadrature point to project to 
//! \param[in] qcoeff The coefficients of q
//! \param[out] q The projected value 
//!              \f$ q = \sum_{i} \phi_{i}\left(\xi_{qp}\right) \hat{coeff}_{i} \f$
//****************************************************************************80
  template<class qtype>
  void ProjectToFaceQP(const intT& side, const intT& nfld, const int& qp, 
		       const Array2D<realT>& qcoeff, Array1D<qtype>& q)
  {
    
    //---> Initialize q array to zero
    q.set_value(0.0);
    for (intT f = 0; f < nfld; f++) { // pde_loop
      for(intT i = 0; i < ndof_; i++) { // mode_loop 
	//---> Do projection as a sum
	q(f) += qcoeff(f,i)*phi_face_(side, qp, i);
	
      }// End pde_loop
    }// End mode_loop 
    
  }// End ProjectToFaceQP

//****************************************************************************80
//!
//! \brief ProjectGradToFaceQP : Takes in coefficients and does a projection 
//!                          of the gradients to the quadrature point specified.
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] side The "side" number of the face, side is same as local face
//! \param[in] nfld The number of fields
//! \param[in] qp The quadrature point index
//! \param[in] qcoeff The coefficients of 
//! \param[out] dq The resulting gradient
//****************************************************************************80
  template<class qtype>
  void ProjectGradToFaceQP(const intT& side, const intT& nfld, const int& qp, 
		       const Array2D<realT>& qcoeff, Array2D<qtype>& dq)
  {
    //---> Initialize dq array to zero
    dq.set_value(0.0);
    for (intT f = 0; f < nfld; f++) { //pde_loop 
      for(intT d = 0; d < ndim_; d++) {//dimension_loop
	for(intT i = 0; i < ndof_; i++) { //mode_loop 
	  dq(f,d) += qcoeff(f,i)*dphi_dxi_face_(side, qp, d, i);
	} // End dimension_loop
      } // End pde_loop
    } // mode_loop
    
  }// End ProjectGradToFaceQP


//****************************************************************************80
//!
//! \brief void ComputeJinv : Computes element mapping data 
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] qp The quadrature point 
//! \param[in] xcoeff Mapping coefficients
//! \param[out] Jinv \f$ \frac{d x}{d \xi}\left(\xi_{qp}\right)  \f$ 
//! \param[out] DetJ \f$ Det\left(\frac{d x}{d \xi}\left(\xi_{qp}\right)\right) \f$
//****************************************************************************80
  void ComputeJinv(const intT& qp, const Array2D<realT>& xcoeff, 
                   Array2D<realT>& Jinv, realT& DetJ)
  {
    //---> Local Variables
    Array2D<realT> J;
    //---> Allocate and initialize J matrix
    J.initialize(ndim_, ndim_);
    J.set_value(0.0);

    for (intT dr = 0; dr < ndim_; dr++) { //pde_loop 
      for(intT dc = 0; dc < ndim_; dc++) {//dimension_loop
        for(intT i = 0; i < ndof_map_; i++) { //mode_loop
          J(dr,dc) += xcoeff(dc,i)*dphi_map_(qp,dr,i);
        } // End dimension_loop
      } // End pde_loop
    } // mode_loop

    switch (ndim_) {// Determinant computation
      case 1: //---> 1-D
        DetJ = J(0,0);
        Jinv(0,0) = 1.0/J(0,0);
        break;
      case 2: //---> 2-D
        DetJ = J(0,0)*J(1,1) - J(1,0)*J(0,1);
        //---> Row 0
        Jinv(0,0) = J(1,1)/DetJ;
        Jinv(0,1) = -J(0,1)/DetJ;

        //---> Row 1
        Jinv(1,0) = -J(1,0)/DetJ;
        Jinv(1,1) = J(0,0)/DetJ;
        break;
      case 3:
        DetJ =
            J(0,0)*(J(1,1)*J(2,2) - J(2,1)*J(1,2)) -
            J(0,1)*(J(1,0)*J(2,2) - J(1,2)*J(2,0)) +
            J(0,2)*(J(1,0)*J(2,1) - J(1,1)*J(2,0));

        //---> Row 0
        Jinv(0,0) = (J(1,1)*J(2,2) - J(2,1)*J(1,2))/DetJ;
        Jinv(0,1) = (J(0,2)*J(2,1) - J(2,2)*J(0,1))/DetJ;
        Jinv(0,2) = (J(0,1)*J(1,2) - J(1,1)*J(0,2))/DetJ;

        //---> Row 1
        Jinv(1,0) = (J(1,2)*J(2,0) - J(2,2)*J(1,0))/DetJ;
        Jinv(1,1) = (J(0,0)*J(2,2) - J(2,0)*J(0,2))/DetJ;
        Jinv(1,2) = (J(0,2)*J(1,0) - J(1,2)*J(0,0))/DetJ;

        //---> Row 2
        Jinv(2,0) = (J(1,0)*J(2,1) - J(2,0)*J(1,1))/DetJ;
        Jinv(2,1) = (J(0,1)*J(2,0) - J(2,1)*J(0,0))/DetJ;
        Jinv(2,2) = (J(0,0)*J(1,1) - J(1,0)*J(0,1))/DetJ;

        break;
      default :
        DetJ = -1.0;
        break;
    }// End Determinant computation

    return;
  }// End ComputeJinv

//****************************************************************************80
//!
//! \brief void ComputeJinvFace : Computes element mapping data at a face 
//!        quadrature point.  
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] side The "side" of the element, side is same as local face
//! \param[in] qp The quadrature point 
//! \param[in] xcoeff Mapping coefficients
//! \param[out] Jinv \f$ \frac{d x}{d \xi}\left(\xi_{qp}\right)  \f$ 
//! \param[out] DetJ \f$ Det\left(\frac{d x}{d \xi}\left(\xi_{qp}\right)\right) \f$
//****************************************************************************80
  void ComputeJinvFace(const int& side, const intT& qp, 
		     const Array2D<realT>& xcoeff, 
		     Array2D<realT>& Jinv, realT& DetJ)
  {
    //---> Local Variables
    Array2D<realT> J;
    //---> Allocate and initialize J matrix
    J.initialize(ndim_, ndim_);
    J.set_value(0.0);
    
    for (intT dr = 0; dr < ndim_; dr++) { //pde_loop 
      for(intT dc = 0; dc < ndim_; dc++) {//dimension_loop
	for(intT i = 0; i < ndof_map_; i++) { //mode_loop 
	  J(dr,dc) += xcoeff(dc,i)*dphi_map_face_(side, qp, dr, i);
	} // End dimension_loop
      } // End pde_loop
    } // mode_loop
 
    switch (ndim_) {// Determinant computation
    case 1: //---> 1-D
      DetJ = J(0,0);
      Jinv(0,0) = 1.0/J(0,0);
      break;
    case 2: //---> 2-D
      DetJ = J(0,0)*J(1,1) - J(1,0)*J(0,1);
      //---> Row 0
      Jinv(0,0) = J(1,1)/DetJ;
      Jinv(0,1) = -J(0,1)/DetJ;

      //---> Row 1
      Jinv(1,0) = -J(1,0)/DetJ;
      Jinv(1,1) = J(0,0)/DetJ;
      break;
    case 3:
      DetJ = J(0,0)*(J(1,1)*J(2,2) - J(2,1)*J(1,2)) - 
	J(0,1)*(J(1,0)*J(2,2) - J(1,2)*J(2,0)) + 
	J(0,2)*(J(1,0)*J(2,1) - J(1,1)*J(2,0));
            
      //---> Row 0 
      Jinv(0,0) = (J(1,1)*J(2,2) - J(2,1)*J(1,2))/DetJ;
      Jinv(0,1) = (J(0,2)*J(2,1) - J(2,2)*J(0,1))/DetJ;
      Jinv(0,2) = (J(0,1)*J(1,2) - J(1,1)*J(0,2))/DetJ;
      
      //---> Row 1
      Jinv(1,0) = (J(1,2)*J(2,0) - J(2,2)*J(1,0))/DetJ;
      Jinv(1,1) = (J(0,0)*J(2,2) - J(2,0)*J(0,2))/DetJ;
      Jinv(1,2) = (J(0,2)*J(1,0) - J(1,2)*J(0,0))/DetJ;
      
      //---> Row 2
      Jinv(2,0) = (J(1,0)*J(2,1) - J(2,0)*J(1,1))/DetJ;
      Jinv(2,1) = (J(0,1)*J(2,0) - J(2,1)*J(0,0))/DetJ;
      Jinv(2,2) = (J(0,0)*J(1,1) - J(1,0)*J(0,1))/DetJ;
      
      break;
    default :
      DetJ = -1.0;
      break;
    }// End Determinant computation

    return;
  }// End ComputeJinvFace

//****************************************************************************80
//!
//! \brief ComputeFaceVectors : Computes the mapping of a face at a face 
//!        quadrature point. 
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] side The "side" a.k.a local face number for this
//! \param[in] qp The quadrature point 
//! \param[in] xcoeff Mapping coefficients of an element
//! \param[out] norm Face normal vector
//! \param[out] tang Face tangent vector 
//****************************************************************************80
  void ComputeFaceVectors(const intT& side, const intT& qp, const 
                          Array2D<realT>& xcoeff, Array1D<realT>& norm,
                          Array1D<realT>& tang)
  {
    Array2D<realT> dxds(ndim_, ndim_ - 1);
    dxds.set_value(0.0);
    for (intT dr = 0; dr < ndim_; dr++) { //dim_loop 
      for(intT dc = 0; dc < ndim_ - 1; dc++) {//dimension_loop
        for(intT i = 0; i < ndof_face_; i++) { //mode_loop
          dxds(dr, dc) += xcoeff(dr, face_dof_map_(side, i) ) *
              dphi_ds_(qp, dc, i);
        } // End dimension_loop
      } // End pde_loop
    } // mode_loop

    switch(ndim_){ // normal, tangent
      case 1:
        if( side == 0 ) { norm(0) = -1.0; }
        else if( side == 1 ){ norm(0) = 1.0; }
        tang(0) = 0.0;
        break;
      case 2:
        //---> Tangent vector
        tang(0) = dxds(0,0);
        tang(1) = dxds(1,0);
        //---> Normal vector
        norm(0) = tang(1);
        norm(1) = -tang(0);
        break;
      case 3:
        //---> Tangent vector
        tang(0) = dxds(0,0);
        tang(1) = dxds(1,0);
        tang(2) = dxds(2,0);

        norm(0) =  (dxds(1,0)*dxds(2,1) - dxds(2,0)*dxds(1,1));
        norm(1) = -(dxds(0,0)*dxds(2,1) - dxds(2,0)*dxds(0,1));
        norm(2) =  (dxds(0,0)*dxds(1,1) - dxds(1,0)*dxds(0,1));
        break;
    }// End normal, tangent

  }// End ComputeFaceVectors

//****************************************************************************80
//!
//! \brief IntegrateFunction : Integrates and abitrary function func(q)
//! \details  This function integrates some function func over the element. 
//!           From an objecti oriented point of view the element owns the 
//!           integrate but not neccesarily the function, hence we bring the
//!           function in with a function pointer.  I
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] nfld The number of fields 
//! \param[in] qcoeff The coefficients of q
//! \param[in] xcoeff The coefficients of the mapping \f$ x(\xi) \f$
//! \param[in] func The function pointer of the function being integrated
//! \param[out] integral The resulting of integrating the function i.e 
//! \f$ \int_{-1}^{1} func(q(\xi)) Det(J) d\xi \f$ 
//****************************************************************************80
  void IntegrateFunction(const intT& nfld,
                      const Array2D<realT>& qcoeff,
                      const Array2D<realT>& xcoeff,
                      void(*func)(const Array1D<realT>&, Array1D<realT>& ),
                      Array1D<realT>& integral )
  {
    //---> Local Variables
    Array1D<realT> q(nfld);
    Array2D<realT> Jinv(ndim_,ndim_);
    realT DetJ;
    Array1D<realT> res(nfld);
    
    //---> Set res and integral to zero
    integral.set_value(0);

    for (intT j = 0; j < nqp_ ; j++){// quad_loop 
      //---> Obtain function argument from coefficients and basis at quad_points
      ProjectToQP(nfld, j, qcoeff, q);

      //---> Evaluate function at quadrature point j
      func(q, res);

      //Evaluate mapping at the quadrature point
      ComputeJinv(j, xcoeff, Jinv, DetJ);

      //---> Loop over the fields and integrate the function for each field
      for (intT f = 0; f < nfld; f++){// field_loop 
        integral(f) += res(f)*wq_(j)*DetJ;
      } // End field_loop
    } // End quad_loop 
    return;
  } // End IntegrateFunction

//****************************************************************************80
//!
//! \brief CGVolResid : Computes the continuous Galerkin Volume residual
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] PDE The partial differential equation class you are solving
//! \param[in] qcoeff The state vector coefficients 
//! \param[in] xcoeff The element geometry coefficients
//! \param[out] resid The integrated residual per degree of freedom
//****************************************************************************80
  template < template<class, class>  class PDE_type >
  void CGVolResid(PDE_type<intT,realT>& PDE, 
		  const Array2D<realT>& qcoeff, 
		  const Array2D<realT>& xcoeff, Array2D<realT>& resid) {
    
    // PDE_type<intT,realT> PDE;
    intT nf = PDE.get_nfld();
    Array1D<realT> q(nf);
    Array1D<realT> flux(nf);
    Array1D<realT> norm(ndim_);
    Array1D<realT> vg(ndim_);
    Array1D<realT> dq(nf);
    Array1D<realT> divF(nf);
    Array2D<realT> dqdxi(nf,ndim_);
    Array2D<realT> Jinv(ndim_,ndim_);
    SquareMatrix< realT> tauinv(nf);
    realT DetJ;
    
    //---> Set integral Value to zero
    resid.set_value(0.0);
    vg.set_value(0.0);
          
    for (intT qp = 0; qp < nqp_ ; qp++){// quad_loop
      //---> Obtain function argument from coefficients and basis at quad_points
      ProjectToQP<realT>(nf, qp, qcoeff, q);
   
      //---> Evaluate mapping at the quadrature point
      ComputeJinv(qp, xcoeff, Jinv, DetJ);
      
      //---> Compute dq/d(xi) at a the quadrature point 
      ProjectGradToQP<realT>(nf, qp, qcoeff, dqdxi);

      //---> Set Divergence of F to zero
      divF.set_value(0.0);

      //---------------------- (-d(phi)/d(xi_{d})*F_{d}) -----------------------
      
      for (intT d = 0; d < ndim_; d++){// Dimension loop
	
	//---> Put the col of Jinv into normal vector
	for (intT d1 = 0; d1 < ndim_; d1++){ // Dimension loop-2
	  norm(d1) = Jinv(d1,d);
	}// End Dimension loop-2
	
	//---> Set flux to zero
	flux.set_value(0.0);
		
	//---> Get residual for PDE
	PDE.template FluxDotVec<realT>(norm, vg, q, flux);

	//---> Assemble - d(dphi)/d(xi_{d})*F_{d}
	for (intT f = 0; f < nf; f++) {// Field loop 
	  for(intT dof = 0; dof < ndof_; dof++){// DoF loop 
	    resid(f,dof) -= (dphi_dxi_(qp,d,dof)*flux(f))*DetJ*wq_(qp);
	  } // End Dof loop 
	} // End Field loop 
	
	for (intT f = 0; f < nf; f++) {// Field loop  
	  dq(f) = dqdxi(f,d);
	}
	
	//---> Get the divergence term 
	PDE.template DerivFluxDotVecTimesX<realT>(norm, vg, q, dq, divF);
	
      }// End Dimension loop
          
      //------------------------ [tau]*(Rt(q,dq)) -----------------------------
      
      //---> Form [tau]^{-1}
      tauinv.set_value(0.0);
      
      for( intT dof = 0; dof < ndof_; dof++){ // DoF loop 
	//---> Form put grad(phi_dof) into norm by the following loop
	for( intT d = 0; d < ndim_; d++) {
	  norm(d) = 0.0;
	  for( intT d1 = 0; d1 < ndim_; d1++){
	    norm(d) += Jinv(d,d1)*dphi_dxi_(qp, d1, dof);
	  }
	}

	PDE.template EigenSystem<realT>(norm, vg, q, tauinv);
	
      }// End DoF loop 
      
      flux.set_value(0.0);
      //---> [tau^{-1}]^{-1}*divF, store in flux;
      tauinv.squareSolve(divF, flux);
           
      //---------------- (d(phi)/d(xi_{d})*d(F_{d})/d(q)*tau*Rt(q,dq)) ---------
   
      for (intT d = 0; d < ndim_; d++){// Dimension loop-2
	//---> Put the col of Jinv into normal vector
	for (intT d1 = 0; d1 < ndim_; d1++){ // Dimension loop
	  norm(d1) = Jinv(d1,d);
	}// End Dimension loop

	divF.set_value(0.0);
	PDE.template DerivFluxDotVecTimesX<realT>(norm, vg, q, flux, divF);
	//---> Assemble  d(dphi)/d(xi_{d})*d(F_{d})(dq)*[tau]*R(q,dqdx);
	for (intT f = 0; f < nf; f++) {// Field loop 
	  for(intT dof = 0; dof < ndof_; dof++) { // Dof loop 
	    resid(f,dof) += (dphi_dxi_(qp,d,dof)*divF(f))*DetJ*wq_(qp);
	    
	  }// End Field loop 
	}// End DoF loop 
	
      } // End Dimension loop-2 
      
    }// End Quadrature loop
    return;
  } // End CGVolResid

//****************************************************************************80
//!
//! \brief CGBCResid : Computes the boundary residual for Continuous Galerkin
//!        methods.  
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] PDE The partial differential equations your are solving
//! \param[in] bc_type The boundary condition type
//! \param[in] side The "side" a.k.a local face index
//! \param[in] qcoeff The coefficients for the state vector
//! \param[in] xcoeff The geometry coefficients
//! \param[in] resid The residual for all dofs of the element
//****************************************************************************80
   template < template<class, class>  class PDE_type >
   void CGBCResid(PDE_type<intT,realT>& PDE, const intT& bc_type, 
		   const intT& side, const Array2D<realT>& qcoeff, 
		   const Array2D<realT>& xcoeff, Array2D<realT>& resid) 
  {
    intT nf = PDE.get_nfld();
    Array1D<realT> q(nf), qb(nf);
    Array1D<realT> flux(nf);
    Array1D<realT> norm(ndim_);
    Array1D<realT> tang(ndim_);
    Array1D<realT> vg(ndim_);
    
    //---> Set grid velocity to zero
    resid.set_value(0.0); 
    vg.set_value(0.0);
    
    for (intT qp = 0; qp < nqp_face_ ; qp++){// quad_loop
           
      //---> Obtain function argument from coefficients and basis at quad_points
      ProjectToFaceQP<realT>(side, nf, qp, qcoeff, q);
      
      //---> Evaluate face mapping at the quadrature point
      ComputeFaceVectors(side, qp, xcoeff, norm, tang);
      
      //---> Obtain boundary state-vector
      PDE.template get_bc_q<realT>(bc_type, norm, tang, vg, q, qb);
      
      //---> Evauate \f$\vec{F} \cdot \vec{n}$\f
      flux.set_value(0.0);
      PDE.template FluxDotVec<realT>(norm, vg, qb, flux);
       
      for (intT f = 0; f < nf; f++) { // field_loop
	for(intT dof = 0; dof < ndof_; dof++) { // dof_loop 
	  /*---> Note norm is not a unit vector magnitude of face area term
	   present in the flux vector*/
	  resid(f,dof) += phi_face_(side, qp, dof)*flux(f)*wq_face_(qp);
	} // End dof_loop
      }// End field_loop
      
    } // End quad_loop
    return;

  }// End CGBCResid

//****************************************************************************80
//!
//! \brief CGVolJac : Computes the continuous Galerkin Volume residual and
//!                   Jacobian. 
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] PDE The partial differential equation class you are solving
//! \param[in] qcoeff The state vector coefficients 
//! \param[in] xcoeff The element geometry coefficients
//! \param[out] resid The integrated residual per degree of freedom
//! \param[out] fjac The integrated jacobian matrix 
//****************************************************************************80
  template < template<class, class>  class PDE_type >
  void CGVolJac(PDE_type<intT,realT>& PDE, 
		  const Array2D<realT>& qcoeff, 
		const Array2D<realT>& xcoeff, Array2D<realT>& resid,
		Array4D<realT>& fjac) {
    
    static const intT nf = PDE_type<intT,realT>::nfld;
    Array1D< Surreal<realT, nf> > q(nf);
    Array1D< Surreal<realT, nf> > flux(nf);
    Array1D<realT> norm(ndim_);
    Array1D<realT> vg(ndim_);
    Array1D< Surreal<realT, nf> > dq(nf);
    Array1D< Surreal<realT, nf> > divF(nf);
    Array2D< Surreal<realT, nf> > dqdxi(nf,ndim_);
    Array2D<realT> Jinv(ndim_,ndim_);
    SquareMatrix< Surreal<realT, nf> > tauinv(nf);
    realT DetJ;
       
    //---> Set integral Value to zero
    resid.set_value(0.0);
    fjac.set_value(0.0);
    vg.set_value(0.0);
          
    for (intT qp = 0; qp < nqp_ ; qp++){// quad_loop
      //---> Obtain function argument from coefficients and basis at quad_points
      ProjectToQP< Surreal<realT,nf> >(nf, qp, qcoeff, q);
   
      //---> Evaluate mapping at the quadrature point
      ComputeJinv(qp, xcoeff, Jinv, DetJ);
      
      //---> Compute dq/d(xi) at a the quadrature point 
      ProjectGradToQP< Surreal<realT,nf> >(nf, qp, qcoeff, dqdxi);

      //---> Set Divergence of F to zero
      divF.set_value(0.0);

      //---> Set derivative of q
      for(intT f = 0; f < nf; f++){// Deriv init
	q(f).set_deriv(f,1.0);
      }// End Deriv init
      //----------------------- Derivative w.r.t. q ----------------------------
      //---------------------- (-d(phi)/d(xi_{d})*F_{d}) -----------------------
      for (intT d = 0; d < ndim_; d++){// Dimension loop
	
	//---> Put the col of Jinv into normal vector
	for (intT d1 = 0; d1 < ndim_; d1++){ // Dimension loop-2
	  norm(d1) = Jinv(d1,d);
	}// End Dimension loop-2
	
	//---> Set flux to zero
	flux.set_value(0.0);
		
	//---> Get residual for PDE
	PDE.template FluxDotVec< Surreal<realT,nf> >(norm, vg, q, flux);

	//---> Assemble - d(dphi)/d(xi_{d})*F_{d}
	for (intT f = 0; f < nf; f++) {// Field loop 
	  for(intT dof = 0; dof < ndof_; dof++){// DoF loop 
	    resid(f,dof) -= (dphi_dxi_(qp,d,dof)*flux(f).Value())*DetJ*wq_(qp);
	  } // End Dof loop 
	} // End Field loop
	  
	//---> Jacobian;
	for (intT f = 0; f < nf; f++) {// Field loop 
	  for(intT fc = 0; fc < nf; fc++){// Col Field loop 
	    for(intT dof = 0; dof < ndof_; dof++){ // DoF loop
	      for(intT dofc = 0; dofc < ndof_; dofc++){// Col DoF loop 
		fjac(f,fc,dof,dofc) -= dphi_dxi_(qp,d,dof)*flux(f).Deriv(fc) * 
		  phi_(qp,dofc)*DetJ*wq_(qp);
	      }// End Col DoF loop
	    }// End Dof Loop
	  } // End Col Field loop
	} // End Field loop
	
	for (intT f = 0; f < nf; f++) {// Field loop  
	  dq(f) = dqdxi(f,d);
	}
	
	//---> Get the divergence term 
	PDE.template DerivFluxDotVecTimesX< Surreal<realT, nf> >
	  (norm, vg, q, dq, divF);
	
      }// End Dimension loop
          
      //------------------------ [tau]*(Rt(q,dq)) -----------------------------
      
      //---> Form [tau]^{-1}
      tauinv.set_value(0.0);
      
      for( intT dof = 0; dof < ndof_; dof++){ // DoF loop 
	//---> Form put grad(phi_dof) into norm by the following loop
	for( intT d = 0; d < ndim_; d++) {
	  norm(d) = 0.0;
	  for( intT d1 = 0; d1 < ndim_; d1++){
	    norm(d) += Jinv(d,d1)*dphi_dxi_(qp, d1, dof);
	  }
	}

	PDE.template EigenSystem< Surreal<realT, nf> >(norm, vg, q, tauinv);
	
      }// End DoF loop 
      
      flux.set_value(0.0);
      //---> [tau^{-1}]^{-1}*divF, store in flux;
      tauinv.squareSolve(divF, flux);
           
      //---------------- (d(phi)/d(xi_{d})*d(F_{d})/d(q)*tau*Rt(q,dq)) ---------
   
      for (intT d = 0; d < ndim_; d++){// Dimension loop-2
	//---> Put the col of Jinv into normal vector
	for (intT d1 = 0; d1 < ndim_; d1++){ // Dimension loop
	  norm(d1) = Jinv(d1,d);
	}// End Dimension loop

	divF.set_value(0.0);
	PDE.template DerivFluxDotVecTimesX< Surreal<realT, nf> >
	  (norm, vg, q, flux, divF);
	//---> Assemble  d(dphi)/d(xi_{d})*d(F_{d})(dq)*[tau]*R(q,dqdx);
	for (intT f = 0; f < nf; f++) {// Field loop 
	  for(intT dof = 0; dof < ndof_; dof++) { // Dof loop 
	    resid(f,dof) += (dphi_dxi_(qp,d,dof)*divF(f).Value())*DetJ*wq_(qp);
	  }// End Field loop 
	}// End DoF loop 
	
	//---> Jacobian;
	for (intT f = 0; f < nf; f++) {// Field loop 
	  for(intT fc = 0; fc < nf; fc++){// Col Field loop 
	    for(intT dof = 0; dof < ndof_; dof++){ // DoF loop
	      for(intT dofc = 0; dofc < ndof_; dofc++){// Col DoF loop 
		fjac(f,fc,dof,dofc) -= dphi_dxi_(qp,d,dof)*divF(f).Deriv(fc) * 
		  phi_(qp,dofc)*DetJ*wq_(qp);
	      }// End Col DoF loop
	    }// End Dof Loop
	  } // End Col Field loop
	} // End Field loop
	
      } // End Dimension loop-2 
      
      //------------------------ Derivative w.r.t. dqdxi -----------------------
      //---> Set derivative of q
      for(intT f = 0; f < nf; f++){// Deriv init
	q(f).set_deriv(f,0.0);
      }// End Deriv init

      //---> Loop over dimensions
      for(intT dim = 0; dim < ndim_; dim++){//Out dimension loop
	
	//---> Initalize dqdxi array to value only
	for(intT f = 0; f < nf; f++) { // Field init
	  for(intT d = 0; d < ndim_; d++){ // Dimension init
	    dqdxi(f,d) = dqdxi(f,d).Value();
	  } // End Dimension init
	} // End Field init
	
       	//---> Set derivative of dqdxi(:,dim) to identity
	for(intT f = 0; f < nf; f++){// Deriv init
	  dqdxi(f,dim).set_deriv(f,1.0);
	}// End Deriv init
	
	for (intT d = 0; d < ndim_; d++){// Dimension loop
	  //---> Put the col of Jinv into normal vector
	  for (intT d1 = 0; d1 < ndim_; d1++){ // Dimension loop-2
	    norm(d1) = Jinv(d1,d);
	  }// End Dimension loop-2
	  
	  //---> Set flux to zero
	  flux.set_value(0.0);
	  
	  //---> Get residual for PDE
	  //PDE.template FluxDotVec< Surreal<realT,nf> >(norm, vg, q, flux);
	  //---> Viscous backend to go here, flux(q,dqdxi);
	  
	  //---> Jacobian;
	  for (intT f = 0; f < nf; f++) {// Field loop 
	    for(intT fc = 0; fc < nf; fc++){// Col Field loop 
	      for(intT dof = 0; dof < ndof_; dof++){ // DoF loop
		for(intT dofc = 0; dofc < ndof_; dofc++){// Col DoF loop 
		  fjac(f,fc,dof,dofc) -= dphi_dxi_(qp,d,dof)*flux(f).Deriv(fc)* 
		    dphi_dxi_(qp,dim,dofc)*DetJ*wq_(qp);
		}// End Col DoF loop
	      }// End Dof Loop
	    } // End Col Field loop
	  } // End Field loop
	
	  for (intT f = 0; f < nf; f++) {// Field loop  
	    dq(f) = dqdxi(f,d);
	  }
	
	  //---> Get the divergence term 
	  PDE.template DerivFluxDotVecTimesX< Surreal<realT, nf> >
	    (norm, vg, q, dq, divF);
	
	}// End Dimension loop

	//----------------------- [tau]*(Rt(q,dq)) -----------------------------
	
	//---> Form [tau]^{-1}
	tauinv.set_value(0.0);
      
	for( intT dof = 0; dof < ndof_; dof++){ // DoF loop 
	  //---> Form put grad(phi_dof) into norm by the following loop
	  for( intT d = 0; d < ndim_; d++) {
	    norm(d) = 0.0;
	    for( intT d1 = 0; d1 < ndim_; d1++){
	      norm(d) += Jinv(d,d1)*dphi_dxi_(qp, d1, dof);
	    }
	  }

	  PDE.template EigenSystem< Surreal<realT, nf> >(norm, vg, q, tauinv);
	
	}// End DoF loop 
      
	flux.set_value(0.0);
	//---> [tau^{-1}]^{-1}*divF, store in flux;
	tauinv.squareSolve(divF, flux);
           
	//--------------- (d(phi)/d(xi_{d})*d(F_{d})/d(q)*tau*Rt(q,dq)) --------
   
	for (intT d = 0; d < ndim_; d++){// Dimension loop-2
	  //---> Put the col of Jinv into normal vector
	  for (intT d1 = 0; d1 < ndim_; d1++){ // Dimension loop
	    norm(d1) = Jinv(d1,d);
	  }// End Dimension loop

	  divF.set_value(0.0);
	  PDE.template DerivFluxDotVecTimesX< Surreal<realT, nf> >
	    (norm, vg, q, flux, divF);
	  //---> Jacobian;
	  for (intT f = 0; f < nf; f++) {// Field loop 
	    for(intT fc = 0; fc < nf; fc++){// Col Field loop 
	      for(intT dof = 0; dof < ndof_; dof++){ // DoF loop
		for(intT dofc = 0; dofc < ndof_; dofc++){// Col DoF loop 
		  fjac(f,fc,dof,dofc) -= dphi_dxi_(qp,d,dof)*divF(f).Deriv(fc)* 
		    dphi_dxi_(qp, dim, dofc)*DetJ*wq_(qp);
		}// End Col DoF loop
	      }// End Dof Loop
	    } // End Col Field loop
	  } // End Field loop
	
	} // End Dimension loop-2 
	
      }// End Outer dimension loop 
      
    }// End Quadrature loop
    return;
  } // End CGVolJac

//****************************************************************************80
//!
//! \brief CGBCResid : Computes the boundary residual for Continuous Galerkin
//!        methods.  
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] PDE The partial differential equations your are solving
//! \param[in] bc_type The boundary condition type
//! \param[in] side The "side" a.k.a local face index
//! \param[in] qcoeff The coefficients for the state vector
//! \param[in] xcoeff The geometry coefficients
//! \param[in] resid The residual for all dofs of the element
//****************************************************************************80
   template < template<class, class>  class PDE_type >
   void CGBCJac(PDE_type<intT,realT>& PDE, const intT& bc_type, 
		   const intT& side, const Array2D<realT>& qcoeff, 
		const Array2D<realT>& xcoeff, Array2D<realT>& resid, 
		Array4D<realT>& fjac) 
  {
    static const intT nf = PDE_type<intT,realT>::nfld;
    Array1D< Surreal<realT,nf> > q(nf), qb(nf);
    Array1D< Surreal<realT,nf> > flux(nf);
    Array1D<realT> norm(ndim_);
    Array1D<realT> tang(ndim_);
    Array1D<realT> vg(ndim_);
    
    //---> Set grid velocity to zero
    resid.set_value(0.0); 
    fjac.set_value(0.0);
    vg.set_value(0.0);
    
    for (intT qp = 0; qp < nqp_face_ ; qp++){// quad_loop
           
      //---> Obtain function argument from coefficients and basis at quad_points
      ProjectToFaceQP< Surreal<realT,nf> >(side, nf, qp, qcoeff, q);
      
      //---> Evaluate face mapping at the quadrature point
      ComputeFaceVectors(side, qp, xcoeff, norm, tang);
      
      //---> Set derivative of q
      for(intT f = 0; f < nf; f++){// Deriv init
	q(f).set_deriv(f,1.0);
      }// End Deriv init
      
      //---> Obtain boundary state-vector
      PDE.template get_bc_q< Surreal<realT,nf> >(bc_type, norm, tang, vg, q, 
						 qb);
      
      //---> Evauate \f$\vec{F} \cdot \vec{n}$\f
      flux.set_value(0.0);
      PDE.template FluxDotVec< Surreal<realT,nf> >(norm, vg, qb, flux);
       
      for (intT f = 0; f < nf; f++) { // field_loop
	for(intT dof = 0; dof < ndof_; dof++) { // dof_loop 
	  /*---> Note norm is not a unit vector magnitude of face area term
	   present in the flux vector*/
	  resid(f,dof) += phi_face_(side, qp, dof)*flux(f).Value()*wq_face_(qp);
	} // End dof_loop
      }// End field_loop
      
      //---> Jacobian;
      for (intT f = 0; f < nf; f++) {// Field loop 
	for(intT fc = 0; fc < nf; fc++){// Col Field loop 
	  for(intT dof = 0; dof < ndof_; dof++){ // DoF loop
	    for(intT dofc = 0; dofc < ndof_; dofc++){// Col DoF loop 
	      fjac(f, fc, dof, dofc) += phi_face_(side, qp, dof) *
		flux(f).Deriv(fc)*phi_face_(side, qp, dofc)*wq_face_(qp);
	    }// End Col DoF loop
	  }// End Dof Loop
	} // End Col Field loop
      } // End Field loop
      
    } // End quad_loop
    return;

  }// End CGBCJac

//****************************************************************************80
//!
//! \brief Volume : Computes the volume of an element given the mapping 
//!                 coefficients.
//! \details
//! \author Nick Burgess
//! \version $Rev$
//! \date $Date$
//! \param[in] xcoeff The mapping coefficients
//! \return vol The element volume
//****************************************************************************80
  realT Volume(const Array2D<realT>& xcoeff)
  {
    //---> Local Variables
    Array2D<realT> Jinv(ndim_, ndim_);
    realT DetJ;
    realT vol;
    
    //---> Set volume to zero
    vol = 0.0;
    for (intT j = 0; j < nqp_ ; j++){// quad_loop 
      //Evaluate mapping at the quadrature point
      ComputeJinv(j, xcoeff, Jinv, DetJ);
      
      vol += wq_(j)*DetJ;
    } // End quad_loop 
    
    return(vol);
  } // End Volume

//****************************************************************************80
//!
//! \brief FaceArea : Computes the volume of an element given the mapping 
//!                 coefficients.
//! \details
//! \author Nick Burgess
//! \version $Rev$
//! \date $Date$
//! \param[in] xcoeff The mapping coefficients
//! \return vol The element volume
//****************************************************************************80
  realT FaceArea(const intT& side, const Array2D<realT>& xcoeff)
  {
    //---> Local Variables
    Array1D<realT> norm(ndim_);
    Array1D<realT> tang(ndim_);
    realT Area;
    
    //---> Set volume to zero
    Area = 0.0;
    for (intT j = 0; j < nqp_face_ ; j++){// quad_loop 
      //Evaluate mapping at the quadrature point
      ComputeFaceVectors(side, j, xcoeff, norm, tang);
      realT DetJ = 0.0;
      for(intT d = 0; d < ndim_; d++) {DetJ += norm(d)*norm(d);}
      DetJ = sqrt(DetJ);
      
      Area += wq_face_(j)*DetJ;
    } // End quad_loop 
    
    return(Area);
  } // End integrate_func

//****************************************************************************80
//! \brief get_ndof : Returns the number of degrees of freedom for the element
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  inline intT get_ndof() {return ndof_;} // End get_ndof

}; // End class Element
#endif
