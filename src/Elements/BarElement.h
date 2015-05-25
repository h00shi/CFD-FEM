// -*-c++-*-

#ifndef BARELEMENT_H
#define BARELEMENT_H
//****************************************************************************80
//! \file BarElement.h 
//! \class BarElement BarElement.h
//! \brief This is the header file defining the class BarElement, 
//! which defines operators and data for a 1-D bar element.  
//! \details Ok, for some reason unknown to me C++ requires the "this->" on 
//! all member data inherited from base class Element.  I don't know why just 
//! keep  that in mind...if you miss one you'll get a compiler error telling 
//! you that that variable is undefined.  
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! 
//****************************************************************************80
#include "my_incl.h"
#include "consts.h"
#include "Array1D.h"
#include "Array2D.h"
#include "Array3D.h"
#include "Polynomials.h"
#include "Element.h"

template< typename intT, typename realT> 
class BarElement : public Element<intT, realT> {

private:
  //+++++++++++++++++++++++++++++++ PRIVATE STUFF ++++++++++++++++++++++++++++++
  Array1D<intT> mode2p_map; /*!< Map of mode number of p for 1-D Bar 
				  elements */

//****************************************************************************80
//!
//! \brief comp_ndof : Computes the number of dofs for a given order p 
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] p The degree of polynomial approximation
//****************************************************************************80
  intT comp_ndof(const intT& p)
  {
    return( (intT)(p + 1));
  }// End comp_ndof

//****************************************************************************80
//! \brief barH1 : Computes specified H1 1-D bar function at specified point 
//! \details This routine computes the 1-D H1 lobatto function 
//! \f$ l_{k}(\xi) \f$.  The general formula for the lobattor function is
//! \f$ l_{k}(\xi) = \frac{1-xi}{2}\frac{1+xi}{2}\psi_{k}(\xi) \f$ where
//! \f$ \psi_{k}(\xi) = P_{k -2}^{1, 1} \f$
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] k The function index you want i.e. you want \formula{l_{k}}
//! \param[in] xi The \f$ \xi \f$ coordinate \f$ \xi \in [-1,1] \f$ 
//! \return l The basis function at point xi \formula{l_{k}(\xi) } 
//****************************************************************************80

realT barH1(const intT& k, const realT& xi)
  {
    //---> Local Variables
    intT n;
    realT L0;
    realT L1;
    
    //---> Function result variable
    realT l;
    
    //---> Vertex functions
    L0 = (1.0 - xi)*half;
    L1 = (1.0 + xi)*half;
     
    //---> Based on the specified function  number requested
    switch(k) { // func_select
    case  0: // User is asking for mode 0
      l = L0;
      break;
    case 1: // Unser is asking for mode 1
      l = L1;
      break;
    default: // User is asking for mode > 1
      /*---> In 1D the order is k, since there a 2 vertex func + 1
	edge bubble */
      l = edgepoly(k, L0, L1);

      break;
    }// end func_select
    
    //---> Resturn result of fuction 
    return(l);
    
  } // End barH1

//****************************************************************************80
//!
//! \brief barH1D : Computes the derivative of specified H1 1-D bar function at
//!   specified point w.r.t. to \f$ \xi \f$ 
//! \details This routine computes the derivative of the 1-D H1 lobatto 
//!  function 
//! \f$ \frac{d l_{k}(\xi)}{d \xi} \f$.  The general formula for the lobattor 
//!  function is
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] k The function index you want i.e. you want \formula{l_{k}}
//! \param[in] xi The \f$ \xi \f$ coordinate \f$ \xi \in [-1,1] \f$ 
//! \return dl The basis function at point xi \f$ l_{k}(\xi) \f$ 
//****************************************************************************80
  realT barH1D(const intT& k, const realT& xi)
  {
    //---> Local Varialbes
    realT L0;
    realT L1;
    realT dL0dxi;
    realT dL1dxi;
    realT dldL0;
    realT dldL1;
    
    //---> Return variable
    realT dl;

    //---> Vertex functions
    L0 = (1.0 - xi)*half;
    L1 = (1.0 + xi)*half;
						
    dL0dxi = -half;
    dL1dxi = half;

    //---> Based on the specified function  number requested
    switch(k) { // func_select
    case  0: // User is asking for mode 0
      dl = dL0dxi;
      break;
    case 1: // Unser is asking for mode 1
      dl = dL1dxi;
      break;
    default: // User is asking for mode > 1
      
      edgepolyD(k, L0,  L1, dldL0, dldL1);
      
      dl = dldL0*dL0dxi + dldL1*dL1dxi;
      //dL0dxi*L1*psi + L0*dL1dxi*psi + half*L1*dpsi;
    }
    
    //---> Return the derivate to user
    return(dl);
  }// end barH1D

protected: 

  Array1D<intT> p2nqp; /*!<Map of number of quadrature points for a 
			     polynomial of degree deg */ 
//****************************************************************************80
//!
//! \brief bar_gauss_points : Tables of the bar gauss points, will populate 
//!                           The gauss points given a number of points
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] n The number of quadrature points
//! \param[out] xiq The coordinates of the quadrature points
//! \param[out] wq The quadrature weights
//****************************************************************************80
  void bar_gauss_points(const intT& n, Array2D<realT>& xiq, Array1D<realT>& wq)
  {
    //---> Use a switch case block to tablulate the gauss points
    switch (n) { // gauss_point_select
    case 1:
      //---> Points
      //Element<intT,realT>::xiq(0,0) = 0.0;
      xiq(0,0) = 0.0;
      //---> Weights
      //Element<intT,realT>::wq(0) = 2.0;
      wq(0) = 2.0;
      break;
    case 2: 
      //---> Points 
      // Element<intT,realT>::xiq(0,0) = -0.577350269189626;
      //Element<intT,realT>::xiq(1,0) =  0.577350269189626;
      xiq(0,0) = -0.577350269189626;
      xiq(1,0) =  0.577350269189626;
      //---> Weights
      //Element<intT,realT>::wq(0) = 1.0;
      //Element<intT,realT>::wq(1) = 1.0;
      wq(0) = 1.0;
      wq(1) = 1.0;
      break;
    case 3:
      //---> Points 
      //Element<intT,realT>::xiq(0,0) = -0.774596669241483;
      //Element<intT,realT>::xiq(1,0) =  0.000000000000000;
      //Element<intT,realT>::xiq(2,0) =  0.774596669241483;
      xiq(0,0) = -0.774596669241483;
      xiq(1,0) =  0.000000000000000;
      xiq(2,0) =  0.774596669241483;
      //---> Weights
      //Element<intT,realT>::wq(0) =  0.555555555555556;
      //Element<intT,realT>::wq(1) =  0.888888888888889;
      //Element<intT,realT>::wq(2) =  0.555555555555556;
      wq(0) =  0.555555555555556;
      wq(1) =  0.888888888888889;
      wq(2) =  0.555555555555556;
      break;
      
    }// End gauss_point_select
    
  } // End bar_gauss_points

//****************************************************************************80
//!
//! \brief edgepoly : Compute the value of a lobatto polynomial on an edge
//!                   given the affine coordinates L0, L1.  
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] p The degree of polynomial you want: MUST BE > 2.  
//! \param[in] L0 The first affine coordinate or "left" one
//! \param[in] L1 The second affine coorindat or "right" one
//! \return phi The edge polynomial value at evaluated at L1 - L0;
//****************************************************************************80
  realT edgepoly( const intT& p, const realT& L0, const realT& L1)
  {
    //---> Return to user the formula for an edge polynomial
    return( L0*L1*Polynomials::lobatto_kern<intT, realT>(p - 2, L1 - L0) );
  } // End edgepoly

//****************************************************************************80
//!
//! \brief edgepolyD : Compute the derivative of the edge polynomial w.r.t
//!                    affine coordintes L0, L1
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] p The degree of polynomila you : MUST BE > 2.
//! \param[in] L0 The first affine coordinate or "left" one
//! \param[in] L1 The second affine coorindat or "right" one
//! \param[out] dphidL0 \f$ \frac{ \partial phi}{\partial L0} \f$
//! \param[out] dphidL1 \f$ \frac{ \partial phi}{\partial L1} \f$
//****************************************************************************80
  void edgepolyD(const intT& p, const realT& L0, const realT& L1, 
		 realT& dphidL0, realT& dphidL1)
  {
    //---> Compute value of kernel function
    realT psi = Polynomials::lobatto_kern<intT, realT>(p - 2, L1 - L0);
    
    //---> Compute value of derivative of kernel function w.r.t (L1 - L0);
    realT dpsi = Polynomials::lobatto_kernD<intT, realT>(p - 2, L1 - L0);
    
    //---> Form partial derivative w.r.t. L0
    dphidL0 = L1*psi - L1*L0*dpsi;
    //---> Form partial derivative w.r.t. L1
    dphidL1 = L0*psi + L1*L0*dpsi;
    return;
  }


public:

//****************************************************************************80
//!
//! \brief BarElement : The constructor for the class BarElement
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! 
//****************************************************************************80
  BarElement() : Element<intT,realT>::Element()
  {
    //---> Initialize basis map data member
    mode2p_map.initialize(11);
    //---> p = 1
    mode2p_map(0) = 1;
    mode2p_map(1) = 1;
    
    //---> p = 2;
    mode2p_map(2) = 2;
    
    //---> p = 3;
    mode2p_map(3) = 3;
    
    //---> p = 4;
    mode2p_map(4) = 4;
    
    //---> p = 5;
    mode2p_map(5) = 5;
    
    //---> p = 6;
    mode2p_map(6) = 6;
    
    //---> p = 7;
    mode2p_map(7) = 7;

    //---> p = 8;
    mode2p_map(8) = 8;
    
    //---> p = 9;
    mode2p_map(9) = 9;
    
    //---> p = 10;
    mode2p_map(10) = 10;
    
    /* Initialize the array that keeps track of how many points required to 
       integrate a polynomail of degree ( ) */
    p2nqp.initialize(22);
    
    p2nqp(0) = 1;
    p2nqp(1) = 1;
    p2nqp(2) = 2;
    p2nqp(3) = 2;
    p2nqp(4) = 3;
    p2nqp(5) = 3;
    p2nqp(6) = 4;
    p2nqp(7) = 4;
    p2nqp(8) = 5;
    p2nqp(9) = 5;
    p2nqp(10) = 6;
    p2nqp(11) = 6;
    p2nqp(12) = 7;
    p2nqp(13) = 7;
    p2nqp(14) = 8;
    p2nqp(15) = 8;
    p2nqp(16) = 9;
    p2nqp(17) = 9;
    p2nqp(18) = 10;
    p2nqp(19) = 10;
    p2nqp(20) = 11;
    p2nqp(21) = 11;

  }; // End BarElement
  
//****************************************************************************80
//!
//! \brief ~BarElement : The destructor for the class BarElement
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! 
//****************************************************************************80
  ~BarElement()
  {
    
    
  } // End ~End

//****************************************************************************80
//!
//! \brief eval_basis : A public interface to evaluate a specified polynomial
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] n The polyomial you want to evaluate
//! \param[in] xi The standard element coordinate \f$\xi\f$ where you want the basis function value.
//! \return phi Value of basis functio n \f$ \phi_{n} \left( \xi\right) \f$
//****************************************************************************80
  realT eval_basis(const intT& n, const realT& xi)
  {
    //---> Evaluate the basis function 
    realT phi = barH1(n, xi);
    //---> Return the value of the polynomial specified by n at point xi
    return(phi);
  } // End eval_basis

//****************************************************************************80
//!
//! \brief eval_basisD : A public interface to evaluate a the derivate of 
//!                      a specified polynomail. 
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] n The polyomial you want to evaluate
//! \param[in] xi The standard element coordinate \f$\xi\f$ where you want the basis function derivate value.
//! \return phi Value of basis functio n \f$ \frac{ d \phi_{n}}{d \xi} \left( \xi\right) \f$
//****************************************************************************80
  realT eval_basisD(const intT& n, const realT& xi)
  {
    //---> Evaluate the basis function derivative
    realT dphi = barH1D(n, xi);
    //---> Return the value of the derivate of the basis function n at point xi
    return(dphi);
  }// End eval_basisD


//****************************************************************************80
//!
//! \brief initialize: Function to initialize an instance of the class
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] p_in The order of polynomials you want to integrate
//! \param[in] pmap_in The order of mapping polynomials
//! \param[in] deg_in The degree of integration required, recommend > 2*p_in.
//****************************************************************************80
  void initialize(const intT& p_in, const intT& pmap_in, const intT& deg_in) 
  {
    
    Element<intT,realT>::ndim = 1;

    //---> First set p to p_in 
    Element<intT,realT>::p = p_in;
    Element<intT,realT>::pmap = pmap_in;
    
    //---> Now set deg to deg_in
    Element<intT,realT>::deg = deg_in;

    //---> Compute ndof
    Element<intT,realT>::ndof = comp_ndof(Element<intT,realT>::p);
    Element<intT,realT>::ndof_map = comp_ndof(Element<intT,realT>::pmap);
    
    //---> Get number of quadrature points
    Element<intT,realT>::nqp = p2nqp(Element<intT,realT>::deg);
    Element<intT,realT>::nqp_face = 1;
    Element<intT,realT>::ndof_face = 1;
    
    //---> Now initialize quadrature vectors
    Element<intT,realT>::xiq.initialize(Element<intT,realT>::nqp, 
					Element<intT,realT>::ndim);

    Element<intT,realT>::wq.initialize(Element<intT,realT>::nqp);
    bar_gauss_points(Element<intT,realT>::nqp, Element<intT,realT>::xiq,
		     Element<intT,realT>::wq );
    
    Element<intT,realT>::xiq_face.initialize(1,1);
    Element<intT,realT>::wq_face.initialize(1);
    Element<intT,realT>::xiq_face(0,0) = 1.0;
    Element<intT,realT>::wq_face(0) = 1.0;
    
    //---> Now initialize the basis function vectors
    Element<intT,realT>::phi.initialize(Element<intT,realT>::nqp, 
					Element<intT,realT>::ndof);
    
    Element<intT,realT>::dphi_dxi.initialize(Element<intT,realT>::nqp, 
					     Element<intT,realT>::ndim, 
					     Element<intT,realT>::ndof);
    Element<intT,realT>::phi_face.initialize(2, 
					     Element<intT,realT>::nqp_face, 
					     Element<intT,realT>::ndof);
    Element<intT,realT>::dphi_dxi_face.initialize(2, 
						 Element<intT,realT>::nqp_face,
						 Element<intT,realT>::ndim,
						 Element<intT,realT>::ndof);
    //---> Now initialize the mapping basis functions
    Element<intT,realT>::phi_map.initialize(Element<intT,realT>::nqp, 
					    Element<intT,realT>::ndof_map);
    
    Element<intT,realT>::dphi_map.initialize(Element<intT,realT>::nqp, 
					     Element<intT,realT>::ndim, 
					     Element<intT,realT>::ndof_map);
    Element<intT,realT>::phi_map_face.initialize(2, 
						 Element<intT,realT>::nqp_face, 
						 Element<intT,realT>::ndof_map);
    Element<intT,realT>::dphi_map_face.initialize(2, 
						 Element<intT,realT>::nqp_face,
						 Element<intT,realT>::ndim,
						 Element<intT,realT>::ndof_map);
    Element<intT,realT>::dphi_ds.initialize(
					    Element<intT,realT>::nqp_face, 
					    1, 
					    Element<intT,realT>::ndof_face);
    //------------------------- Element Quadrature rules -----------------------

    for (intT i = 0; i < Element<intT,realT>::nqp; i++) { //qp_loop 
      for (intT j = 0; j < Element<intT,realT>::ndof; j++) { //dof_loop

	//---> Evalute basis dof: j at qp: i
	Element<intT,realT>::phi(i,j) =
	  eval_basis(j, Element<intT,realT>::xiq(i,0));
	//---> Evalute basisD dof: i at qp: j
	Element<intT,realT>::dphi_dxi(i,0,j) =
	  eval_basisD(j, Element<intT,realT>::xiq(i,0));
		
      }// End dof_loop 

      for (intT j = 0; j < Element<intT,realT>::ndof_map; j++){//dof_map_loop 
	//---> Evaluate mapping basis dof j at qp: i
	Element<intT,realT>::phi_map(i,j) = 
	  eval_basis(j, Element<intT,realT>::xiq(i,0));
	
	//---> Evalute mapping basisD dof: j at qp: i;
	Element<intT,realT>::dphi_map(i,0,j) = 
	  eval_basisD(j, Element<intT,realT>::xiq(i,0));
	
      } // End dof_map_loop 
    }// End qp_loop 
    
    //------------------------- Face Quadrature Rules --------------------------
    
    for (intT j = 0; j < Element<intT,realT>::ndof; j++){ //dof_loop
      //----> Evaluate basis at xi = -1;
      Element<intT,realT>::phi_face(0,0,j) = eval_basis(j,-1.0);
      Element<intT,realT>::dphi_dxi_face(0,0,0,j) = eval_basisD(j,-1.0);
    
      //----> Evaluate basis at xi = 1;
      Element<intT,realT>::phi_face(1,0,j) = eval_basis(j,1.0);
      Element<intT,realT>::dphi_dxi_face(1,0,0,j) = eval_basisD(j,1.0);
      
    }
    
    for (intT j = 0; j < Element<intT,realT>::ndof_map; j++){ //dof_loop
      //----> Evaluate basis at xi = -1;
      Element<intT,realT>::phi_map_face(0,0,j) = eval_basis(j,-1.0);
      Element<intT,realT>::dphi_map_face(0,0,0,j) = eval_basisD(j,-1.0);
      
      //----> Evaluate basis at xi = 1;
      Element<intT,realT>::phi_map_face(1,0,j) = eval_basis(j,1.0);
      Element<intT,realT>::dphi_map_face(1,0,0,j) = eval_basisD(j,1.0);
    }
  
    //------------------------ Face Parametrization Functions ------------------
    for (intT i = 0; i < Element<intT,realT>::nqp_face; i++) { //qp_loop
      for (intT j = 0; j < Element<intT,realT>::ndof_face; j++) { //dof_loop
	Element<intT,realT>::dphi_ds(i,0,j) = 1.0;
      }// End dof_loop
    }// End qp_loop
    
    //---> Face dof map.  For 1-D It's always 0 left side and 1 on right side.
    Element<intT,realT>::face_dof_map.initialize(2,1);
    Element<intT,realT>::face_dof_map(0,0) = 0;
    Element<intT,realT>::face_dof_map(1,0) = 1;
    return;
  }// End initialize
  
}; // End class BarElement

#endif
