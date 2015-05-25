// -*-c++-*-
#ifndef TRIELEMENT_H
#define TRIELEMENT_H
//****************************************************************************80
//! \file TriElement.h
//! \class TriElement TriElement.h
//! \brief This the header file defining the class TriElement, which defines operators and data for the a 2-D triangular element.
//! \details
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
#include "BarElement.h"
#include "Element.h"

template< typename intT, typename realT> 
class TriElement : public BarElement<intT, realT> {

 private:
 //+++++++++++++++++++++++++++++++ PRIVATE STUFF ++++++++++++++++++++++++++++++
  Array1D<intT>  mode2p_map; /*!< Map of mode number of p for 2-D Tri elements
			      */
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
    return( (p + 1)*(p + 2)/2);
  }// End comp_ndof
  
//****************************************************************************80
//!
//! \brief dof_type : A function to determine what type of heirarchical 
//!        function the specified dof k comes from.  
//! \details The functions that evaluate the basis only take in a dof number
//! called k and coordinates.  It is up to us to 
//! determine what type of heirarchical basis function dof k is.  For 
//! example for a p = 2 triangle mode 4 is an edge mode.  This function helps 
//! us determine that this is the case. 
//! For degree p we add a number of each type of dof as follows
//! \verbatim
//! There are 3 vertex dofs per triangle, these are first 3 0,1,2;
//! For p = 2 and higher we have (p - 1)*3 edge modes.  
//! For a degree p we add 3 edges modes to the previous basis set.  This is 
//! seen by nedge = (p - 1)*3, for p we have  ((p - 1) + 1) -1)*3 edge dofs,
//! which after some math is (p - 1 - 1)*3 + 3, i.e. 3 more edge dofs over 
//! what was in p - 1. 
//! \endverbatim 
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] k The dof number
//! \param[in] p The degree that mode k is from
//! \return mtype The mode type as an intT
//****************************************************************************80
  intT dof_type(const intT& k, const intT& p)
  {
    //---> Function return variable
    intT dtype;
    //---> First compute the number of dofs in space p - 1
    intT ndofpm1 = comp_ndof(p - 1);
    
    //---> Get local mode number
    intT li = k - ndofpm1;
    
    /*---> We know that we add 3 edge dofs per polynomial degree.  So if 
      li - 3 is < 0 then k is edge degree of freedom type 1. */  
    if( li - 3 < 0 ) {// determine_dtype
      //---> We found and edge
      dtype = 1;
    }
    else if( (li - 3) - (p - 2) < 0) {
      //---> We found a bubble mode
      dtype = 2;
    }
    else {
      //---> Oops something is wrong return a negative type to propogate error
      dtype = -1;
    } // end determin_type

    //---> Return dtype to user
    return(dtype);
  }// End dof_type

//****************************************************************************80
//!
//! \brief triH1 : Computes the value of a specified H1 basis function at
//!                a specified point \f$(\xi, \eta)\f$, in a triangle.  
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] k The index of the basis function
//! \param[in] xi The coordinate where you want the basis function evaluated
//! \param[in] eta The coorindate where you want the basis function evaluated
//! \return phi The basis function: \f$\phi(\xi,\eta)\f$
//****************************************************************************80
  realT triH1(const intT& k, const realT& xi, const realT& eta)
  {
    //---> Local Variables
    realT L[3]; 
    
    //---> Function return variable
    realT phi;
    
    //---> We'll need the vertex functions for any basis we have so specify them
    L[0] = -half*(xi + eta);
    L[1] = half*(xi + 1.0);
    L[2] = half*(eta + 1.0);
    
    //---> Get the discretization order based on specified mode
    intT p = mode2p_map(k);
    intT ndofpm1 = comp_ndof(p - 1);
 
    if( p == 1 ) { // order_check
      
      /*---> If p = 1 then the vertex functions are 
	encoded such that L[k] is the correct vertex 
	function. */
      phi = L[k];
    }
    else {

      intT dtype = dof_type(k,p);
      	
      //---> Compute edge index:
      intT eindex = k - ndofpm1;
     
      switch (dtype) { // Dof_type
      case 1: //----------------------- Edge dof, use edge -------------------- 
	switch(eindex) {// Edge_type
	case -1 :
	  //---> Error;
	  phi = -9.9e99;
	  break;
	case 0 : //---> Edge 0 of verticies 1 to 2 (0 to 1)
	  phi = BarElement<intT, realT>::edgepoly(p, L[0], L[1]);
	  break;
	case 1 : //---> Edge 1 of verticies 2 to 3 (1 to 2)
	  phi = BarElement<intT, realT>::edgepoly(p, L[1], L[2]);
	  break;
	case 2 : //---> Edge 2 or verticies 3 to 1 ( 2 to 0)
	  phi = BarElement<intT, realT>::edgepoly(p, L[2], L[0]);
	  break;
	} // End Edge_type
	
	break;
      case 2: //---> Bubble functions
	/*---> Find local bubble dof number (starting at 1) which is 
	  which is found by taking the current dof k and subtrating all the 
	  dofs up to now, which is -(ndofpm1 + 3), the 3 is for the 
	  3 edge dofs we've already added) */ 
	intT ibub = (k - ndofpm1) - 3 + 1;
	phi = tribubpoly(ibub, p, L[0], L[1], L[2]);

	break;
      } // End Dof_type

    }
    
    //---> Return value to user
    return(phi);
    
  } // End triH1;


//****************************************************************************80
//!
//! \brief triH1D : Derivative of H1 triangular basis function
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] k The index of the basis function
//! \param[in] xi The \f$ \xi \f$ coordinate \f$ \xi \in [-1,1] \f$ 
//! \param[in] eta The \f$ \eta \f$ coordinate  \f$ \eta \in [-1,1] \f$ 
//! \param[out] dphidxi The value: \f$ \frac{\partial \phi}{\partial \xi}\f$ at
//! \f$ (\xi, \eta) \f$
//! \param[out] dphideta The value: \f$ \frac{\partial \phi}{\partial \eta } \f$ at \f$ (\xi, \eta) \f$
//****************************************************************************80
  void triH1D(const intT& k, const realT& xi, const realT& eta, 
	      realT& dphidxi, realT& dphideta)
  {
    //---> Local Variables
    realT L[3];
    realT dLdxi[3];
    realT dLdeta[3];
    realT dphidL0; 
    realT dphidL1;
    realT dphidL2;
    
    //---> We'll need the vertex functions for any basis we have so specify them
    L[0] = -half*(xi + eta);
    L[1] = half*(xi + 1.0);
    L[2] = half*(eta + 1.0);
    
    dLdxi[0] = -half;
    dLdxi[1] = half;
    dLdxi[2] = 0.0;
    
    dLdeta[0] = -half;
    dLdeta[1] = 0.0;
    dLdeta[2] = half;
        
    //---> Get the discretization order based on specified mode
    intT p = mode2p_map(k);
    /*---> Compute number dofs in p - 1 discretization 
      (p - 1 + 1)*(p - 1 + 2)/2  */
    intT ndofpm1 = comp_ndof(p - 1);
    
    if( p == 1 ) { // order_check
      
      /*---> If p = 1 then the vertex functions are 
	encoded such that d(L[k])/(d \xi, \eta) is the correct vertex 
	function. */
      dphidxi = dLdxi[k];
      dphideta = dLdeta[k];
    }
    else {
      intT dtype = dof_type(k,p);
      	
      //---> Compute edge index:
      intT eindex = k - ndofpm1;
      switch (dtype) { // Dof_type
      case 1: //----------------------- Edge dof, use edge --------------------
	
	switch(eindex) { // Edge_type
	case -1 :
	  //---> Error;
	  dphidxi = -9.9e99;
	  dphideta = dphidxi;
	  
	  break;
	case 0 : //---> Edge 0 of verticies 0 to 1 (1 to 2)
	  BarElement<intT, realT>::edgepolyD(p, L[0], L[1], dphidL0, dphidL1);
	  //---> d(phi)/d(xi)
	  dphidxi = dphidL0*dLdxi[0] + dphidL1*dLdxi[1];
	  //---> d(phi)/d(eta)
	  dphideta = dphidL0*dLdeta[0] + dphidL1*dLdeta[1];
	  
	  break;
	case 1 : //---> Edge 1 of verticies 2 to 3 (1 to 2)
	  BarElement<intT, realT>::edgepolyD(p, L[1], L[2], dphidL0, dphidL1);
	  //---> d(phi)/d(xi)
	  dphidxi = dphidL0*dLdxi[1] + dphidL1*dLdxi[2];
	  //---> d(phi)/d(eta)
	  dphideta = dphidL0*dLdeta[1] + dphidL1*dLdeta[2];
       
	break;
      case 2 : //---> Edge 2 or verticies 3 to 1 ( 2 to 0)
	BarElement<intT, realT>::edgepolyD(p, L[2], L[0], dphidL0, dphidL1);
	//---> d(phi)/d(xi)
	dphidxi = dphidL0*dLdxi[2] + dphidL1*dLdxi[0];
	//---> d(phi)/d(eta)
	dphideta = dphidL0*dLdeta[2] + dphidL1*dLdeta[0];
	break;
	}// End Edge_type
	
	break;
      case 2: //---> Bubble functions
	/*---> Find local bubble dof number (starting at 1) which is 
	  which is found by taking the current dof k and subtrating all the 
	  dofs up to now, which is -(ndofpm1 + 3), the 3 is for the 
	  3 edge dofs we've already added) */ 
	intT ibub = (k - ndofpm1) - 3 + 1;
	tribubpolyD(ibub, p, L[0], L[1], L[2], dphidL0, dphidL1, dphidL2);
	
	//---> d(phi)/d(xi)
	dphidxi = dphidL0*dLdxi[0] + dphidL1*dLdxi[1] + dphidL2*dLdxi[2];
	
	//---> d(phi)/d(eta)
	dphideta =  dphidL0*dLdeta[0] + dphidL1*dLdeta[1] + dphidL2*dLdeta[2];
	
	break;
      
      } // End dof_type
      
    }// End order_check
    
    //---> Return nothing result goes out in dphi as a pass by reference. 
    return;
  } // End triH1D 
protected:
  
  Array1D<intT> p2nqp; /*!<Map of number of quadrature points for a 
			     polynomial of degree deg */ 
  
//****************************************************************************80
//!
//! \brief Tables of the triangle gauss points. Populates base class arrays
//!        given a number of points
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] n The number of points required
//! \param[out] xiq The quadrature point coordinates
//! \param[out] wq The quadrature weights
//****************************************************************************80
  void tri_gauss_points(const intT& n, Array2D<realT>& xiq, Array1D<realT>& wq)
  {
    //---> Use a switch case block to tablulate the gauss points
    switch (n) { // gauss_point_select
    case 1:
      //---> p=1
      //---> Points
      xiq(0,0) = -0.333333333333333; 
      xiq(0,1) = -0.333333333333333;
      
      //---> Weights
      wq(0) = 2.0;
      
      break;
    case 3: 
      //---> p=2
      //---> Points 
      xiq(0,0) = -0.666666666666667;
      xiq(0,1) = -0.666666666666667;
	
      xiq(1,0) = -0.666666666666667;
      xiq(1,1) =  0.333333333333333;
      
      xiq(2,0) =  0.333333333333333;
      xiq(2,1) = -0.666666666666667;
      
      //---> Weights
      wq(0) = 0.666666666666667;
      wq(1) = 0.666666666666667;
      wq(2) = 0.666666666666667;
      
      break;

    case 4:
      //---> p=3
      //---> Points 
      xiq(0,0) = -0.333333333333333;
      xiq(0,1) = -0.333333333333333;

      xiq(1,0) = -0.600000000000000;
      xiq(1,1) = -0.600000000000000;
      
      xiq(2,0) = -0.600000000000000;
      xiq(2,1) =  0.200000000000000;
      
      xiq(3,0) =  0.200000000000000;
      xiq(3,1) = -0.600000000000000;
      
      //---> Weights
      wq(0) = -1.125000000000000;
      wq(1) =  1.041666666666667;
      wq(2) =  1.041666666666667;
      wq(3) =  1.041666666666667;
     
      break;	

    case 6:
      //---> p=4
      //---> Points 
      xiq(0,0) = -0.108103018168070;
      xiq(0,1) = -0.108103018168070;
      
      xiq(1,0) = -0.108103018168070;
      xiq(1,1) = -0.783793963663860;
      
      xiq(2,0) = -0.783793963663860;
      xiq(2,1) = -0.108103018168070;
      
      xiq(3,0) = -0.816847572980458;
      xiq(3,1) =  0.816847572980458;
     
      xiq(4,0) = -0.816847572980458;
      xiq(4,1) =  0.633695145960918;
	
      xiq(5,0) =  0.633695145960918; 
      xiq(5,1) = -0.816847572980458;
	
      //---> Weights
      wq(0) = 0.446763179356022;
      wq(1) = 0.446763179356022;
      wq(2) = 0.446763179356022;
      wq(3) = 0.219903487310644;
      wq(4) = 0.219903487310644;
      wq(5) = 0.219903487310644;
      
      break;
      
    }// End gauss_point_select
    
  } // End bar_gauss_points

//****************************************************************************80
//!
//! \brief map_face_to_elem : Maps the face quadrature point coordinate to 
//!                           element coordinates
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] side The side of the element we are on 
//! \param[in] u The edge coordinate for the line segement
//! \param[out] xi The xi-direction coordinate in the triangle
//! \param[out] eta The eta-direction coordinate in the triangle
//****************************************************************************80
  void map_face_to_elem(const intT& side, const realT& u, realT& xi, realT& eta)
  {
    realT xi0;
    realT xi1;

    realT eta0;
    realT eta1;
    
    switch (side) { // Pick side
    case 0:
      xi0 = -1.0;
      xi1 = 1.0;
      
      eta0 = -1.0;
      eta1 = -1.0;
      break;
    case 1:
      xi0 = 1.0;
      xi1 = -1.0;
      
      eta0 = -1.0;
      eta1 = 1.0;
      break;
    case 2:
      xi0 = -1.0;
      xi1 = -1.0;
      
      eta0 = 1.0;
      eta1 = -1.0;
      break;
    }// Pick side

    xi = xi0*(1.0 - u)*half + xi1*(1.0 + u)*half;
    eta = eta0*(1.0 - u)*half + eta1*(1.0 + u)*half;

  } // End map_face_to_elem
 
//****************************************************************************80
//!
//! \brief tribubpoly : Computes the face functions for a triangle
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] ibub Specify which tri bubble function you want
//! \param[in] p The degree of bubble polynomial, there can be more than 1 per 
//!            degree
//! \param[in] L0 First affine coordinate
//! \param[in] L1 Second affine coordinate
//! \param[in] L2 Third affine coordinate
//! \return phi The tri bubble polynomial request
//****************************************************************************80
  realT tribubpoly(const intT& ibub, const intT& p, const realT& L0, 
		   const realT& L1, const realT& L2)
  {
    //---> Local Variables
    intT n1 = ibub;
    intT n2 = (p - 1) - n1;

    // //---> Now evaluate kernel function
	// realT psi1 = 
	//   Polynomials::lobatto_kern<intT, realT>(n1 - 1, L[1] - L[0]);
	// realT psi2 = 
	//   Polynomials::lobatto_kern<intT, realT>(n2 - 1, L[0] - L[2]);
	
	// phi = L[2]*L[0]*L[1]*psi1*psi2;
    return( L0*L1*L2*
	    Polynomials::lobatto_kern<intT, realT>(n1 - 1, L1 - L0)*
	    Polynomials::lobatto_kern<intT, realT>(n2 - 1, L0 - L2)
	    );
  }// End tribubpoly

//****************************************************************************80
//!
//! \brief tribubpolyD : Computes the derivative of tribubpoly w.r.t. all 3
//!        affine coordinates.  
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] ibub Specify which tri bubble function you want
//! \param[in] p The degree of bubble polynomial, there can be more than 1 per 
//!            degree
//! \param[in] L0 First affine coordinate
//! \param[in] L1 Second affine coordinate
//! \param[in] L2 Third affine coordinate
//! \param[out] dphidL0 \f$ \frac{ \partial phi}{\partial L0} \f$
//! \param[out] dphidL1 \f$ \frac{ \partial phi}{\partial L1} \f$
//! \param[out] dphidL2 \f$ \frac{ \partial phi}{\partial L2} \f$
//****************************************************************************80
  void tribubpolyD( const intT& ibub, const intT& p, const realT& L0, 
		    const realT& L1, const realT& L2, realT& dphidL0, 
		    realT& dphidL1, realT& dphidL2)
  {
    //---> Local Variables
    intT n1 = ibub;
    intT n2 = (p - 1) - n1;

    //---> Kernel functions
    realT psi1 = Polynomials::lobatto_kern<intT, realT>(n1 - 1, L1 - L0);
    realT psi2 = Polynomials::lobatto_kern<intT, realT>(n2 - 1, L0 - L2);

    //---> Derivative of Kernel functions w.r.t. their arguments
    realT dpsi1 = Polynomials::lobatto_kernD<intT, realT>(n1 - 1, L1 - L0);
    realT dpsi2 = Polynomials::lobatto_kernD<intT, realT>(n2 - 1, L0 - L2);
   
    dphidL0 = L1*L2*psi1*psi2 + L0*L1*L2*(-dpsi1*psi2 + psi1*dpsi2);
    dphidL1 = L0*L2*psi1*psi2 + L0*L1*L2*(dpsi1*psi2);
    dphidL2 = L0*L1*psi1*psi2 + L0*L1*L2*(-psi1*dpsi2);
    //---> Return nothing by value 
    return;
  }

public:

//****************************************************************************80
//!
//! \brief TriElement : The constructor for class tri-element
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! 
//****************************************************************************80
  TriElement() : BarElement<intT,realT>::BarElement()
  {
    //---> Initialize basis map data member
    mode2p_map.initialize(21);
    //---> p = 1;
    mode2p_map(0) = 1;
    mode2p_map(1) = 1;
    mode2p_map(2) = 1;
    
    //---> p = 2;
    mode2p_map(3) = 2;
    mode2p_map(4) = 2;
    mode2p_map(5) = 2;

    //---> p = 3;
    mode2p_map(6) = 3;
    mode2p_map(7) = 3;
    mode2p_map(8) = 3;
    mode2p_map(9) = 3;

    //---> p = 4;
    mode2p_map(10) = 4;
    mode2p_map(11) = 4;
    mode2p_map(12) = 4;
    mode2p_map(13) = 4;
    mode2p_map(14) = 4;
    
    //---> p = 5;
    mode2p_map(15) = 5;
    mode2p_map(16) = 5;
    mode2p_map(17) = 5;
    mode2p_map(18) = 5;
    mode2p_map(19) = 5;
    mode2p_map(20) = 5;
    
    //---> Initialize p2nqp
    p2nqp.initialize(6);
    p2nqp(0) = 1;
    p2nqp(1) = 1;
    p2nqp(2) = 3;
    p2nqp(3) = 4;
    p2nqp(4) = 6;
    p2nqp(5) = 7;
    
  }// End TriElement

//****************************************************************************80
//!
//! \brief ~TriElement : Destructor for class TriElement
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! 
//****************************************************************************80
  ~TriElement()
  {

  } // End ~TriElement


//****************************************************************************80
//! 
//! \brief A function to evaluate the selected basis function. 
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] n The basis function number you want
//! \param[in] xi The standard element coordinates \f$ \xi\f$
//! \param[in] eta The standard element coordinates \f$ \eta \f$
//! \return phi Value of basis function n \f$ \phi_{n}\left( (\xi,\eta)\right) 
//! \f$
//****************************************************************************80
  realT eval_basis(const intT& n, const realT&  xi, const realT& eta) 
  {
    //---> Evaluate basis function
    realT phi = triH1(n, xi, eta);
    
    //---> Return to user
    return(phi);
    
  } // End eval_basis

//****************************************************************************80
//!
//! \brief A function to evalute the selected basis function derivative
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! 
//****************************************************************************80
  void eval_basisD(const int& n, const realT& xi, const realT& eta, 
		   realT& dphidxi, realT& dphideta)
  {
    //---> Evaluate basis function derivatives
    triH1D(n, xi, eta, dphidxi, dphideta);
    
    return;
    
  } // End eval_basisD

//****************************************************************************80
//!
//! \brief Function to initialize an instance of the class
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
    Element<intT,realT>::ndim = 2;
    //---> First set p to p_in 
    Element<intT,realT>::p = p_in;
    Element<intT,realT>::pmap = pmap_in;
    
    //---> Now set deg to deg_in
    Element<intT,realT>::deg = deg_in;

    //---> Compute ndof
    Element<intT,realT>::ndof = comp_ndof(Element<intT,realT>::p);
    Element<intT,realT>::ndof_map = comp_ndof(Element<intT,realT>::pmap);
    Element<intT,realT>::ndof_face = Element<intT,realT>::pmap + 1;
 
    //---> Get number of quadrature points
    Element<intT,realT>::nqp = p2nqp(Element<intT,realT>::deg);
    Element<intT,realT>::nqp_face = 
      BarElement<intT,realT>::p2nqp(Element<intT,realT>::deg); 
							       
      
    //---> Now initialize quadrature vectors
    Element<intT,realT>::xiq.initialize(Element<intT,realT>::nqp, 
					Element<intT,realT>::ndim);

    Element<intT,realT>::wq.initialize(Element<intT,realT>::nqp);
    
    tri_gauss_points(Element<intT,realT>::nqp, 
		     Element<intT,realT>::xiq,
		     Element<intT,realT>::wq);
    
    Element<intT,realT>::xiq_face.initialize(Element<intT,realT>::nqp_face, 1);
    
    Element<intT,realT>::wq_face.initialize(Element<intT,realT>::nqp_face);
    
    BarElement<intT,realT>::bar_gauss_points(Element<intT,realT>::nqp_face, 
					     Element<intT,realT>::xiq_face, 
					     Element<intT,realT>::wq_face);
    
    //---> Now initialize the basis function vectors
    Element<intT,realT>::phi.initialize(Element<intT,realT>::nqp, 
					Element<intT,realT>::ndof);
    Element<intT,realT>::dphi_dxi.initialize(Element<intT,realT>::nqp, 
					     Element<intT,realT>::ndim, 
					     Element<intT,realT>::ndof);
    
    Element<intT,realT>::phi_face.initialize(3, 
					     Element<intT,realT>::nqp_face,
					     Element<intT,realT>::ndof);
    Element<intT,realT>::
      dphi_dxi_face.initialize(3, 
			       Element<intT,realT>::nqp_face, 
			       Element<intT,realT>::ndim,
			       Element<intT,realT>::ndof);
    
    //---> Now initialize the mapping basis functions
    Element<intT,realT>::phi_map.initialize(Element<intT,realT>::nqp, 
					    Element<intT,realT>::ndof_map);
    Element<intT,realT>::dphi_map.initialize(Element<intT,realT>::nqp, 
					     Element<intT,realT>::ndim, 
					     Element<intT,realT>::ndof_map);
    
    Element<intT,realT>::phi_map_face.initialize(3, 
						 Element<intT,realT>::nqp_face, 
						 Element<intT,realT>::ndof_map);
    Element<intT,realT>::
      dphi_map_face.initialize(3, 
			       Element<intT,realT>::nqp_face, 
			       Element<intT,realT>::ndim, 
			       Element<intT,realT>::ndof_map);
    
    Element<intT,realT>::dphi_ds.initialize(Element<intT,realT>::nqp_face, 
    				Element<intT,realT>::ndim - 1, 
    				Element<intT,realT>::ndof_face);
    
    //------------------------- Element Quadrature rules -----------------------
   
    for (intT i = 0; i < Element<intT,realT>::nqp; i++) { //qp_loop 
      for (intT j = 0; j < Element<intT,realT>::ndof; j++) { //dof_loop

	//---> Evalute basis dof: j at qp:i
	Element<intT,realT>::phi(i,j) = 
	  eval_basis(j, Element<intT,realT>::xiq(i,0), 
		     Element<intT,realT>::xiq(i,1) );

	//---> Evalute basisD dof: j at qp: i
	eval_basisD(j, Element<intT,realT>::xiq(i,0), 
		    Element<intT,realT>::xiq(i,1), 
		    Element<intT,realT>::dphi_dxi(i,0,j), 
		    Element<intT,realT>::dphi_dxi(i,1,j) );
	
      }// End dof_loop 

      for (intT j = 0; j < Element<intT,realT>::ndof_map; j++){//dof_map_loop 
	//---> Evaluate mapping basis dof j at qp: i
	Element<intT,realT>::phi_map(i,j) = 
	  eval_basis(j, Element<intT,realT>::xiq(i,0), 
		     Element<intT,realT>::xiq(i,1) );
	
	//---> Evalute mapping basisD dof: i at qp: j;
	eval_basisD(j, Element<intT,realT>::xiq(j,0), 
		    Element<intT,realT>::xiq(j,1), 
		    Element<intT,realT>::dphi_map(i,0,j), 
		    Element<intT,realT>::dphi_map(i,1,j));
      } // End dof_map_loop 
      
    }// End qp_loop 

    //------------------------- Face Quadrature Rules --------------------------
   
    for(intT side = 0; side < 3; side++){// side_loop
      for (intT i = 0; i < Element<intT,realT>::nqp_face; i++) { //qp_loop
	realT xi, eta;
	map_face_to_elem(side, Element<intT,realT>::xiq_face(i,0), xi, eta);
	
	for (intT j = 0; j < Element<intT,realT>::ndof; j++) { //dof_loop
	  
	  //---> Evalute basis dof: j at qp:i
	   Element<intT,realT>::phi_face(side, i, j) = eval_basis(j, xi, eta);
	   eval_basisD(j, xi, eta, 
		       Element<intT,realT>::dphi_dxi_face(side, i, 0, j),
		       Element<intT,realT>::dphi_dxi_face(side, i, 1, j) );
	}  //End dof_loop
	
	for (intT j = 0; j < Element<intT,realT>::ndof_map; j++){//dof_map_loop 
	  //---> Evaluate mapping basis dof j at qp: i
	  Element<intT,realT>::phi_map_face(side, i, j) = eval_basis(j,xi,eta);
	  eval_basisD(j, xi, eta, 
		      Element<intT,realT>::dphi_map_face(side, i, 0, j),
		      Element<intT,realT>::dphi_map_face(side, i, 1, j) );
	  
	}// End dof_map_loop
      } // End qp_loop
    }// End side_loop
    
    //------------------------ Face Parametrization Functions ------------------
    for (intT i = 0; i < Element<intT,realT>::nqp_face; i++) { //qp_loop
      for (intT j = 0; j < Element<intT,realT>::ndof_face; j++) { //dof_loop
	Element<intT,realT>::dphi_ds(i,0,j) = 
	  BarElement<intT,realT>::
	  eval_basisD(j, Element<intT,realT>::xiq_face(i,0));
      }// End dof_loop
    }// End qp_loop
    //---> Face dof map. In 2-D for triangles it is as follows
    Element<intT,realT>::face_dof_map.initialize(3,2); //---> P = 1 only right now
    //---> Side 0
    Element<intT,realT>::face_dof_map(0,0) = 0;
    Element<intT,realT>::face_dof_map(0,1) = 1;
    
    //---> Side 1
    Element<intT,realT>::face_dof_map(1,0) = 1;
    Element<intT,realT>::face_dof_map(1,1) = 2;
    
    //---> Side 2
    Element<intT,realT>::face_dof_map(2,0) = 2;
    Element<intT,realT>::face_dof_map(2,1) = 0;
    return;
  }// End initialize
  
}; // End class TriElement

#endif
