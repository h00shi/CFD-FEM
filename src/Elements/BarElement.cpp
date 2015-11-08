//-*-c++-*-
#include "Elements/BarElement.h"
//****************************************************************************80
BarElement::BarElement() : Element::Element()
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
     integrate a polynomial of degree ( ) */
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
//****************************************************************************8
 BarElement::~BarElement()
 {
  
  } // End ~BarElement

//****************************************************************************80
intT BarElement::CompNdof(const intT& p)
{
  return( (intT)(p + 1));
}// End CompNdof

//****************************************************************************80
realT BarElement::BarH1(const intT& k, const realT& xi)
{
  //---> Local Variables
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
  case 1: // User is asking for mode 1
    l = L1;
    break;
  default: // User is asking for mode > 1
    /*---> In 1D the order is k, since there a 2 vertex func + 1
      edge bubble */
    l = EdgePoly(k, L0, L1);
    
    break;
  }// end func_select
  
  //---> Return result of function
  return(l);
  
} // End BarH1

//****************************************************************************80
realT BarElement::BarH1D(const intT& k, const realT& xi)
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
  case 1: // User is asking for mode 1
    dl = dL1dxi;
    break;
  default: // User is asking for mode > 1
    
    EdgePolyD(k, L0,  L1, dldL0, dldL1);
    
    dl = dldL0*dL0dxi + dldL1*dL1dxi;
    //dL0dxi*L1*psi + L0*dL1dxi*psi + half*L1*dpsi;
  }
  
  //---> Return the derivative to user
  return(dl);
}// end BarH1D

//****************************************************************************80
void BarElement::BarGaussPoints(const intT& n, Array2D<realT>& xiq,
				Array1D<realT>& wq)
{
  //---> Use a switch case block to tabulate the Gauss points
  switch (n) { // gauss_point_select
  case 1:
    //---> Points
    //Element::xiq(0,0) = 0.0;
    xiq(0,0) = 0.0;
    //---> Weights
    //Element::wq(0) = 2.0;
    wq(0) = 2.0;
    break;
  case 2: 
    //---> Points 
    // Element::xiq(0,0) = -0.577350269189626;
    //Element::xiq(1,0) =  0.577350269189626;
    xiq(0,0) = -0.577350269189626;
    xiq(1,0) =  0.577350269189626;
    //---> Weights
    //Element::wq(0) = 1.0;
    //Element::wq(1) = 1.0;
    wq(0) = 1.0;
    wq(1) = 1.0;
    break;
  case 3:
    //---> Points 
    //Element::xiq(0,0) = -0.774596669241483;
    //Element::xiq(1,0) =  0.000000000000000;
    //Element::xiq(2,0) =  0.774596669241483;
    xiq(0,0) = -0.774596669241483;
    xiq(1,0) =  0.000000000000000;
    xiq(2,0) =  0.774596669241483;
    //---> Weights
    //Element::wq(0) =  0.555555555555556;
    //Element::wq(1) =  0.888888888888889;
    //Element::wq(2) =  0.555555555555556;
    wq(0) =  0.555555555555556;
    wq(1) =  0.888888888888889;
    wq(2) =  0.555555555555556;
    break;
      
  }// End gauss_point_select
    
} // End BarGaussPoints
//****************************************************************************80
realT BarElement::EdgePoly( const intT& p, const realT& L0, const realT& L1)
{
  //---> Return to user the formula for an edge polynomial
  return( L0*L1*Polynomials::LobattoKern(p - 2, L1 - L0) );

} // End EdgePoly

//****************************************************************************80
void BarElement::EdgePolyD(const intT& p, const realT& L0, const realT& L1, 
			   realT& dphidL0, realT& dphidL1)
{
  //---> Compute value of kernel function
  realT psi = Polynomials::LobattoKern(p - 2, L1 - L0);
  
  //---> Compute value of derivative of kernel function w.r.t (L1 - L0);
  realT dpsi = Polynomials::LobattoKernD(p - 2, L1 - L0);
  
  //---> Form partial derivative w.r.t. L0
  dphidL0 = L1*psi - L1*L0*dpsi;
  //---> Form partial derivative w.r.t. L1
  dphidL1 = L0*psi + L1*L0*dpsi;
  return;
} // End EdgePolyD

//****************************************************************************80
realT BarElement::EvalBasis(const intT& n, const realT& xi)
  {
    //---> Evaluate the basis function 
    realT phi = BarH1(n, xi);
    //---> Return the value of the polynomial specified by n at point xi
    return(phi);
  } // End EvalBasis

//****************************************************************************80
realT BarElement::EvalBasisD(const intT& n, const realT& xi)
{
  //---> Evaluate the basis function derivative
  realT dphi = BarH1D(n, xi);
  //---> Return the value of the derivative of the basis function n at point xi
  return(dphi);
}// End EvalBasisD

//****************************************************************************80
void BarElement::initialize(const intT& p_in, const intT& pmap_in, 
			    const intT& deg_in) 
{
    
  Element::ndim_ = 1;
    
  //---> First set p to p_in 
  Element::p_ = p_in;
  Element::pmap_ = pmap_in;
    
  //---> Now set deg to deg_in
  Element::deg_ = deg_in;
    
  //---> Compute ndof
  Element::ndof_ = CompNdof(Element::p_);
  Element::ndof_map_ = CompNdof(Element::pmap_);
    
  //---> Get number of quadrature points
  Element::nqp_ = p2nqp(Element::deg_);
  Element::nqp_face_ = 1;
  Element::ndof_face_ = 1;
    
  //---> Now initialize quadrature vectors
  Element::xiq_.initialize(Element::nqp_, 
			  Element::ndim_);
    
  Element::wq_.initialize(Element::nqp_);
  BarGaussPoints(Element::nqp_, Element::xiq_,
		   Element::wq_ );
    
  Element::xiq_face_.initialize(1,1);
  Element::wq_face_.initialize(1);
  Element::xiq_face_(0,0) = 1.0;
  Element::wq_face_(0) = 1.0;
    
  //---> Now initialize the basis function vectors
  Element::phi_.initialize(Element::nqp_, 
			  Element::ndof_);
    
  Element::dphi_dxi_.initialize(Element::nqp_, 
			       Element::ndim_, 
			       Element::ndof_);
  Element::phi_face_.initialize(2, 
			       Element::nqp_face_, 
			       Element::ndof_);
  Element::dphi_dxi_face_.initialize(2, 
				    Element::nqp_face_,
				    Element::ndim_,
				    Element::ndof_);
  //---> Now initialize the mapping basis functions
  Element::phi_map_.initialize(Element::nqp_, 
			      Element::ndof_map_);
    
  Element::dphi_map_.initialize(Element::nqp_, 
			       Element::ndim_, 
			       Element::ndof_map_);
  Element::phi_map_face_.initialize(2, 
				   Element::nqp_face_, 
				   Element::ndof_map_);
  Element::dphi_map_face_.initialize(2, 
				    Element::nqp_face_,
				    Element::ndim_,
				    Element::ndof_map_);
  Element::dphi_ds_.initialize(
			      Element::nqp_face_, 
			      1, 
			      Element::ndof_face_);
  //------------------------- Element Quadrature rules -----------------------
    
  for (intT i = 0; i < Element::nqp_; i++) { //qp_loop 
    for (intT j = 0; j < Element::ndof_; j++) { //dof_loop
	
      //---> Evalute basis dof: j at qp: i
      Element::phi_(i,j) =
	EvalBasis(j, Element::xiq_(i,0));
      //---> Evalute basisD dof: i at qp: j
      Element::dphi_dxi_(i,0,j) =
	EvalBasisD(j, Element::xiq_(i,0));
	
    }// End dof_loop 
      
    for (intT j = 0; j < Element::ndof_map_; j++){//dof_map_loop 
      //---> Evaluate mapping basis dof j at qp: i
      Element::phi_map_(i,j) = 
	EvalBasis(j, Element::xiq_(i,0));
	
      //---> Evalute mapping basisD dof: j at qp: i;
      Element::dphi_map_(i,0,j) = 
	EvalBasisD(j, Element::xiq_(i,0));
	
    } // End dof_map_loop 
  }// End qp_loop 
    
  //------------------------- Face Quadrature Rules --------------------------
    
  for (intT j = 0; j < Element::ndof_; j++){ //dof_loop
    //----> Evaluate basis at xi = -1;
    Element::phi_face_(0,0,j) = EvalBasis(j,-1.0);
    Element::dphi_dxi_face_(0,0,0,j) = EvalBasisD(j,-1.0);
      
    //----> Evaluate basis at xi = 1;
    Element::phi_face_(1,0,j) = EvalBasis(j,1.0);
    Element::dphi_dxi_face_(1,0,0,j) = EvalBasisD(j,1.0);
      
  }
    
  for (intT j = 0; j < Element::ndof_map_; j++){ //dof_loop
    //----> Evaluate basis at xi = -1;
    Element::phi_map_face_(0,0,j) = EvalBasis(j,-1.0);
    Element::dphi_map_face_(0,0,0,j) = EvalBasisD(j,-1.0);
      
    //----> Evaluate basis at xi = 1;
    Element::phi_map_face_(1,0,j) = EvalBasis(j,1.0);
    Element::dphi_map_face_(1,0,0,j) = EvalBasisD(j,1.0);
  }
    
  //------------------------ Face Parametrization Functions ------------------
  for (intT i = 0; i < Element::nqp_face_; i++) { //qp_loop
    for (intT j = 0; j < Element::ndof_face_; j++) { //dof_loop
      Element::dphi_ds_(i,0,j) = 1.0;
    }// End dof_loop
  }// End qp_loop
    
  //---> Face dof map.  For 1-D It's always 0 left side and 1 on right side.
  Element::face_dof_map_.initialize(2,1);
  Element::face_dof_map_(0,0) = 0;
  Element::face_dof_map_(1,0) = 1;
  return;
}// End initialize
