#include "TriElement.h"
//****************************************************************************80
TriElement::TriElement() : BarElement::BarElement()
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
TriElement::~TriElement()
{
} // End ~TriElement

//****************************************************************************80
intT TriElement::CompNdof(const intT& p)
{
  return( (p + 1)*(p + 2)/2);
}// End CompNdof

//****************************************************************************80
intT TriElement::DofType(const intT& k, const intT& p)
{
    //---> Function return variable
  intT dtype;
  //---> First compute the number of dofs in space p - 1
  intT ndofpm1 = CompNdof(p - 1);
  
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
}// End DofType

//****************************************************************************80
realT TriElement::TriH1(const intT& k, const realT& xi, const realT& eta)
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
  intT ndofpm1 = CompNdof(p - 1);
 
  if( p == 1 ) { // order_check
      
    /*---> If p = 1 then the vertex functions are 
      encoded such that L[k] is the correct vertex 
      function. */
    phi = L[k];
  }
  else {

    intT dtype = DofType(k,p);
      	
    //---> Compute edge index:
    intT eindex = k - ndofpm1;
     
    switch (dtype) { // DofType
    case 1: //----------------------- Edge dof, use edge -------------------- 
      switch(eindex) {// Edge_type
      case -1 :
	//---> Error;
	phi = -9.9e99;
	break;
      case 0 : //---> Edge 0 of verticies 1 to 2 (0 to 1)
	phi = BarElement::EdgePoly(p, L[0], L[1]);
	break;
      case 1 : //---> Edge 1 of verticies 2 to 3 (1 to 2)
	phi = BarElement::EdgePoly(p, L[1], L[2]);
	break;
      case 2 : //---> Edge 2 or verticies 3 to 1 ( 2 to 0)
	phi = BarElement::EdgePoly(p, L[2], L[0]);
	break;
      } // End Edge_type
	
      break;
    case 2: //---> Bubble functions
      /*---> Find local bubble dof number (starting at 1) which is 
	which is found by taking the current dof k and subtrating all the 
	dofs up to now, which is -(ndofpm1 + 3), the 3 is for the 
	3 edge dofs we've already added) */ 
      intT ibub = (k - ndofpm1) - 3 + 1;
      phi = TriBubPoly(ibub, p, L[0], L[1], L[2]);

      break;
    } // End DofType

  }
    
  //---> Return value to user
  return(phi);
    
} // End TriH1;

//****************************************************************************80
void TriElement::TriH1D(const intT& k, const realT& xi, const realT& eta, 
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
  intT ndofpm1 = CompNdof(p - 1);
    
  if( p == 1 ) { // order_check
      
    /*---> If p = 1 then the vertex functions are 
      encoded such that d(L[k])/(d \xi, \eta) is the correct vertex 
      function. */
    dphidxi = dLdxi[k];
    dphideta = dLdeta[k];
  }
  else {
    intT dtype = DofType(k,p);
      	
    //---> Compute edge index:
    intT eindex = k - ndofpm1;
    switch (dtype) { // DofType
    case 1: //----------------------- Edge dof, use edge --------------------
	
      switch(eindex) { // Edge_type
      case -1 :
	//---> Error;
	dphidxi = -9.9e99;
	dphideta = dphidxi;
	  
	break;
      case 0 : //---> Edge 0 of verticies 0 to 1 (1 to 2)
	BarElement::EdgePolyD(p, L[0], L[1], dphidL0, dphidL1);
	//---> d(phi)/d(xi)
	dphidxi = dphidL0*dLdxi[0] + dphidL1*dLdxi[1];
	//---> d(phi)/d(eta)
	dphideta = dphidL0*dLdeta[0] + dphidL1*dLdeta[1];
	  
	break;
      case 1 : //---> Edge 1 of verticies 2 to 3 (1 to 2)
	BarElement::EdgePolyD(p, L[1], L[2], dphidL0, dphidL1);
	//---> d(phi)/d(xi)
	dphidxi = dphidL0*dLdxi[1] + dphidL1*dLdxi[2];
	//---> d(phi)/d(eta)
	dphideta = dphidL0*dLdeta[1] + dphidL1*dLdeta[2];
       
	break;
      case 2 : //---> Edge 2 or verticies 3 to 1 ( 2 to 0)
	BarElement::EdgePolyD(p, L[2], L[0], dphidL0, dphidL1);
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
      TriBubPolyD(ibub, p, L[0], L[1], L[2], dphidL0, dphidL1, dphidL2);
	
      //---> d(phi)/d(xi)
      dphidxi = dphidL0*dLdxi[0] + dphidL1*dLdxi[1] + dphidL2*dLdxi[2];
	
      //---> d(phi)/d(eta)
      dphideta =  dphidL0*dLdeta[0] + dphidL1*dLdeta[1] + dphidL2*dLdeta[2];
	
      break;
      
    } // End DofType
      
  }// End order_check
    
  //---> Return nothing result goes out in dphi as a pass by reference. 
  return;
} // End TriH1D 

//****************************************************************************80
void TriElement::TriGaussPoints(const intT& n, Array2D<realT>& xiq, 
				Array1D<realT>& wq)
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
    
} // End TriGaussPoints

//****************************************************************************80
void TriElement::MapFaceToElem(const intT& side, const realT& u, realT& xi, 
			       realT& eta)
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
  
} // End MapFaceToElem
//****************************************************************************80
realT TriElement::TriBubPoly(const intT& ibub, const intT& p, const realT& L0, 
			     const realT& L1, const realT& L2)
{
  //---> Local Variables
  intT n1 = ibub;
  intT n2 = (p - 1) - n1;
  
  // //---> Now evaluate kernel function
  // realT psi1 = 
  //   Polynomials::LobattoKern(n1 - 1, L[1] - L[0]);
  // realT psi2 = 
  //   Polynomials::LobattoKern(n2 - 1, L[0] - L[2]);
  
  // phi = L[2]*L[0]*L[1]*psi1*psi2;
  return( L0*L1*L2*
	  Polynomials::LobattoKern(n1 - 1, L1 - L0)*
	  Polynomials::LobattoKern(n2 - 1, L0 - L2)
	  );
}// End TriBubPoly

//****************************************************************************80
void TriElement::TriBubPolyD( const intT& ibub, const intT& p, const realT& L0, 
			      const realT& L1, const realT& L2, realT& dphidL0, 
			      realT& dphidL1, realT& dphidL2)
{
  //---> Local Variables
  intT n1 = ibub;
  intT n2 = (p - 1) - n1;
  
  //---> Kernel functions
  realT psi1 = Polynomials::LobattoKern(n1 - 1, L1 - L0);
  realT psi2 = Polynomials::LobattoKern(n2 - 1, L0 - L2);
  
  //---> Derivative of Kernel functions w.r.t. their arguments
  realT dpsi1 = Polynomials::LobattoKernD(n1 - 1, L1 - L0);
  realT dpsi2 = Polynomials::LobattoKernD(n2 - 1, L0 - L2);
  
  dphidL0 = L1*L2*psi1*psi2 + L0*L1*L2*(-dpsi1*psi2 + psi1*dpsi2);
  dphidL1 = L0*L2*psi1*psi2 + L0*L1*L2*(dpsi1*psi2);
  dphidL2 = L0*L1*psi1*psi2 + L0*L1*L2*(-psi1*dpsi2);
  //---> Return nothing by value 
  return;
}

//****************************************************************************80
realT TriElement::EvalBasis(const intT& n, const realT&  xi, const realT& eta) 
{
  //---> Evaluate basis function
  realT phi = TriH1(n, xi, eta);
  
  //---> Return to user
  return(phi);
  
} // End EvalBasis

//****************************************************************************80
void TriElement::EvalBasisD(const int& n, const realT& xi, const realT& eta, 
			    realT& dphidxi, realT& dphideta)
{
  //---> Evaluate basis function derivatives
  TriH1D(n, xi, eta, dphidxi, dphideta);
  
  return;
  
} // End EvalBasisD

//****************************************************************************80
void TriElement::initialize(const intT& p_in, const intT& pmap_in, 
			    const intT& deg_in) 
{
  Element::ndim_ = 2;
  //---> First set p to p_in 
  Element::p_ = p_in;
  Element::pmap_ = pmap_in;
    
  //---> Now set deg to deg_in
  Element::deg_ = deg_in;

  //---> Compute ndof
  Element::ndof_ = CompNdof(Element::p_);
  Element::ndof_map_ = CompNdof(Element::pmap_);
  Element::ndof_face_ = Element::pmap_ + 1;
 
  //---> Get number of quadrature points
  Element::nqp_ = p2nqp(Element::deg_);
  Element::nqp_face_ = BarElement::p2nqp(Element::deg_); 
							       
      
  //---> Now initialize quadrature vectors
  Element::xiq_.initialize(Element::nqp_, 
			  Element::ndim_);

  Element::wq_.initialize(Element::nqp_);
    
  TriGaussPoints(Element::nqp_, 
		   Element::xiq_,
		   Element::wq_);
    
  Element::xiq_face_.initialize(Element::nqp_face_, 1);
    
  Element::wq_face_.initialize(Element::nqp_face_);
    
  BarElement::BarGaussPoints(Element::nqp_face_, 
			     Element::xiq_face_, 
			     Element::wq_face_);
    
  //---> Now initialize the basis function vectors
  Element::phi_.initialize(Element::nqp_, 
			   Element::ndof_);
  Element::dphi_dxi_.initialize(Element::nqp_, 
				Element::ndim_, 
				Element::ndof_);
    
  Element::phi_face_.initialize(3, 
			       Element::nqp_face_,
			       Element::ndof_);
  Element::dphi_dxi_face_.initialize(3, 
				     Element::nqp_face_, 
				     Element::ndim_,
				     Element::ndof_);
    
  //---> Now initialize the mapping basis functions
  Element::phi_map_.initialize(Element::nqp_, 
			       Element::ndof_map_);
  Element::dphi_map_.initialize(Element::nqp_, 
				Element::ndim_, 
				Element::ndof_map_);
    
  Element::phi_map_face_.initialize(3, Element::nqp_face_, Element::ndof_map_);
  Element::dphi_map_face_.initialize(3, Element::nqp_face_, Element::ndim_, 
				     Element::ndof_map_);
    
  Element::dphi_ds_.initialize(Element::nqp_face_, 
			       Element::ndim_ - 1, 
			       Element::ndof_face_);
    
  //------------------------- Element Quadrature rules -----------------------
   
  for (intT i = 0; i < Element::nqp_; i++) { //qp_loop 
    for (intT j = 0; j < Element::ndof_; j++) { //dof_loop

      //---> Evalute basis dof: j at qp:i
      Element::phi_(i,j) = EvalBasis(j, Element::xiq_(i,0), Element::xiq_(i,1));

      //---> Evalute basisD dof: j at qp: i
      EvalBasisD(j, Element::xiq_(i,0), Element::xiq_(i,1), 
		 Element::dphi_dxi_(i,0,j), 
		 Element::dphi_dxi_(i,1,j) );
      
    }// End dof_loop 

    for (intT j = 0; j < Element::ndof_map_; j++){//dof_map_loop 
      //---> Evaluate mapping basis dof j at qp: i
      Element::phi_map_(i,j) = EvalBasis(j, Element::xiq_(i,0), 
					 Element::xiq_(i,1) );
	
      //---> Evalute mapping basisD dof: i at qp: j;
      EvalBasisD(j, Element::xiq_(i,0), 
		  Element::xiq_(i,1), 
		  Element::dphi_map_(i,0,j), 
		  Element::dphi_map_(i,1,j));
    } // End dof_map_loop 
      
  }// End qp_loop 

  //------------------------- Face Quadrature Rules --------------------------
   
  for(intT side = 0; side < 3; side++){// side_loop
    for (intT i = 0; i < Element::nqp_face_; i++) { //qp_loop
      realT xi, eta;
      MapFaceToElem(side, Element::xiq_face_(i,0), xi, eta);
	
      for (intT j = 0; j < Element::ndof_; j++) { //dof_loop
	  
	//---> Evalute basis dof: j at qp:i
	Element::phi_face_(side, i, j) = EvalBasis(j, xi, eta);
	EvalBasisD(j, xi, eta, Element::dphi_dxi_face_(side, i, 0, j),
		   Element::dphi_dxi_face_(side, i, 1, j) );
      }  //End dof_loop
	
      for (intT j = 0; j < Element::ndof_map_; j++){//dof_map_loop 
	//---> Evaluate mapping basis dof j at qp: i
	Element::phi_map_face_(side, i, j) = EvalBasis(j,xi,eta);
	EvalBasisD(j, xi, eta, 
		   Element::dphi_map_face_(side, i, 0, j),
		   Element::dphi_map_face_(side, i, 1, j) );
	  
      }// End dof_map_loop
    } // End qp_loop
  }// End side_loop
    
  //------------------------ Face Parametrization Functions ------------------
  for (intT i = 0; i < Element::nqp_face_; i++) { //qp_loop
    for (intT j = 0; j < Element::ndof_face_; j++) { //dof_loop
      Element::dphi_ds_(i,0,j) = BarElement::
	EvalBasisD(j, Element::xiq_face_(i,0));
    }// End dof_loop
  }// End qp_loop
  //---> Face dof map. In 2-D for triangles it is as follows
  Element::face_dof_map_.initialize(3,2); /* P = 1 only 
					    right now */
  //---> Side 0
  Element::face_dof_map_(0,0) = 0;
  Element::face_dof_map_(0,1) = 1;
    
  //---> Side 1
  Element::face_dof_map_(1,0) = 1;
  Element::face_dof_map_(1,1) = 2;
    
  //---> Side 2
  Element::face_dof_map_(2,0) = 2;
  Element::face_dof_map_(2,1) = 0;
  return;
}// End initialize
