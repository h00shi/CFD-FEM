#include "TetElement.h"
//****************************************************************************80
TetElement::TetElement() : TriElement::TriElement()
{
  //---> Initialize basis map data member
  mode2p_map.initialize(56);
    
  //---> p = 1;
  mode2p_map(0) = 1;
  mode2p_map(1) = 1;
  mode2p_map(2) = 1;
  mode2p_map(3) = 1;

  //---> p = 2;
  mode2p_map(4) = 2;
  mode2p_map(5) = 2;
  mode2p_map(6) = 2;
  mode2p_map(7) = 2;
  mode2p_map(8) = 2;
  mode2p_map(9) = 2;

  //---> p = 3;
  mode2p_map(10) = 3;
  mode2p_map(11) = 3;
  mode2p_map(12) = 3;
  mode2p_map(13) = 3;
  mode2p_map(14) = 3;
  mode2p_map(15) = 3;
  mode2p_map(16) = 3;
  mode2p_map(17) = 3;
  mode2p_map(18) = 3;
  mode2p_map(19) = 3;
    
  //---> p = 4;
  mode2p_map(20) = 4;
  mode2p_map(21) = 4;
  mode2p_map(22) = 4;
  mode2p_map(23) = 4;
  mode2p_map(24) = 4;
  mode2p_map(25) = 4;
  mode2p_map(26) = 4;
  mode2p_map(27) = 4;
  mode2p_map(28) = 4;
  mode2p_map(29) = 4;
  mode2p_map(30) = 4;
  mode2p_map(31) = 4;
  mode2p_map(32) = 4;
  mode2p_map(33) = 4;
  mode2p_map(34) = 4;

  //---> p = 5;
  mode2p_map(35) = 5;
  mode2p_map(36) = 5;
  mode2p_map(37) = 5;
  mode2p_map(38) = 5;
  mode2p_map(39) = 5;
  mode2p_map(40) = 5;
  mode2p_map(41) = 5;
  mode2p_map(42) = 5;
  mode2p_map(43) = 5;
  mode2p_map(44) = 5;
  mode2p_map(45) = 5;
  mode2p_map(46) = 5;
  mode2p_map(47) = 5;
  mode2p_map(48) = 5;
  mode2p_map(49) = 5;
  mode2p_map(50) = 5;
  mode2p_map(51) = 5;
  mode2p_map(52) = 5;
  mode2p_map(53) = 5;
  mode2p_map(54) = 5;
  mode2p_map(55) = 5; 

  //---> Initialize p2nqp
  p2nqp.initialize(6);
  p2nqp(0) = 1;
  p2nqp(1) = 1;
  p2nqp(2) = 4;
  p2nqp(3) = 5;
  p2nqp(4) = 11;
  p2nqp(5) = 14;
       
}//End TetElement

//****************************************************************************80
TetElement::~TetElement()
{
  
} // End ~TetElement

//****************************************************************************80
intT TetElement::CompNdof(const intT& p)
{
  return( (p + 1)*(p + 2)*(p + 3)/6);
}// End CompNdof
//****************************************************************************80
intT TetElement::DofType(const intT& k, const intT& p)
{
  //---> Function return variable
  intT dtype;
  //---> First compute the number of dofs in space p - 1
  intT ndofpm1 = CompNdof(p - 1);
  
  //---> Get local mode number
  intT li = k - ndofpm1;
  
  intT ne = 6; //Number of edge dofs per degree p
  intT nf = (p - 2)*4; // Number of face dofs per degree p
  intT nb = (p - 2)*(p - 3)/2; // Number of bubble dofs per degree p
  
  /*---> We know that we add 3 edge dofs per polynomial degree.  So if 
    li - 3 is < 0 then k is edge degree of freedom type 1. */  
  if( li - ne < 0 ) {// determine_dtype
    //---> We found and edge dof
    dtype = 1;
  }
  else if( li - ne - nf < 0) {
    //---> We found a face dof
    dtype = 2;
  }
  else if( li - ne - nf - nb < 0 ) {
    //---> We found a bubble mode
    dtype = 3; 
  }
  else {
    //---> Oops something is wrong return a negative type to propogate error
    dtype = -1;
  } // end determin_type
  
    //---> Return dtype to user
  return(dtype);

}// End DofType

//****************************************************************************80
void TetElement::FaceData(const intT& k, const intT& p, intT& face, intT& ibub)
{
  //---> Compute the necessary temps
  intT ndofpm1 = CompNdof(p - 1);
  //---> Local index is how many dofs past the p-1 dofs and all 6 edges we are
  intT li = k - ndofpm1 - 6;
  /*---> There are  (p - 2) more dofs per face so face number face is as 
    follows */
  face =(intT)ceil(li/(p - 2));
  //---> Bubble mode is given as 
  ibub = li - (face)*(p - 2) + 1;
  
} // End FaceData

//****************************************************************************80
realT TetElement::TetH1(const intT& k, const realT& xi, const realT& eta, 
			const realT& zeta)
{
  //---> Local Variables
  realT L[4];
  
  //---> Function return variable
  realT phi;
  
  //---> We'll need the vertex functions for any basis we have so specify them
  L[0] = -(1.0 + xi + eta + zeta)*half;
  L[1] = (xi + 1.0)*half;
  L[2] = (eta + 1.0)*half;
  L[3] = (zeta + 1.0)*half;
  
  //---> Get the polynomial order to which this mode belongs
  intT p = mode2p_map(k);
  intT ndofpm1 = CompNdof(p - 1);
  
  if( p == 1 ) {// Order_check
    phi = L[k];
    
  }
  else {
    intT dtype = DofType(k, p);
    
    //---> Compute edge index
    intT eindex = k - ndofpm1;
    switch (dtype) {// DofType
    case 1: //---> Edge dof
	
      switch(eindex) {// Edge_type
      case -1 : 
	//---> Error;
	phi = -9.9e99;
	break;
      case 0 : //---> Edge 0 made of verticies 0 to 1 (1 to 2)
	phi = BarElement::EdgePoly(p, L[0], L[1]);
	break;
      case 1 : //---> Edge 1 of verticies 1 to 2 (2 to 3)
	phi = BarElement::EdgePoly(p, L[1], L[2]);
	break;
      case 2 : //---> Edge 2 or verticies 2 to 0 ( 3 to 1)
	phi = BarElement::EdgePoly(p, L[2], L[0]);
	break;
      case 3 :  //---> Edge 3 of verticies 0 to 3 (1 to 4)
	phi = BarElement::EdgePoly(p, L[0], L[3]);
	break; 
      case 4 :  //---> Edge 4 of verticies 1 to 3 (2 to 4)
	phi = BarElement::EdgePoly(p, L[1], L[3]);
	break; 
      case 5 :  //---> Edge 5 of verticies 2 to 3 (3 to 4)
	phi = BarElement::EdgePoly(p, L[2], L[3]);
	break; 
      } // End Edge_type
      break;
    case 2: //---> Face dof
      //---> Declare 2 face indentifying things we will need
      intT face;
      intT fbub;
	
      //---> Get some face data
      FaceData(k, p, face, fbub);
	
      switch (face) { // Face_type
      case 0: //---> Face 0 : nodes 0, 1, 3 ( 1, 2, 4) 
	phi = TriElement::TriBubPoly(fbub, p, L[0], L[1], L[3]);
	break;
      case 1: //---> Face 1 : nodes 1, 2, 3 ( 2, 3, 4)
	phi = TriElement::TriBubPoly(fbub, p, L[1], L[2], L[3]);
	break;
      case 2: //---> Face 2 : nodes 2, 0, 3 ( 3, 1, 4)
	phi = TriElement::TriBubPoly(fbub, p, L[2], L[0], L[3]);
	break;
      case 3: //---> Face 3 : nodes 0, 2, 1 ( 1, 3, 2)
	phi = TriElement::TriBubPoly(fbub, p, L[0], L[2], L[1]);
	break;
      } // End Face_type

      break;
    case 3: //---> Bubble dof
      //---> Get the local bubble index 
      intT ibub = k - ndofpm1 - 6 - (p - 2)*4  + 1;
      phi = TetBubPoly(ibub, p, L[0], L[1], L[2], L[3]);
	
	
    }// End DofType
      
  }// End order check
    
  //---> Resturn result of fuction 
  return(phi);
    
} // End TetH1
//****************************************************************************80
void TetElement::TetH1D(const intT& k, const realT& xi , const realT& eta, 
			const realT& zeta, realT& dphidxi, realT& dphideta, 
			realT& dphidzeta)
{
  //---> Local Variables
  realT L[4];
  realT dLdxi[4];
  realT dLdeta[4];
  realT dLdzeta[4];
  realT dphidL0;
  realT dphidL1;
  realT dphidL2;
  realT dphidL3;
    
  //---> We'll need the vertex functions for any basis we have so specify them
  L[0] = -(1.0 + xi + eta + zeta)*half;
  L[1] = (xi + 1.0)*half;
  L[2] = (eta + 1.0)*half;
  L[3] = (zeta + 1.0)*half;

  dLdxi[0] = -half;
  dLdxi[1] = half;
  dLdxi[2] = 0.0;
  dLdxi[3] = 0.0;
    
  dLdeta[0] = -half;
  dLdeta[1] = 0.0;
  dLdeta[2] = half;
  dLdeta[3] = 0.0;
    
  dLdzeta[0] = -half;
  dLdzeta[1] = 0.0;
  dLdzeta[2] = 0.0;
  dLdzeta[3] = half;

  //---> Get the discretization order based on specified mode
  intT p = mode2p_map(k);
  intT ndofpm1 = CompNdof(p - 1);
    
  //---> Check the order of the request dof
  if ( p == 1 ) { // Order_check
    //---> For p = 1 the derivatives are really simple since we use these
    //     for constructing other coordinates
    dphidxi = dLdxi[k];
    dphideta = dLdeta[k];
    dphidzeta = dLdzeta[k];
  }
  else {
    //---> Now for every thing more complicated
    dphidxi = 0.0;
    dphideta = 0.0;
    dphidzeta = 0.0;
    
    intT dtype = DofType(k, p);
      
    //---> Compute edge index
    intT eindex = k - ndofpm1;
    switch (dtype) {// DofType
    case 1: //---> Edge dof

      switch(eindex) {// Edge_type
      case -1 : 
	//---> Error;
	dphidxi = -9.9e99;
	dphideta = -9.9e99;
	dphidzeta = -9,9e99;
	break;
      case 0 : //---> Edge 0 made of verticies 0 to 1 (1 to 2)
	BarElement::EdgePolyD(p, L[0], L[1], dphidL0, dphidL1);
	//---> d(phi)/d(xi)
	dphidxi = dphidL0*dLdxi[0] + dphidL1*dLdxi[1];
	//---> d(phi)/d(eta)
	dphideta = dphidL0*dLdeta[0] + dphidL1*dLdeta[1];
	//---> d(phi)/d(zeta)
	dphidzeta = dphidL0*dLdzeta[0] + dphidL1*dLdzeta[1];
	break;
      case 1 : //---> Edge 1 of verticies 1 to 2 (2 to 3)
	BarElement::EdgePolyD(p, L[1], L[2], dphidL1, dphidL2);
	//---> d(phi)/d(xi)
	dphidxi = dphidL1*dLdxi[1] + dphidL2*dLdxi[2];
	//---> d(phi)/d(eta)
	dphideta = dphidL1*dLdeta[1] + dphidL2*dLdeta[2];
	//---> d(phi)/d(zeta)
	dphidzeta = dphidL1*dLdzeta[1] + dphidL2*dLdzeta[2];
	break;
      case 2 : //---> Edge 2 or verticies 2 to 0 ( 3 to 1)
	BarElement::EdgePolyD(p, L[2], L[0], dphidL2, dphidL0);
	//---> d(phi)/d(xi)
	dphidxi = dphidL2*dLdxi[2] + dphidL0*dLdxi[0];
	//---> d(phi)/d(eta)
	dphideta = dphidL2*dLdeta[2] + dphidL0*dLdeta[0];
	//---> d(phi)/d(zeta)
	dphidzeta = dphidL2*dLdzeta[2] + dphidL0*dLdzeta[0];
	break;
      case 3 :  //---> Edge 3 of verticies 0 to 3 (1 to 4)
	BarElement::EdgePolyD(p, L[0], L[3], dphidL0, dphidL3);
	//---> d(phi)/d(xi)
	dphidxi = dphidL0*dLdxi[0] + dphidL3*dLdxi[3];
	//---> d(phi)/d(eta)
	dphideta = dphidL0*dLdeta[0] + dphidL3*dLdeta[3];
	//---> d(phi)/d(zeta)
	dphidzeta = dphidL0*dLdzeta[0] + dphidL3*dLdzeta[3];
	break; 
      case 4 :  //---> Edge 4 of verticies 1 to 3 (2 to 4)
	BarElement::EdgePolyD(p, L[1], L[3], dphidL1, dphidL3);
	//---> d(phi)/d(xi)
	dphidxi = dphidL1*dLdxi[1] + dphidL3*dLdxi[3];
	//---> d(phi)/d(eta)
	dphideta = dphidL1*dLdeta[1] + dphidL3*dLdeta[3];
	//---> d(phi)/d(zeta)
	dphidzeta = dphidL1*dLdzeta[1] + dphidL3*dLdzeta[3];
	break; 
      case 5 :  //---> Edge 5 of verticies 2 to 3 (3 to 4)
	BarElement::EdgePolyD(p, L[2], L[3], dphidL2, dphidL3);
	//---> d(phi)/d(xi)
	dphidxi = dphidL2*dLdxi[2] + dphidL3*dLdxi[3];
	//---> d(phi)/d(eta)
	dphideta = dphidL2*dLdeta[2] + dphidL3*dLdeta[3];
	//---> d(phi)/d(zeta)
	dphidzeta = dphidL2*dLdzeta[2] + dphidL3*dLdzeta[3];
	break; 
      } // End Edge_type
      break;
    case 2: //---> Face dof
      //---> Declare 2 face indentifying things we will need
      intT face;
      intT fbub;
	
      //---> Get some face data
      FaceData(k, p, face, fbub);
	
      switch (face) { // Face_type
      case 0: //---> Face 0 : nodes 0, 1, 3 ( 1, 2, 4) 
	TriElement::TriBubPolyD(fbub, p, L[0], L[1], L[3], 
					     dphidL0, dphidL1, dphidL3);
	//---> d(phi)/d(xi)
	dphidxi = dphidL0*dLdxi[0] + dphidL1*dLdxi[1] + dphidL3*dLdxi[3];
	//---> d(phi)/d(eta)
	dphideta = dphidL0*dLdeta[0] + dphidL1*dLdeta[1] + dphidL3*dLdeta[3];
	//---> d(phi)/d(zeta)
	dphidzeta = dphidL0*dLdzeta[0] + dphidL1*dLdzeta[1] + 
	  dphidL3*dLdzeta[3];
	break;
      case 1: //---> Face 1 : nodes 1, 2, 3 ( 2, 3, 4)
	TriElement::TriBubPolyD(fbub, p, L[1], L[2], L[3],
					     dphidL1, dphidL2, dphidL3);
	//---> d(phi)/d(xi)
	dphidxi = dphidL1*dLdxi[1] + dphidL2*dLdxi[2] + dphidL3*dLdxi[3];
	//---> d(phi)/d(eta)
	dphideta = dphidL1*dLdeta[1] + dphidL2*dLdeta[2] + dphidL3*dLdeta[3];
	//---> d(phi)/d(zeta)
	dphidzeta = dphidL1*dLdzeta[1] + dphidL2*dLdzeta[2] + 
	  dphidL3*dLdzeta[3];
	break;
      case 2: //---> Face 2 : nodes 2, 0, 3 ( 3, 1, 4)
	TriElement::TriBubPolyD(fbub, p, L[2], L[0], L[3], 
					     dphidL2, dphidL0, dphidL3);
	//---> d(phi)/d(xi)
	dphidxi = dphidL2*dLdxi[2] + dphidL0*dLdxi[0] + dphidL3*dLdxi[3];
	//---> d(phi)/d(eta)
	dphideta = dphidL2*dLdeta[2] + dphidL0*dLdeta[0] + dphidL3*dLdeta[3];
	//---> d(phi)/d(zeta)
	dphidzeta = dphidL2*dLdzeta[2] + dphidL0*dLdzeta[0] + 
	  dphidL3*dLdzeta[3];
	break;
      case 3: //---> Face 3 : nodes 0, 2, 1 ( 1, 3, 2)
	TriElement::TriBubPolyD(fbub, p, L[0], L[2], L[1], 
					     dphidL0, dphidL2, dphidL1);
	//---> d(phi)/d(xi)
	dphidxi = dphidL0*dLdxi[0] + dphidL2*dLdxi[2] + dphidL1*dLdxi[1];
	//---> d(phi)/d(eta)
	dphideta = dphidL0*dLdeta[0] + dphidL2*dLdeta[2] + dphidL1*dLdeta[1];
	//---> d(phi)/d(zeta)
	dphidzeta = dphidL0*dLdzeta[0] + dphidL2*dLdzeta[2] + 
	  dphidL1*dLdzeta[1];
	break;
      } // End Face_type

      break;
    case 3: //---> Bubble dof
      //---> Get the local bubble index 
      intT ibub = k - ndofpm1 - 6 - (p - 2)*4  + 1;
      TetBubPolyD(ibub, p, L[0], L[1], L[2], L[3], dphidL0, dphidL1, dphidL2, 
		  dphidL3);
	
      //---> d(phi)/d(xi)
      dphidxi = dphidL0*dLdxi[0] + dphidL1*dLdxi[1] + dphidL2*dLdxi[2] + 
	dphidL3*dLdxi[3];
      //---> d(phi)/d(eta)
      dphideta = dphidL0*dLdeta[0] + dphidL1*dLdeta[1] + dphidL2*dLdeta[2] + 
	dphidL3*dLdeta[3];
      //---> d(phi)/d(zeta)
      dphidzeta = dphidL0*dLdzeta[0] + dphidL1*dLdzeta[1] + 
	dphidL2*dLdzeta[2] + dphidL3*dLdzeta[3];
	
	
    }// End DofType
      
  }// End order check
   
} // End TetH1D

void TetElement::TetGaussPoints(const intT& n, Array2D<realT>& xiq, Array1D<realT>& wq)
{
  //---> Use a switch case block to tablulate the gauss points
  switch (n) { // gauss_point_select
  case 1:
    //---> Order p=1
    //---> Points
    xiq(0,0) = -0.500000000000000;
    xiq(0,1) = -0.500000000000000;
    xiq(0,2) = -0.500000000000000;

    //---> Weights
    wq(0) = 1.333333333333333;

    break;
  case 4:
    //---> Order p=2
    //---> Points 
    xiq(0,0) = -0.723606797749979;
    xiq(0,1) = -0.723606797749979;
    xiq(0,2) = -0.723606797749979;

    xiq(1,0) = 0.170820393249937;
    xiq(1,1) = -0.723606797749979;
    xiq(1,2) = -0.723606797749979;

    xiq(2,0) = -0.723606797749979;
    xiq(2,1) = 0.170820393249937;
    xiq(2,2) = -0.723606797749979;

    xiq(3,0) = -0.723606797749979;
    xiq(3,1) = -0.723606797749979;
    xiq(3,2) = 0.170820393249937;

    //---> Weights
    wq(0) = 0.333333333333333;
    wq(1) = 0.333333333333333;
    wq(2) = 0.333333333333333;
    wq(3) = 0.333333333333333;

    break;
  case 5:
    //---> Order p=3
    //---> Points 
    xiq(0,0) = -0.500000000000000;
    xiq(0,1) = -0.500000000000000;
    xiq(0,2) = -0.500000000000000;

    xiq(1,0) = -0.666666666666667;
    xiq(1,1) = -0.666666666666667;
    xiq(1,2) = -0.666666666666667;

    xiq(2,0) = -0.666666666666667;
    xiq(2,1) = -0.666666666666667;
    xiq(2,2) = 0.000000000000000;

    xiq(3,0) = -0.666666666666667;
    xiq(3,1) = 0.000000000000000;
    xiq(3,2) = -0.666666666666667;

    xiq(4,0) = 0.000000000000000;
    xiq(4,1) = -0.666666666666667;
    xiq(4,2) = -0.666666666666667;

    //---> Weights
    wq(0) = -1.066666666666667;
    wq(1) = 0.600000000000000;
    wq(2) = 0.600000000000000;
    wq(3) = 0.600000000000000;
    wq(4) = 0.600000000000000;

    break;
  case 11:
    //---> Order p=4
    //---> Points 
    xiq(0,0) = -0.500000000000000;
    xiq(0,1) = -0.500000000000000;
    xiq(0,2) = -0.500000000000000;

    xiq(1,0) = -0.857142857142857;
    xiq(1,1) = -0.857142857142857;
    xiq(1,2) = -0.857142857142857;

    xiq(2,0) = -0.857142857142857;
    xiq(2,1) = -0.857142857142857;
    xiq(2,2) = 0.571428571428571;

    xiq(3,0) = -0.857142857142857;
    xiq(3,1) = 0.571428571428571;
    xiq(3,2) = -0.857142857142857;

    xiq(4,0) = 0.571428571428571;
    xiq(4,1) = -0.857142857142857;
    xiq(4,2) = -0.857142857142857;

    xiq(5,0) = -0.201192847666402;
    xiq(5,1) = -0.201192847666402;
    xiq(5,2) = -0.798807152333598;

    xiq(6,0) = -0.201192847666402;
    xiq(6,1) = -0.798807152333598;
    xiq(6,2) = -0.201192847666402;

    xiq(7,0) = -0.798807152333598;
    xiq(7,1) = -0.201192847666402;
    xiq(7,2) = -0.201192847666402;

    xiq(8,0) = -0.201192847666402;
    xiq(8,1) = -0.798807152333598;
    xiq(8,2) = -0.798807152333598;

    xiq(9,0) = -0.798807152333598;
    xiq(9,1) = -0.201192847666402;
    xiq(9,2) = -0.798807152333598;

    xiq(10,0) = -0.798807152333598;
    xiq(10,1) = -0.798807152333598;
    xiq(10,2) = -0.201192847666402;

    //---> Weights
    wq(0) = -0.105244444444444;
    wq(1) = 0.060977777777778;
    wq(2) = 0.060977777777778;
    wq(3) = 0.060977777777778;
    wq(4) = 0.060977777777778;
    wq(5) = 0.199111111111111;
    wq(6) = 0.199111111111111;
    wq(7) = 0.199111111111111;
    wq(8) = 0.199111111111111;
    wq(9) = 0.199111111111111;
    wq(10) = 0.199111111111111;

    break;
  case 14:
    //---> Order p=5
    //---> Points 
    xiq(0,0) = -0.814529499378218;
    xiq(0,1) = -0.814529499378218;
    xiq(0,2) = -0.814529499378218;

    xiq(1,0) = 0.443588498134653;
    xiq(1,1) = -0.814529499378218;
    xiq(1,2) = -0.814529499378218;

    xiq(2,0) = -0.814529499378218;
    xiq(2,1) = 0.443588498134653;
    xiq(2,2) = -0.814529499378218;

    xiq(3,0) = -0.814529499378218;
    xiq(3,1) = -0.814529499378218;
    xiq(3,2) = 0.443588498134653;

    xiq(4,0) = -0.378228161473399;
    xiq(4,1) = -0.378228161473399;
    xiq(4,2) = -0.378228161473399;

    xiq(5,0) = -0.865315515579804;
    xiq(5,1) = -0.378228161473399;
    xiq(5,2) = -0.378228161473399;

    xiq(6,0) = -0.378228161473399;
    xiq(6,1) = -0.865315515579804;
    xiq(6,2) = -0.378228161473399;

    xiq(7,0) = -0.378228161473399;
    xiq(7,1) = -0.378228161473399;
    xiq(7,2) = -0.865315515579804;

    xiq(8,0) = -0.091007408251299;
    xiq(8,1) = -0.091007408251299;
    xiq(8,2) = -0.908992591748701;

    xiq(9,0) = -0.091007408251299;
    xiq(9,1) = -0.908992591748701;
    xiq(9,2) = -0.091007408251299;

    xiq(10,0) = -0.908992591748701;
    xiq(10,1) = -0.091007408251299;
    xiq(10,2) = -0.091007408251299;

    xiq(11,0) = -0.091007408251299;
    xiq(11,1) = -0.908992591748701;
    xiq(11,2) = -0.908992591748701;

    xiq(12,0) = -0.908992591748701;
    xiq(12,1) = -0.091007408251299;
    xiq(12,2) = -0.908992591748701;

    xiq(13,0) = -0.908992591748701;
    xiq(13,1) = -0.908992591748701;
    xiq(13,2) = -0.091007408251299;

    //---> Weights
    wq(0) = 0.097990724155149;
    wq(1) = 0.097990724155149;
    wq(2) = 0.097990724155149;
    wq(3) = 0.097990724155149;
    wq(4) = 0.150250567624021;
    wq(5) = 0.150250567624021;
    wq(6) = 0.150250567624021;
    wq(7) = 0.150250567624021;
    wq(8) = 0.056728027702775;
    wq(9) = 0.056728027702775;
    wq(10) = 0.056728027702775;
    wq(11) = 0.056728027702775;
    wq(12) = 0.056728027702775;
    wq(13) = 0.056728027702775;

    break;
  case 24:
    //---> Order p=6
    //---> Points 
    xiq(0,0) = -0.570794257481696;
    xiq(0,1) = -0.570794257481696;
    xiq(0,2) = -0.570794257481696;

    xiq(1,0) = -0.287617227554912;
    xiq(1,1) = -0.570794257481696;
    xiq(1,2) = -0.570794257481696;

    xiq(2,0) = -0.570794257481696;
    xiq(2,1) = -0.287617227554912;
    xiq(2,2) = -0.570794257481696;

    xiq(3,0) = -0.570794257481696;
    xiq(3,1) = -0.570794257481696;
    xiq(3,2) = -0.287617227554912;

    xiq(4,0) = -0.918652082930777;
    xiq(4,1) = -0.918652082930777;
    xiq(4,2) = -0.918652082930777;

    xiq(5,0) = 0.755956248792332;
    xiq(5,1) = -0.918652082930777;
    xiq(5,2) = -0.918652082930777;

    xiq(6,0) = -0.918652082930777;
    xiq(6,1) = 0.755956248792332;
    xiq(6,2) = -0.918652082930777;

    xiq(7,0) = -0.918652082930777;
    xiq(7,1) = -0.918652082930777;
    xiq(7,2) = 0.755956248792332;

    xiq(8,0) = -0.355324219715449;
    xiq(8,1) = -0.355324219715449;
    xiq(8,2) = -0.355324219715449;

    xiq(9,0) = -0.934027340853653;
    xiq(9,1) = -0.355324219715449;
    xiq(9,2) = -0.355324219715449;

    xiq(10,0) = -0.355324219715449;
    xiq(10,1) = -0.934027340853653;
    xiq(10,2) = -0.355324219715449;

    xiq(11,0) = -0.355324219715449;
    xiq(11,1) = -0.355324219715449;
    xiq(11,2) = -0.934027340853653;

    xiq(12,0) = -0.872677996249965;
    xiq(12,1) = -0.872677996249965;
    xiq(12,2) = -0.460655337083368;

    xiq(13,0) = -0.872677996249965;
    xiq(13,1) = -0.460655337083368;
    xiq(13,2) = -0.872677996249965;

    xiq(14,0) = -0.872677996249965;
    xiq(14,1) = -0.872677996249965;
    xiq(14,2) = 0.206011329583298;

    xiq(15,0) = -0.872677996249965;
    xiq(15,1) = 0.206011329583298;
    xiq(15,2) = -0.872677996249965;

    xiq(16,0) = -0.872677996249965;
    xiq(16,1) = -0.460655337083368;
    xiq(16,2) = 0.206011329583298;

    xiq(17,0) = -0.872677996249965;
    xiq(17,1) = 0.206011329583298;
    xiq(17,2) = -0.460655337083368;

    xiq(18,0) = -0.460655337083368;
    xiq(18,1) = -0.872677996249965;
    xiq(18,2) = -0.872677996249965;

    xiq(19,0) = -0.460655337083368;
    xiq(19,1) = -0.872677996249965;
    xiq(19,2) = 0.206011329583298;

    xiq(20,0) = -0.460655337083368;
    xiq(20,1) = 0.206011329583298;
    xiq(20,2) = -0.872677996249965;

    xiq(21,0) = 0.206011329583298;
    xiq(21,1) = -0.872677996249965;
    xiq(21,2) = -0.460655337083368;

    xiq(22,0) = 0.206011329583298;
    xiq(22,1) = -0.872677996249965;
    xiq(22,2) = -0.872677996249965;

    xiq(23,0) = 0.206011329583298;
    xiq(23,1) = -0.460655337083368;
    xiq(23,2) = -0.872677996249965;

    //---> Weights
    wq(0) = 0.053230333677557;
    wq(1) = 0.053230333677557;
    wq(2) = 0.053230333677557;
    wq(3) = 0.053230333677557;
    wq(4) = 0.013436281407094;
    wq(5) = 0.013436281407094;
    wq(6) = 0.013436281407094;
    wq(7) = 0.013436281407094;
    wq(8) = 0.073809575391540;
    wq(9) = 0.073809575391540;
    wq(10) = 0.073809575391540;
    wq(11) = 0.073809575391540;
    wq(12) = 0.064285714285714;
    wq(13) = 0.064285714285714;
    wq(14) = 0.064285714285714;
    wq(15) = 0.064285714285714;
    wq(16) = 0.064285714285714;
    wq(17) = 0.064285714285714;
    wq(18) = 0.064285714285714;
    wq(19) = 0.064285714285714;
    wq(20) = 0.064285714285714;
    wq(21) = 0.064285714285714;
    wq(22) = 0.064285714285714;
    wq(23) = 0.064285714285714;

    break;

  }// End gauss_point_select

} // End TetGaussPoints

//****************************************************************************80
void TetElement::MapFaceToElem(const intT& side, const realT& u, const realT& v,
			       realT& xi, realT& eta, realT& zeta)
{
  realT xi0;
  realT xi1;
  realT xi2;
    
  realT eta0;
  realT eta1;
  realT eta2;
   
  realT zeta0;
  realT zeta1;
  realT zeta2;
    
  switch(side) { //Pick side
  case 0:
    xi0 = -1.0;
    xi1 =  1.0;
    xi2 = -1.0;
      
    eta0 = -1.0;
    eta1 = -1.0;
    eta2 = -1.0;
      
    zeta0 = -1.0;
    zeta1 = -1.0;
    zeta2 =  1.0;
    break;
  case 1:
    xi0 =  1.0;
    xi1 = -1.0;
    xi2 = -1.0;
      
    eta0 = -1.0;
    eta1 =  1.0;
    eta2 = -1.0;
      
    zeta0 = -1.0;
    zeta1 = -1.0;
    zeta2 =  1.0;
    break;
  case 2:
    xi0 = -1.0;
    xi1 = -1.0;
    xi2 = -1.0;
      
    eta0 =  1.0;
    eta1 = -1.0;
    eta2 = -1.0;
      
    zeta0 = -1.0;
    zeta1 = -1.0;
    zeta2 =  1.0;
    break;
  case 3:
    xi0 = -1.0;
    xi1 = -1.0;
    xi2 =  1.0;
      
    eta0 = -1.0;
    eta1 =  1.0;
    eta2 = -1.0;
      
    zeta0 = -1.0;
    zeta1 = -1.0;
    zeta2 = -1.0;
    break;
  } // End Pick side

    //---> Interpolate using affine coordinates;
  xi = -half*(u + v)*xi0 + half*(u + 1.0)*xi1 + half*(v + 1.0)*xi2;
  eta = -half*(u + v)*eta0 + half*(u + 1.0)*eta1 + half*(v + 1.0)*eta2;
  zeta = -half*(u + v)*zeta0 + half*(u + 1.0)*zeta1 + half*(v + 1.0)*zeta2;
    
} // End MapFaceToElem
//****************************************************************************80
realT TetElement::TetBubPoly(const intT& ibub, const intT& p, const realT& L0, 
			     const realT& L1, const realT& L2, const realT& L3)
{
  //---> Local Variables
  intT np = p - 4 + 1;
  intT k = 0;
  intT n1 = 0;
  intT n2 = 0;
  intT n3 = 0;
  
  //---> Based on value of ibub get the indicies n1, n2, n3;  
  /*---> Ok this is a bit strange...but the index triples are layed out 
    in a triangular pattern of "points".  This triangular patter is such 
    that the max value of any n1,n2,n3 is np with is p - 3.  We can actually 
    represent these values exactly using linear affine(barycentric 
    coodrinates).  We can also rep phyiscal coorindates by i,j and using 
    all this gives us n1, n2, n3.  */
  k = 0;
  for( intT j = 0; j < np; j++) {//i-loop 
    for( intT i = 0; i < np - j; i++) {//i-loop
      k += 1;
      if( k == ibub) {
	n1 = np - i - j;
	n2 = 1 + i;
	n3 = 1 + j;
      }
    } // End i-loop
  } // End j-loop
   
  return(L0*L1*L2*L3*
	 Polynomials::LobattoKern(n1 - 1, L2 - L0)*
	 Polynomials::LobattoKern(n2 - 1, L1 - L0)*
	 Polynomials::LobattoKern(n3 - 1, L3 - L0));
} // End TetBubPoly

//****************************************************************************80
void TetElement::TetBubPolyD(const intT& ibub, const intT& p, const realT& L0, 
			     const realT& L1, const realT& L2, const realT& L3, 
			     realT& dphidL0, realT& dphidL1, realT& dphidL2, 
			     realT& dphidL3) 
{
  //---> Local Variables
  intT np = p - 4 + 1;
  intT k = 0;
  intT n1 = 0;
  intT n2 = 0;
  intT n3 = 0;
  
  //---> Based on value of ibub get the indicies n1, n2, n3;  
    /*---> Ok this is a bit strange...but the index triples are layed out 
      in a triangular pattern of "points".  This triangular patter is such 
      that the max value of any n1,n2,n3 is np with is p - 3.  We can actually 
      represent these values exactly using linear affine(barycentric 
      coodrinates).  We can also rep phyiscal coorindates by i,j and using 
      all this gives us n1, n2, n3.  */
    k = 0;
    for( intT j = 0; j < np; j++) {//i-loop 
      for( intT i = 0; i < np - j; i++) {//i-loop
	k += 1;
	if( k == ibub) {
	  n1 = np - i - j;
	  n2 = 1 + i;
	  n3 = 1 + j;
	}
      } // End i-loop
    } // End j-loop
    
    realT psi1 = Polynomials::LobattoKern(n1 - 1, L2 - L0);
    realT psi2 = Polynomials::LobattoKern(n2 - 1, L1 - L0);
    realT psi3 = Polynomials::LobattoKern(n3 - 1, L3 - L0);

    realT dpsi1 = Polynomials::LobattoKernD(n1 - 1, L2 - L0);
    realT dpsi2 = Polynomials::LobattoKernD(n2 - 1, L1 - L0);
    realT dpsi3 = Polynomials::LobattoKernD(n3 - 1, L3 - L0);

    dphidL0 = L1*L2*L3*psi1*psi2*psi3 + L0*L1*L2*L3*(-dpsi1)*(-dpsi2)*(-dpsi3);
    dphidL1 = L0*L2*L3*psi1*psi2*psi3 + L0*L1*L2*L3*(dpsi2);
    dphidL2 = L0*L1*L3*psi1*psi2*psi3 + L0*L1*L2*L3*(dpsi1);
    dphidL3 = L0*L1*L2*psi1*psi2*psi3 + L0*L1*L2*L3*(dpsi3);

  }// End TetBubPolyD

//****************************************************************************80
realT TetElement::EvalBasis(const intT& n, const realT&  xi, const realT& eta, 
			    const realT& zeta)
{
  //---> Evaluate basis function
  realT phi = TetH1(n, xi, eta, zeta);
  
  //---> Return to user
  return(phi);
  
} // End EvalBasis

//****************************************************************************80
void TetElement::EvalBasisD(const int& n, const realT& xi, const realT& eta, 
			    const realT& zeta, realT& dphidxi, realT& dphideta,
			    realT& dphidzeta)
{
  //---> Evaluate basis function derivatives
  TetH1D(n, xi, eta, zeta, dphidxi, dphideta, dphidzeta);

    return;

} // End EvalBasisD

//****************************************************************************80
void TetElement::initialize(const intT& p_in, const intT& pmap_in, 
			    const intT& deg_in)
{ 
  Element::ndim_ = 3;
  //---> First set p to p_in 
  Element::p_ = p_in;
  Element::pmap_ = pmap_in;

  //---> Now set deg to deg_in
  Element::deg_ = deg_in;

  //---> Compute ndof
  Element::ndof_ = CompNdof(Element::p_);
  
  Element::ndof_map_ = CompNdof(Element::pmap_);
  Element::ndof_face_ = (Element::p_ + 1)*(Element::p_ + 2)/2;
  
  //---> Get number of quadrature points
  Element::nqp_ = p2nqp(Element::deg_);
  Element::nqp_face_ = TriElement::p2nqp(Element::deg_);
    
  //---> Now initialize quadrature vectors
  Element::xiq_.initialize(Element::nqp_, Element::ndim_);
    
  Element::wq_.initialize(Element::nqp_);
  TetGaussPoints(Element::nqp_, Element::xiq_, Element::wq_);
   
  Element::xiq_face_.initialize(Element::nqp_face_, 2);
  Element::wq_face_.initialize(Element::nqp_face_);
    
  TriElement::TriGaussPoints(Element::nqp_face_, 
			     Element::xiq_face_, 
			     Element::wq_face_);
  //---> Now initialize the basis function vectors
  Element::phi_.initialize(Element::nqp_, 
			   Element::ndof_);

  Element::dphi_dxi_.initialize(Element::nqp_, 
				Element::ndim_, 
				Element::ndof_);
    
  Element::phi_face_.initialize(4,
				Element::nqp_face_,
				Element::ndof_);
  Element::dphi_dxi_face_.initialize(4,
				    Element::nqp_face_,
				    Element::ndim_, 
				    Element::ndof_);

  //---> Now initialize the mapping basis functions
  Element::phi_map_.initialize(Element::nqp_,
			       Element::ndof_map_);
  Element::dphi_map_.initialize(Element::nqp_, 
			       Element::ndim_,
			       Element::ndof_map_);
    
  Element::phi_map_face_.initialize(4,
				    Element::nqp_face_,
				    Element::ndof_map_);
  Element::dphi_map_face_.initialize(4,
				     Element::nqp_face_, 
				     Element::ndim_,
				     Element::ndof_map_);
    
  Element::dphi_ds_.initialize(Element::nqp_face_, 
			       Element::ndim_ - 1, 
			       Element::ndof_face_);
  //------------------------- Element Quadrature rules -----------------------
    
  for (intT qp = 0; qp < Element::nqp_; qp++) { //qp_loop 
    for (intT dof = 0; dof < Element::ndof_; dof++) { //dof_loop

      //---> Evalute basis dof at qp
      Element::phi_(qp, dof) = EvalBasis(dof, Element::xiq_(qp,0), 
					 Element::xiq_(qp,1), 
					 Element::xiq_(qp,2) );
	
      //---> Evalute basisD dof: j at qp: i
      EvalBasisD(dof, Element::xiq_(qp,0), Element::xiq_(qp,1), 
		 Element::xiq_(qp,2), 
		 Element::dphi_dxi_(qp,0,dof), 
		 Element::dphi_dxi_(qp,1,dof), 
		 Element::dphi_dxi_(qp,2,dof) );

    }// End dof_loop 

    for (intT dof = 0; dof < Element::ndof_map_; dof++){//dof_map_loop 
      //---> Evaluate mapping basis dof at qp
      Element::phi_map_(qp,dof) = EvalBasis(dof, Element::xiq_(qp,0), 
					   Element::xiq_(qp,1),
					   Element::xiq_(qp,2) );
	
      //---> Evalute mapping basisD dof at qp
      EvalBasisD(dof, Element::xiq_(qp,0), Element::xiq_(qp,1), 
		  Element::xiq_(qp,2),
		  Element::dphi_map_(qp,0,dof), 
		  Element::dphi_map_(qp,1,dof), 
		  Element::dphi_map_(qp,2,dof));
    } // End dof_map_loop 

  }// End qp_loop 
    
  //------------------------- Face Quadrature Rules --------------------------
    
  for(intT side = 0; side < 4; side++){// side_loop
    for (intT qp = 0; qp < Element::nqp_face_; qp++) { //qp_loop
      realT xi, eta, zeta;
      MapFaceToElem(side, Element::xiq_face_(qp,0), 
		    Element::xiq_face_(qp,1), 
		    xi, eta, zeta);
	
      for (intT dof = 0; dof < Element::ndof_; dof++) { //dof_loop
	//---> Evalute basis dof: j at qp:i
	Element::phi_face_(side, qp, dof) = EvalBasis(dof, xi, eta, zeta);
	EvalBasisD(dof, xi, eta, zeta, 
		    Element::dphi_dxi_face_(side, qp, 0, dof),
		    Element::dphi_dxi_face_(side, qp, 1, dof),
		    Element::dphi_dxi_face_(side, qp, 2, dof));

      }// End dof_loop
	
      for (intT dof = 0; dof < Element::ndof_map_; dof++){//dof_map_loop 
	//---> Evaluate mapping basis dof at qp
	Element::phi_map_face_(side, qp, dof) = EvalBasis(dof, xi, eta, zeta);
	EvalBasisD(dof, xi, eta, zeta,
		    Element::dphi_map_face_(side, qp, 0, dof),
		    Element::dphi_map_face_(side, qp, 1, dof),
		    Element::dphi_map_face_(side, qp, 2, dof));

      }// End dof_map_loop
    }// End qp_loop
  }// End side_loop
  //------------------------ Face Parametrization Functions ------------------
  for (intT qp = 0; qp < Element::nqp_face_; qp++) { //qp_loop
    for (intT dof = 0; dof < Element::ndof_face_; dof++) { //dof_loop
      TriElement::EvalBasisD(dof, 
			     Element::xiq_face_(qp,0),
			     Element::xiq_face_(qp,1), 
			     Element::dphi_ds_(qp,0,dof),
			     Element::dphi_ds_(qp,1,dof));
      
    }// End dof_loop
  }// End qp_loop
  
  //---> Face dof map.  In 3-D you have 4 sides and many dofs per side
  Element::face_dof_map_.initialize(4,3); //---> P = 1 only right now;
    
  //---> Side 0;
  Element::face_dof_map_(0,0) = 0;
  Element::face_dof_map_(0,1) = 1;
  Element::face_dof_map_(0,2) = 3;
    
  //---> Side 1;
  Element::face_dof_map_(1,0) = 1;
  Element::face_dof_map_(1,1) = 2;
  Element::face_dof_map_(1,2) = 3;
    
  //---> Side 2;
  Element::face_dof_map_(2,0) = 2;
  Element::face_dof_map_(2,1) = 0;
  Element::face_dof_map_(2,2) = 3;

  //---> Side 3;
  Element::face_dof_map_(3,0) = 0;
  Element::face_dof_map_(3,1) = 2;
  Element::face_dof_map_(3,2) = 1;
    
  return;
}// End initialize
