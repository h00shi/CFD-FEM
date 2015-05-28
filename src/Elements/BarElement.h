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
 
class BarElement : public Element{

private:
  //+++++++++++++++++++++++++++++++ PRIVATE STUFF ++++++++++++++++++++++++++++++
  Array1D<intT> mode2p_map; /*!< Map of mode number of p for 1-D Bar 
				  elements */

//****************************************************************************80
//!
//! \brief CompNdof : Computes the number of dofs for a given order p 
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] p The degree of polynomial approximation
//****************************************************************************80
  intT CompNdof(const intT& p);
 
//****************************************************************************80
//! \brief BarH1 : Computes specified H1 1-D bar function at specified point 
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
  realT BarH1(const intT& k, const realT& xi);
 
//****************************************************************************80
//!
//! \brief BarH1D : Computes the derivative of specified H1 1-D bar function at
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
  realT BarH1D(const intT& k, const realT& xi);
 
protected: 

  Array1D<intT> p2nqp; /*!<Map of number of quadrature points for a 
			     polynomial of degree deg */ 
//****************************************************************************80
//!
//! \brief BarGaussPoints : Tables of the bar gauss points, will populate 
//!                           The gauss points given a number of points
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] n The number of quadrature points
//! \param[out] xiq The coordinates of the quadrature points
//! \param[out] wq The quadrature weights
//****************************************************************************80
  void BarGaussPoints(const intT& n, Array2D<realT>& xiq, Array1D<realT>& wq);
  
//****************************************************************************80
//!
//! \brief EdgePoly : Compute the value of a lobatto polynomial on an edge
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
  realT EdgePoly( const intT& p, const realT& L0, const realT& L1);

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
  void EdgePolyD(const intT& p, const realT& L0, const realT& L1, 
		 realT& dphidL0, realT& dphidL1);
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
  BarElement();
  
//****************************************************************************80
//!
//! \brief ~BarElement : The destructor for the class BarElement
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! 
//****************************************************************************80
  ~BarElement();
  
//****************************************************************************80
//!
//! \brief EvalBasis : A public interface to evaluate a specified polynomial
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] n The polyomial you want to evaluate
//! \param[in] xi The standard element coordinate \f$\xi\f$ where you want the basis function value.
//! \return phi Value of basis functio n \f$ \phi_{n} \left( \xi\right) \f$
//****************************************************************************80
  realT EvalBasis(const intT& n, const realT& xi);
  
//****************************************************************************80
//!
//! \brief EvalBasisD : A public interface to evaluate a the derivate of 
//!                      a specified polynomail. 
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] n The polyomial you want to evaluate
//! \param[in] xi The standard element coordinate \f$\xi\f$ where you want the basis function derivate value.
//! \return phi Value of basis functio n \f$ \frac{ d \phi_{n}}{d \xi} \left( \xi\right) \f$
//****************************************************************************80
  realT EvalBasisD(const intT& n, const realT& xi);
  
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
  void initialize(const intT& p_in, const intT& pmap_in, const intT& deg_in); 
  
}; // End class BarElement

#endif
