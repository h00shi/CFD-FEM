// -*-c++-*-
#ifndef TRIELEMENT_H
#define TRIELEMENT_H
//****************************************************************************80
//! \file TriElement.h
//! \class TriElement TriElement.h
//! \brief This the header file defining the class TriElement, which defines 
//!  operators and data for the a 2-D triangular element.
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! 
//****************************************************************************80

#include "my_incl.h"
#include "consts.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/Array3D.h"
#include "BasisFunctions/Polynomials.h"
#include "Elements/BarElement.h"

class TriElement : public BarElement {

 private:
 //+++++++++++++++++++++++++++++++ PRIVATE STUFF ++++++++++++++++++++++++++++++
  Array1D<intT>  mode2p_map; /*!< Map of mode number of p for 2-D Tri elements
			      */
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
//!
//! \brief DofType : A function to determine what type of heirarchical 
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
  intT DofType(const intT& k, const intT& p);
  
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
  realT TriH1(const intT& k, const realT& xi, const realT& eta);
 
//****************************************************************************80
//!
//! \brief TriH1D : Derivative of H1 triangular basis function
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
  void TriH1D(const intT& k, const realT& xi, const realT& eta, 
	      realT& dphidxi, realT& dphideta);
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
  void TriGaussPoints(const intT& n, Array2D<realT>& xiq, Array1D<realT>& wq);
 
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
  void MapFaceToElem(const intT& side, const realT& u, realT& xi, realT& eta);
  
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
  realT TriBubPoly(const intT& ibub, const intT& p, const realT& L0, 
		   const realT& L1, const realT& L2);


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
  void TriBubPolyD( const intT& ibub, const intT& p, const realT& L0, 
		    const realT& L1, const realT& L2, realT& dphidL0, 
		    realT& dphidL1, realT& dphidL2);
    
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
  TriElement();

//****************************************************************************80
//!
//! \brief ~TriElement : Destructor for class TriElement
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! 
//****************************************************************************80
  virtual ~TriElement();

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
  realT EvalBasis(const intT& n, const realT&  xi, const realT& eta); 

//****************************************************************************80
//!
//! \brief A function to evalute the selected basis function derivative
//! \details
//! \nick 
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! 
//****************************************************************************80
  void EvalBasisD(const int& n, const realT& xi, const realT& eta, 
		  realT& dphidxi, realT& dphideta);
  
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
  void initialize(const intT& p_in, const intT& pmap_in, const intT& deg_in); 
   
}; // End class TriElement

#endif
