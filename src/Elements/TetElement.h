// -*-c++-*-
#ifndef TETELEMENT_H
#define TETELEMENT_H
//****************************************************************************80
//! \file TetElement.h
//! \class TetElement TetElement.h
//! \brief This the header file defining the class TetElement, which defines
//!  operators and data for the a 3-D tetrahedral element.
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! 
//****************************************************************************80

#include "my_incl.h"
#include "consts.h"
#include "Array1D.h"
#include "Array2D.h"
#include "Array3D.h"
#include "Polynomials.h"
#include "TriElement.h"
 
class TetElement : public TriElement{

 private:
 //+++++++++++++++++++++++++++++++ PRIVATE STUFF ++++++++++++++++++++++++++++++
  Array1D<intT>  mode2p_map; /*!< Map of mode number of p for 3-D Tet elements 
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
//! \brief FaceData : Obtains the face id number and bub index corresponding 
//!                   to dof k of degree p
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] k The dof number we are working with
//! \param[in] p The degree of polynomial
//! \param[out] face The face number 0 to 3;
//! \param[out] ibub The face bubble index starting at 1 
//****************************************************************************80
  void FaceData(const intT& k, const intT& p, intT& face, intT& ibub);

//****************************************************************************80
//!
//! \brief TetH1 : Computes the value of a specified H1 basis function at
//!                a specified point \f$(\xi, \eta, \zeta) \f$ in a tetrahedra. 
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] k The basis function number
//! \param[in] xi The \f$ \xi \f$ coordinate \f$ \xi \in [-1,1] \f$ 
//! \param[in] eta The \f$ \eta \f$ coordinate  \f$ \eta \in [-1,1] \f$ 
//! \param[in] zeta The \f$ \zeta \f$ coorindate \f$ \zeta \in [-1,1] \f$
//! \return phi The value of the basis fucion \f$ \phi(\xi,\eta,\zeta) \f$
//****************************************************************************80
  realT TetH1(const intT& k, const realT& xi, const realT& eta, 
	      const realT& zeta);

//****************************************************************************80
//!
//! \brief TetH1D : Computes value of the derivative of a speciefied H1 
//!                 basis function at a specified point \f$(\xi,\eta,\zeta)\f$
//!                 in a Tetrahedra. 
//! \details
//! \nick 
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \param[in] k The basis function number we are evaluating
//! \param[in] xi \f$ \xi \f$ coordinate
//! \param[in] eta \f$ \eta \f$ coordinate 
//! \param[in] zeta \f$ \zeta \f$ coordinate
//! \param[out] dphidxi \f$ \frac{\partial \phi_{k}}{\partial \xi} \f$
//! \param[out] dphideta \f$ \frac{\partial \phi_{k}}{\partial \eta} \f$
//! \param[out] dphidzeta \f$ \frac{\partial \phi_{k}}{\partial \zeta} \f$
//****************************************************************************80
  void TetH1D(const intT& k, const realT& xi , const realT& eta, 
	      const realT& zeta, realT& dphidxi, realT& dphideta, 
	      realT& dphidzeta);

protected:
  
  Array1D<intT> p2nqp; /*!<Map of number of quadrature points for a 
                             polynomial of degree deg */

//****************************************************************************80
//!
//! \brief Tables of the Tetrahedral gauss points. Populates base class arrays
//!        given a number of points
//! \details
//! \ryan
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] n The number of points required
//! \param[out] xiq The quadrature point coordinates
//! \param[out] wq The quadrature weights
//****************************************************************************80
  void TetGaussPoints(const intT& n, Array2D<realT>& xiq, Array1D<realT>& wq);

//****************************************************************************80
//!
//! \brief map_face_to_elem : Maps the face point coordinats to their element
//! coordinates
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] side The side of the element we are on
//! \param[in] u The face coordinate for the triangle 
//! \param[in] v The face coordinate for the triangle
//! \param[out] xi The xi-direction coordinate in the tetrahedra
//! \param[out] eta The eta-direction coordinate in the tetrahedra
//! \param[out] zeta The zeta-direction coordinate in the tetrahedra
//****************************************************************************80
  void MapFaceToElem(const intT& side, const realT& u, const realT& v, 
		     realT& xi, realT& eta, realT& zeta);
    
//****************************************************************************80
//!
//! \brief TetBubPoly : Tetrahedral bubble polynomial
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] ibub The bubble index
//! \param[in] p The degree of polynomial
//! \param[in] L0 The first affine coordinate
//! \param[in] L1 The second affine coordinate
//! \param[in] L2 The third affine coordinate
//! \param[in] L3 The fourth affine coordinate
//! \return phi The bubble basis function
//****************************************************************************80
  realT TetBubPoly(const intT& ibub, const intT& p, const realT& L0, 
		   const realT& L1, const realT& L2, const realT& L3);
 
//****************************************************************************80
//!
//! \brief TetBubPolyD : Derivative of tetrahedral bubble functions w.r.t. 
//!        affine coordinates
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] ibub The bubble index
//! \param[in] p The degree of polynomial
//! \param[in] L0 The first affine coordinate
//! \param[in] L1 The second affine coordinate
//! \param[in] L2 The third affine coordinate
//! \param[in] L3 The fourth affine coordinate
//! \param[out] dphidL0 \f$ \frac{ \partial \phi}{\partial L_{0}} \f$
//! \param[out] dphidL1 \f$ \frac{ \partial \phi}{\partial L_{1}} \f$
//! \param[out] dphidL2 \f$ \frac{ \partial \phi}{\partial L_{2}} \f$
//! \param[out] dphidL3 \f$ \frac{ \partial \phi}{\partial L_{3}} \f$
//****************************************************************************80
  void TetBubPolyD(const intT& ibub, const intT& p, const realT& L0, 
		   const realT& L1, const realT& L2, const realT& L3, 
		   realT& dphidL0, realT& dphidL1, realT& dphidL2, 
		   realT& dphidL3);

public: 
//****************************************************************************80
//!
//! \brief TetElement : The constructor for class tet-element
//! \details
//! \ryan
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//****************************************************************************80
  TetElement();		 
 
//****************************************************************************80
//!
//! \brief ~TetElement : Destructor for class TetElement
//! \details
//! \ryan
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! 
//****************************************************************************80
  virtual ~TetElement();
 
//****************************************************************************80
//! 
//! \brief A function to evaluate the selected basis function. 
//! \details
//! \ryan
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] n The basis function number you want
//! \param[in] xi The standard element coordinates \f$ \xi\f$
//! \param[in] eta The standard element coordinates \f$ \eta \f$
//! \param[in] zeta The standard element coordinates \f$ \zeta \f$
//! \return phi Value of basis function n \f$ \phi_{n}\left( 
//! (\xi,\eta,\zeta)\right) \f$
//****************************************************************************80
  realT EvalBasis(const intT& n, const realT&  xi, const realT& eta, 
		  const realT& zeta);
  
//****************************************************************************80
//!
//! \brief A function to evalute the selected basis function derivative
//! \details
//! \ryan
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! 
//****************************************************************************80
  void EvalBasisD(const int& n, const realT& xi, const realT& eta, 
		  const realT& zeta, realT& dphidxi, realT& dphideta, 
		  realT& dphidzeta);

//****************************************************************************80
//!
//! \brief Function to initialize an instance of the class
//! \details
//! \ryan
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] p_in The order of polynomials you want to integrate
//! \param[in] pmap_in The order of mapping polynomials
//! \param[in] deg_in The degree of integration required, recommend > 2*p_in.
//****************************************************************************80
  void initialize(const intT& p_in, const intT& pmap_in, const intT& deg_in);

}; // End of class TetElement
#endif
