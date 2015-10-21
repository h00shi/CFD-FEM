//-*-c++-*-
#ifndef CGRESIDUAL_H
#define CGRESIDUAL_H

#include "my_incl.h"
#include "consts.h"
#include "Equation.h"
#include "Array1D.h"
#include "Array2D.h"
#include "Array3D.h"
#include "CGMesh.h"
#include "Solution.h"
#include "Element.h"
#include "BarElement.h"
#include "TriElement.h"
#include "TetElement.h"
#include "system_module.h"

//****************************************************************************80
//!
//! \brief A class to compute Galerkin Least Squares Residuals
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \tparam intT Integter data type
//! \tparam realT Real data type
//! \tparam PDET Partial differential equation type 
//****************************************************************************80
template<class PDET>
class CGResidual {

private:
  PDET& PDE; /*!< Reference to PDE class already instantiated */
  Solution<intT,realT>& soln; /*!< Reference to solution class already 
				instantiated. */
  CGMesh& grid; /*!< Reference to grid class already
				instantiated. */
  
//****************************************************************************80
//!
//! \brief  Solver : Default constructor
//! \details This is blocked...so you can't screw up 
//! \nick 
//! \version $Rev$
//! \date $Date$
//!
//****************************************************************************80
  CGResidual() = delete;
public: 

//****************************************************************************80
//!
//! \brief CGResidual : Constructor taking PDE, grid and soln references
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] PDE_ref The reference to the PDE class
//! \param[in] grid_ref The reference to the grid class
//! \param[in] soln_ref The reference to the soln class
//****************************************************************************80
  CGResidual(PDET & PDE_ref, CGMesh& grid_ref) :
      PDE(PDE_ref), grid(grid_ref)
  {
    //---> Set global order local variable
    intT pg = soln.get_pg();
    
    //---> Instantiate Elements
    BarElement<intT, realT> bar;
    bar.initialize(pg, pg, 2*pg);
    TriElement<intT, realT> tri;
    tri.initialize(pg, pg, 2*pg);
    TetElement<intT, realT> tet;
    tet.initialize(pg, pg, 2*pg);
    
    
    //---> Setup the grid for use in computing the residual
    for (intT e = 0; e < grid.nelement; e++){ // element loop
     
      Array2D<realT> xmap(PDE.get_ndim(), grid.get_nnode_on_element(e));
      
      for (intT d = 0; d < PDE.get_ndim(); d++){ // xmap
        for (intT n = 0; n < grid.get_nnode_on_element(e); n++) {
          xmap(d, n) = grid.get_node_coord(grid.get_node_on_element(e, n), d);
        }
      } // End xmap
      

    
    }// End element loop 

  } // End CGResidual

//****************************************************************************80
//!
//! \brief CGResid : Computes and returns the volume residual
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[out] resid The reference to the residual array
//****************************************************************************80
  void CalcResidual(Array2D<realT>& resid)
  {
    //---> Set global order local variable
    intT pg = soln.get_pg();
    
    //---> Instantiate Elements
    BarElement<intT, realT> bar;
    bar.initialize(pg, pg, 2*pg);
    TriElement<intT, realT> tri;
    tri.initialize(pg, pg, 2*pg);
    TetElement<intT, realT> tet;
    tet.initialize(pg, pg, 2*pg);
    resid.set_value(0.0);


    //------------------------------- Volume Loop ------------------------------
    for (intT e = 0; e < grid.nelement; e++) { // Element Resid Loop
      //---> Local Element Data
      Array2D<realT> xmap(PDE.get_ndim(), grid.get_nnode_on_element(e));
      Array2D<realT> qcoeff(PDE.get_nfld(), grid.get_nnode_on_element(e));
      Array2D<realT> rescoeff(PDE.get_nfld(), 
				    grid.get_nnode_on_element(e));
      
      //---> Fill Local Variables
      for (intT f = 0; f < PDE.get_nfld(); f++) { // qcoeff
	for (intT n = 0; n < grid.get_nnode_on_element(e); n++) {  
	  qcoeff(f, n) = soln.get_qnew(grid.get_node_on_element(e, n), f);
	} 
      }

      for (intT d = 0; d < PDE.get_ndim(); d++){ // xmap
	for (intT n = 0; n < grid.get_nnode_on_element(e); n++) {  
	  xmap(d, n) = grid.get_node_coord(grid.get_node_on_element(e, n), d);
	} 
      } // End xmap
      
      switch(grid.get_element_type(e)) {
      case 0:
	bar.template CGVolResid<PDET>(PDE, qcoeff, xmap, rescoeff);
	break;
      case 1:
	tri.template CGVolResid<PDET>(PDE, qcoeff, xmap, rescoeff);
	break;
      case 3:
	tet.template CGVolResid<PDET>(PDE, qcoeff, xmap, rescoeff);
	break;
      }
     
      //---> Fill push residual to dof's
      for (intT f = 0; f < PDE.get_nfld(); f++) { // qcoeff
	for (intT n = 0; n < grid.get_nnode_on_element(e); n++) {  
	  resid(grid.get_node_on_element(e, n), f) += rescoeff(f,n);
	} 
      }
    
    } // End Element Resid Loop 
  
    //------------------------------ BC Face Loop ------------------------------
    for (intT bcf = 0; bcf < grid.nbc_face; bcf++){ // BC Face Loop 
      //---> Get the element attached to the boundary face
      intT e = grid.get_element_on_bc_face(bcf);
      //---> Which "side" of element e does this face correspond to?
      intT side = grid.get_bc_local_face(bcf);
      //---> Get boundary condition type
      intT bctype = grid.get_bcid_type(grid.get_bc_face_id(bcf));
      
      //---> Local Element Data
      Array2D<realT> xmap(PDE.get_ndim(), grid.get_nnode_on_element(e));
      Array2D<realT> qcoeff(PDE.get_nfld(), grid.get_nnode_on_element(e));
      Array2D<realT> rescoeff(PDE.get_nfld(), 
				    grid.get_nnode_on_element(e));
        
      //---> Fill Local Variables
      for (intT f = 0; f < PDE.get_nfld(); f++) { // qcoeff
	for (intT n = 0; n < grid.get_nnode_on_element(e); n++) {  
	  qcoeff(f, n) = soln.get_qnew(grid.get_node_on_element(e, n), f);
	} 
      }

      for (intT d = 0; d < PDE.get_ndim(); d++){ // xmap
	for (intT n = 0; n < grid.get_nnode_on_element(e); n++) {  
	  xmap(d, n) = grid.get_node_coord(grid.get_node_on_element(e, n), d);
	} 
      } // End xmap
    
      switch(grid.get_element_type(e)) {
      case 0:
	bar.template CGBCResid<PDET>(PDE, bctype, side, qcoeff, xmap, 
				     rescoeff);
	break;
      case 1:
	tri.template CGBCResid<PDET>(PDE, bctype, side, qcoeff, xmap, 
				     rescoeff);
	break;
      case 3:
	tet.template CGBCResid<PDET>(PDE, bctype, side, qcoeff, xmap, 
				     rescoeff);
	break;
      }

      //---> Fill push residual to dof's
      for (intT f = 0; f < PDE.get_nfld(); f++) { // Field loop
	for (intT n = 0; n < grid.get_nnode_on_element(e); n++) { // Dof loop 
	  resid(grid.get_node_on_element(e, n), f) += rescoeff(f,n);
	} // End Dof loop 
      } // End field loop
    
    }// End BC Face Loop 

  } // End CGResid
  
//****************************************************************************80
//!
//! \brief CGJac : Computes the Galerkin Least Square discretization
//!                 residual and Jacobian.
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[out] resid The residual
//! \param[out] jac The jacobian
//****************************************************************************80
  void CGJac(Array2D<realT>& resid)
  {
  //---> Set global order local variable
    intT pg = soln.get_pg();
    
    //---> Instantiate Elements
    BarElement<intT, realT> bar;
    bar.initialize(pg, pg, 2*pg);
    TriElement<intT, realT> tri;
    tri.initialize(pg, pg, 2*pg);
    TetElement<intT, realT> tet;
    tet.initialize(pg, pg, 2*pg);
    resid.set_value(0.0);
    //------------------------------- Volume Loop ------------------------------
    for (intT e = 0; e < grid.nelement; e++) { // Element Resid Loop
      //---> Local Element Data
      Array2D<realT> xmap(PDE.get_ndim(), grid.get_nnode_on_element(e));
      Array2D<realT> qcoeff(PDE.get_nfld(), grid.get_nnode_on_element(e));
      Array2D<realT> rescoeff(PDE.get_nfld(), 
				    grid.get_nnode_on_element(e));
      Array4D<realT> jaccoeff(PDE.get_nfld(), PDE.get_nfld(), 
			      grid.get_nnode_on_element(e),
			      grid.get_nnode_on_element(e));
      
      //---> Fill Local Variables
      for (intT f = 0; f < PDE.get_nfld(); f++) { // qcoeff
	for (intT n = 0; n < grid.get_nnode_on_element(e); n++) {  
	  qcoeff(f, n) = soln.get_qnew(grid.get_node_on_element(e, n), f);
	} 
      }

      for (intT d = 0; d < PDE.get_ndim(); d++){ // xmap
	for (intT n = 0; n < grid.get_nnode_on_element(e); n++) {  
	  xmap(d, n) = grid.get_node_coord(grid.get_node_on_element(e, n), d);
	} 
      } // End xmap
      
      switch(grid.get_element_type(e)) {
      case 0:
	bar.template CGVolJac<PDET>(PDE, qcoeff, xmap, rescoeff, jaccoeff);
	break;
      case 1:
	tri.template CGVolJac<PDET>(PDE, qcoeff, xmap, rescoeff, jaccoeff);
	break;
      case 3:
	tet.template CGVolJac<PDET>(PDE, qcoeff, xmap, rescoeff, jaccoeff);
	break;
      }
     
      //---> Fill push residual to dof's
      for (intT f = 0; f < PDE.get_nfld(); f++) { // qcoeff
	for (intT n = 0; n < grid.get_nnode_on_element(e); n++) {  
	  resid(grid.get_node_on_element(e, n), f) += rescoeff(f,n);
	} 
      }
    
    } // End Element Resid Loop 
  
    //------------------------------ BC Face Loop ------------------------------
    for (intT bcf = 0; bcf < grid.nbc_face; bcf++){ // BC Face Loop 
      //---> Get the element attached to the boundary face
      intT e = grid.get_element_on_bc_face(bcf);
      //---> Which "side" of element e does this face correspond to?
      intT side = grid.get_bc_local_face(bcf);
      //---> Get boundary condition type
      intT bctype = grid.get_bcid_type(grid.get_bc_face_id(bcf));
      
      //---> Local Element Data
      Array2D<realT> xmap(PDE.get_ndim(), grid.get_nnode_on_element(e));
      Array2D<realT> qcoeff(PDE.get_nfld(), grid.get_nnode_on_element(e));
      Array2D<realT> rescoeff(PDE.get_nfld(), 
				    grid.get_nnode_on_element(e));
      Array4D<realT> jaccoeff(PDE.get_nfld(), PDE.get_nfld(), 
			      grid.get_nnode_on_element(e),
			      grid.get_nnode_on_element(e));
      //---> Fill Local Variables
      for (intT f = 0; f < PDE.get_nfld(); f++) { // qcoeff
	for (intT n = 0; n < grid.get_nnode_on_element(e); n++) {  
	  qcoeff(f, n) = soln.get_qnew(grid.get_node_on_element(e, n), f);
	} 
      }

      for (intT d = 0; d < PDE.get_ndim(); d++){ // xmap
	for (intT n = 0; n < grid.get_nnode_on_element(e); n++) {  
	  xmap(d, n) = grid.get_node_coord(grid.get_node_on_element(e, n), d);
	} 
      } // End xmap
    
      switch(grid.get_element_type(e)) {
      case 0:
	bar.template CGBCJac<PDET>(PDE, bctype, side, qcoeff, xmap, 
				   rescoeff, jaccoeff);
	break;
      case 1:
	tri.template CGBCJac<PDET>(PDE, bctype, side, qcoeff, xmap, 
				   rescoeff, jaccoeff);
	break;
      case 3:
	tet.template CGBCJac<PDET>(PDE, bctype, side, qcoeff, xmap, 
				   rescoeff, jaccoeff);
	break;
      }

      //---> Fill push residual to dof's
      for (intT f = 0; f < PDE.get_nfld(); f++) { // Field loop
	for (intT n = 0; n < grid.get_nnode_on_element(e); n++) { // Dof loop 
	  resid(grid.get_node_on_element(e, n), f) += rescoeff(f,n);
	} // End Dof loop 
      } // End field loop
    
    }// End BC Face Loop 
    

  } //End CGJac

//****************************************************************************80
//!
//! \brief l2norm : Computes the vector-L2 norm of the residual 
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] resid The residual stored as an Array
//****************************************************************************80
  realT l2norm(Array2D<realT>& resid) 
  {
    //---> Initialize the norm to 0.0
    realT norm = 0.0;
    
    for (intT n = 0; n < grid.nnode; n++){ // element loop 
      for (intT f = 0; f < PDE.get_nfld(); f++) { // Field loop 
	norm += resid(n,f)*resid(n,f);
      }// End Field loop 
    } // End element loop 
    
    return(sqrt(norm));
    
  } // End l2norm

//****************************************************************************80
//!
//! \brief infnorm : Computes and infinity norm of the residual 
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] resid The residual stored as an Array
//****************************************************************************80
  realT infnorm(Array2D<realT>& resid)
  {
    //---> Initialize the norm to 0.0
    realT rmax = 0.0;
    
    for (intT n = 0; n < grid.nnode; n++){ // element loop 
      for (intT f = 0; f < PDE.get_nfld(); f++) { // Field loop 
	rmax = max(rmax, abs(resid(n,f)) ) ;
      }// End Field loop 
    } // End element loop 

    return(rmax);
  }// End infnorm
    
//****************************************************************************80
//!
//! \brief check_nan : Runs diagnostics for NaN in the residual
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] resid The residual vector stored as Array2D
//****************************************************************************80
  void check_nan(const Array2D<realT>& resid) 
  {
    //---> Loop over nodes
    for(intT n = 0; n < grid.nnode; n++) { // node loop 
      for(intT f = 0; f < PDE.get_nfld(); f++){// Field loop
	if( resid(n,f) != resid(n,f)) { 
	  std::cout << "NaN found at Node: " << n << " Field: " 
		    << f << std::endl;
	  //std::cout << "The Physical Location is" << grid.get_

	  system_module::pause();
	}
      }// End Field loop 
    }// End node loop 

  } // End check_nan


};

#endif

