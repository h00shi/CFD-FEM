//-*-c++-*-

#ifndef SOLUTION_H
#define SOLUTION_H

//****************************************************************************80
//! \class Solution 
//! \brief This is the header file defining the class Solution
//! \ryan
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//! \tparam intT Template argument meant to mimic integer
//! \tparam realT Template argument meant to mimic real numbers
//****************************************************************************80

#include "Array1D.h"
#include "Array2D.h"
#include "UnstGrid.h"
#include <fstream>


template <typename intT, typename realT >
class Solution {
private:
  
  realT fsmach;
  UnstGrid<intT, realT>& grid; /*!< Ref. to grid on which solution 
					defined */
  Array2D<realT> qnew; /*!<  The latest update to qvector */
  Array2D<realT> qold; /*!< The old q vector */
  Array2D<realT> deltaq; /*!< Newton solver \f$ \Delta q \f$ */
  intT nfld;
  intT pg;

//****************************************************************************80
//!
//! \brief Block default construction of Solution
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! 
//****************************************************************************80
  Solution() : nfld(1), grid(NULL), qnew(1,1), qold(1,1), deltaq(1,1), pg(1) {}
public:

//****************************************************************************80
//!
//! \brief Solution() : 
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \param[in] InGrid Reference to grid on which we will define solution
//! \param[in] nf The input number of fields
//! \param[in] init_field Integer controlling how solution is initialized
//! \param[in] qinit The initial value of q
//! \param[in] p The input value of max polynomial order
//****************************************************************************80
  Solution(UnstGrid<intT,realT>& InGrid, const intT& nf, 
	   const intT& init_field, const Array1D<realT>& qinit, const intT& p) :
    nfld(nf), grid(InGrid), qnew(InGrid.nnode, nf), qold(InGrid.nnode, nf), 
    deltaq(InGrid.nnode, nf), pg(p)
  {
 
    switch (init_field) {
    case 0: //---> Initialize to free-stream values 
      
      for (int n = 0; n < grid.nnode; n++) { // node loop
        for (int f = 0; f < nfld; f++) {
	  qnew(n,f) = qinit(f);
	  qold(n,f) = qinit(f);
	  deltaq(n,f) = 0.0;
        }
      }

      break;
    }

  }

//****************************************************************************80
//!
//! \brief ~Solution() : Destructor for solution class
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! 
//****************************************************************************80
  ~Solution()
  {
    /*---> No need to freem memory, all components are based on classes that
      have pointers delete in destructors*/
    
  }

//****************************************************************************80
//!
//! \brief get_qnew : Returns a component of qnew
//! \details
//! \ryan
//! \version $Rev$
//! \date $Date$
//! \param[in] dof The degree of freedom we want
//! \param[in] fld The field we want 
//! \return qnew(dof,fld) Value of qnew for (dof, fld)
//****************************************************************************80
  realT get_qnew(const intT& dof, const intT& fld) {
    return qnew(dof, fld);
  } // End get_qnew

//****************************************************************************80
//!
//! \brief get_qold : Returns a component of qold
//! \details
//! \ryan
//! \version $Rev$
//! \date $Date$
//! \param[in] dof The degree of freedom we want
//! \param[in] fld The field we want 
//! \return qold(dof,fld) Value of qnew for (dof, fld)
//****************************************************************************80
  realT get_qold(const intT& dof, const intT& fld) {
    return qold(dof, fld);
  }

//****************************************************************************80
//!
//! \brief get_deltaq : Returns a component of deltaq
//! \details
//! \ryan
//! \version $Rev$
//! \date $Date$
//! \param[in] dof The degree of freedom we want
//! \param[in] fld The field we want 
//! \return qnew(dof,fld) Value of qnew for (dof, fld)
//****************************************************************************80
  realT get_deltaq(const intT& dof, const intT& fld) {
    return deltaq(dof, fld);
  }
  
  // realT updateq(const intT& i, const intT& j) {
  //   qnew(i,j) += deltaq(i,j);
  //   return qnew(i,j);
  // }
  
//****************************************************************************80
//!
//! \brief set_qnew : Sets qnew for dof and field 
//! \details
//! \ryan
//! \version $Rev$
//! \date $Date$
//! \param[in] dof The degree of freedom we want
//! \param[in] fld The field we want 
//! \param[in] val The value value to set
//****************************************************************************80
  void set_qnew(const intT& dof, const intT& fld, const realT& val) {
    qnew(dof, fld) = val;
  }

//****************************************************************************80
//!
//! \brief set_qold : Sets qold for dof and field 
//! \details
//! \ryan
//! \version $Rev$
//! \date $Date$
//! \param[in] dof The degree of freedom we want
//! \param[in] fld The field we want 
//! \param[in] val The value value to set
//****************************************************************************80
  void set_qold(const intT& dof, const intT& fld, const realT& val) {
    qold(dof, fld) = val;
  }

//****************************************************************************80
//!
//! \brief set_deltaq : Sets qold for dof and field 
//! \details
//! \ryan
//! \version $Rev$
//! \date $Date$
//! \param[in] dof The degree of freedom we want
//! \param[in] fld The field we want 
//! \param[in] val The value value to set
//****************************************************************************80
  void set_deltaq(const intT& dof, const intT& fld, const realT& val) {
    deltaq(dof, fld) = val;
  }

//****************************************************************************80
//!
//! \brief write_tecplot : Writes the solution to a tecplot file
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! 
//****************************************************************************80

  void write_tecplot(std::string& Project) {
    std::ofstream tecfile;
    
    tecfile.open( (char*)(Project + std::string("_soln.tec")).c_str() );
    tecfile << "Title = Computational Mesh" << std::endl;
    tecfile << "FILETYPE=FULL" << std::endl;
    switch (grid.ndim){
    case 1:
      tecfile << "Variables = X";
      break;
    case 2:
      tecfile << "Variables = X, Y";
      break;
    case 3:
      tecfile << "Variables = X, Y, Z";
      break;
    }
    for(int f = 0; f < nfld; f++) {
      std::string field;
      std::ostringstream converter;
      converter << f;
      field = converter.str();
      tecfile << std::string(", Q-").c_str() + field; 
    }
    tecfile << std::endl;

    tecfile << "ZONE T = Solution" << std::endl;
    switch (grid.ndim) {
    case 1:
      tecfile << "ZONETYPE=FELINESEG" << std::endl;
      break;
    case 2:
      tecfile << "ZONETYPE=FETRIANGLE" << std::endl;
      break;
    case 3:
      tecfile << "ZONETYPE=FETETRAHEDRON" << std::endl;
      break;
    }
    tecfile << "NODES = " << grid.nnode << std::endl;
    tecfile << "ELEMENTS=" << grid.nelement << std::endl;
    tecfile << "DATAPACKING=POINT" << std::endl;
    for(intT n = 0; n < grid.nnode; n++){// Node loop 
      
      for(intT d = 0; d < grid.ndim; d++){ // Dimension loop 
	tecfile << grid.get_node_coord(n, d) << " ";
      } // End Node loop
      
      for(intT f = 0; f < nfld; f++){ // Field loop 
	tecfile << get_qnew(n,f) << " ";  
      }// End Field loop
      
      tecfile << std::endl;
    } // End Dimension loop 

    for(intT e = 0; e < grid.nelement; e++){// Element loop 
      for(intT n = 0; n < grid.get_nnode_on_element(e); n++) { // Node loop 
	tecfile << grid.get_node_on_element(e,n)+1 << " ";
      } // End Node loop 
      tecfile << std::endl;
    } // End element loop
    tecfile.close();

  }// End write_tecplot

//****************************************************************************80
//!
//! \brief get_pg : Returns the value of pg, the global order for the solution
//!                 instance.
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! \return pg Global order for this solution instance
//****************************************************************************80
  intT get_pg() {return(pg);}

};  // end Solution class
#endif

