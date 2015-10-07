//-*-c++-*-
#ifndef MKLDSS_H
#define MKLDSS_H

#include "my_incl.h"
#include "SystemUtils/SystemModule.h"
#include "SparseMatrix/CSRMatrix.h"
#ifdef MKL_DSS
#include "mkl_types.h"
#include "mkl_dss.h"
#endif

//****************************************************************************80
//! \class MKLDSS
//! \brief MKLDSS : Provides interfaces to MKL direct sparse solver
//!        external libraries.
//!
//! \nick
//! \version $Rev: 5 $
//! \tparam intT Template parameter for type of integer
//! \tparam dataT Template parameter for type of data_
//****************************************************************************80
class MKLDSS{

public:
//****************************************************************************80
//! \brief 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  MKLDSS(CSRMatrix<realT>& matrix);

//****************************************************************************80
//! \brief ~MKLDSS : Constructor for this class.  Requires
//!        number of rows and total number of non-zero entries in CSR Matrix.
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[in] nrow The number of rows in the sparse matrix
//! \param[in] nnz The number of of non-zero entries in the sparse matrix
//****************************************************************************80
  ~MKLDSS();
   
//****************************************************************************80
//! \brief MKLDSSFactorize : For objects of type MatrixCSR perform a sparese LU
//!                    LU factorization using Intel MKL Direct Sparse Solver
//!                    (DSS)
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//! \return Returns boolean true if sucessful, false if unsucessful
//****************************************************************************80
  bool Factorize();
 
//****************************************************************************80
//! \brief  MKLDSSSolve: Given that the matrix has already been factorized...
//!           solve using MKL DSS Solve
//! \details  Solves [A]{x}={b} using MKL DSSS
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[in] b The right hand side or b vector according to the equation
//! \param[out] x The solution vector
//****************************************************************************80
  void Solve(const List2D<double>& b, List2D<double>& x);
  
private:
  CSRMatrix<realT>& matrix_;
#ifdef MKL_DSS
  _MKL_DSS_HANDLE_t mkl_handle_; /*!< MKL Direct Sparse Solver (DSS) handle */

  /*---> We've bumped into another computer scientists are muppets problem.
    Intel in their infinite intelligence and wisdom defined DSS
    interfaces that require certain integers.  They go to the trouble of
    defining certain values of these integers for us.  However, they
    don't do the definition in a manner consistent with the DSS function
    intefaces.  For example dss_create(handle, MKL_DSS_DEFAULTS) won't
    compile, due to a reference lvalue error.  So now I've got to band
    aid the problem with the followign static const ints.  Such a simple
    thing and yet they screw it up.  They really are just muppets.  */

  const MKL_INT mkl_dss_defaults = MKL_DSS_DEFAULTS; /*!< Inteface
                                                       fix to DSS */
  const MKL_INT mkl_dss_non_symmetric = MKL_DSS_NON_SYMMETRIC; /*!< Interface
                                                                 fix to DSS */
  const MKL_INT mkl_dss_auto_order = MKL_DSS_AUTO_ORDER;
  const MKL_INT mkl_dss_indefinite = MKL_DSS_INDEFINITE;
  const MKL_INT mkl_dss_zero_based_indexing = MKL_DSS_ZERO_BASED_INDEXING;
#endif
//****************************************************************************80
//! \brief MKLDSS
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  MKLDSS() = delete;

   
//****************************************************************************80
//! \brief Create : Wrapper for MKL DSS creation.
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//****************************************************************************80
  void Create();
  
//****************************************************************************80
//! \brief  MKLDSSDeleleFactorization : Call the delete function to delete
//!         factorized matrix for MKL DSS Solver.
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//!
//****************************************************************************80
  void Delete();

};// End class MKLDSS



#endif
