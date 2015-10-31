/*
 * PetscSolver.h
 *
 *  Created on: Oct 25, 2015
 *      Author: rabbit
 */

#ifndef PETSCSOLVER_H_
#define PETSCSOLVER_H_
//****************************************************************************80
//! \class PetscSolver
//! \brief  PetscSolver : Provides solvers from Petsc Library.
//! \nick
//! \version $Rev: 5 $
//****************************************************************************80
#include "my_incl.h"
#include "SystemUtils/SystemModule.h"
#include "DataStructures//Array1D.h"
#include "DataStructures/List2D.h"
#include "SparseMatrix/CSRMatrix.h"
#include "LinearSolver/LinearSolver.h"
#include "ParallelComm/Communication.h"
#ifdef PETSC
#include "petsc.h"
#endif


class PetscSolver : public LinearSolver
{
public:
//****************************************************************************80
//! \brief  PetscSolver : Provides solvers from Petsc Library.
//! \nick
//! \version $Rev: 5 $
//! \param[in] A The matrix in CSR format(Parallel too)
//! \param[in] graph_local2global Local 2 global node mapping of the graph used
//! to make the matrix A
//****************************************************************************80
  PetscSolver(CSRMatrix<realT>& A, intT num_its);

//****************************************************************************80
//! \brief  ~PetscSolver : Class destructor for PetscSolver
//! \nick
//! \version $Rev: 5 $
//****************************************************************************80
  ~PetscSolver();
//****************************************************************************80
//!
//! \brief Solve : fGMRES routine for solving the linear system
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
bool Factorize();

//****************************************************************************80
//!
//! \brief Solve : fGMRES routine for solving the linear system
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
void Solve(const List2D<realT>& rhs, List2D<realT>& out );
void DestroyPetscObjects();
private:
  //--->
  CSRMatrix<realT>& A_;//!< Matrix representing preconditioner
#ifdef PETSC
  Mat Amat_; //!< Petsc Matrix that we use to interface;
  Vec xvec_; // !< Petsc x vector
  Vec bvec_; // !< Petsc RHS Vector
  KSP ksp_;
  PC pc_;
#endif

  int nrow_;
  intT num_its_ = 1;//!< Number of iterations
  List2D<intT> colidx_global_;
  //Array1D<intT> petsc_global_row_;
  List2D<intT> petsc_global_row_;
  Array1D<intT> row_start_;
  bool destroy_objects_ = true;
  PetscSolver() = delete;
  void InitializeGlobalArrays();
  void InitializePetscObjects();
};

#endif /* PETSCSOLVER_H_ */
