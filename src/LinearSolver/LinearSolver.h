/*
 * LinearSolver.h
 *
 *  Created on: Oct 25, 2015
 *      Author: rabbit
 */

#ifndef LINEARSOLVER_H_
#define LINEARSOLVER_H_
#include "my_incl.h"
#include "DataStructures/List2D.h"
//****************************************************************************80
//! \brief LinearSolver : A class that represents the interface to a linear
//!                      solver.  [A]{x}={b}.
//! \details Please note that no matrix is required, to solve one need only
//!          specify the child of this class per spec and provide b, x.
//! \nick
//! \version $Rev$
//****************************************************************************80
class LinearSolver
{
public:
//****************************************************************************80
//!
//! \brief LinearSolver : LinearSolver class default tolerances
//! \nick
//! \version $Rev$
//****************************************************************************80
  LinearSolver();
//****************************************************************************80
//!
//! \brief LinearSolver : LinearSolver class destructor
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
  virtual ~LinearSolver() = 0;
//****************************************************************************80
//!
//! \brief LinearSolver : LinearSolver class constructor with variable
//!                       tolerances
//! \nick
//! \version $Rev$
//****************************************************************************80
  LinearSolver(const realT atol, const realT rtol);
//****************************************************************************80
//!
//! \brief Factorize: This function is the stage where all setup related to
//!                    to setting up the solver is done.
//! \nick
//! \version $Rev$
//****************************************************************************80
  virtual bool Factorize() = 0;
//****************************************************************************80
//!
//! \brief Solve: Given a rhs-vector and output vector x solve the system
//!        Return some resemblance of {x}=[A]^{-1}{b}.
//! \nick
//! \version $Rev$
//****************************************************************************80
  virtual void Solve(const List2D<realT>& b, List2D<realT>& x) = 0;

private:
  realT atol_ = 5.0e-15;//!< Absolute tolerance for linear system solution
  realT rtol_ = 5.0e-14;//!< Relative tolerance for linear system solution;
};



#endif /* LINEARSOLVER_H_ */
