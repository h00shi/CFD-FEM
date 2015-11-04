/*
 * Residual.h
 *
 *  Created on: Nov 1, 2015
 *      Author: rabbit
 */

#ifndef RESIDUAL_H_
#define RESIDUAL_H_

//****************************************************************************80
//! \brief An abstract class to represent residual interfaces and data
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//****************************************************************************80
class Residual{
public:
//****************************************************************************80
//! \brief Constructor for abstract residual class
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//****************************************************************************80
  Residual(){}

//****************************************************************************80
//! \brief Destructor for abstract residual class
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//****************************************************************************80
  virtual ~Residual() = 0;

//****************************************************************************80
//! \brief Abstract interface for computing a residual
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[out] resid Residual argument
//****************************************************************************80
  virtual void ComputeResidual(List2D<realT>& resid) = 0;

//****************************************************************************80
//! \brief Abstract interface for computing a residual and jacobian matrix
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[out] resid Residual argument
//! \param[out] jac Jacobian matrix
//****************************************************************************80
  template<class MatrixType>
  virtual void ComputeResidualJacobian(List2D<realT>& resid,
                                       MatrixType& jac) = 0;
//****************************************************************************80
//! \brief Abstract interface for computing a residual and jacobian matrix
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[out] v Vector to multiply by jacobian matrix
//! \param[out] Jv Jacobian times vector result
//****************************************************************************80
  virtual void ComputeJacobianVectorProduct(List2D<realT>& v,
                                            List2D<realT>& Jv) = 0;

private:

};


#endif /* RESIDUAL_H_ */
