//-*-c++-*-
#ifndef SURREAL_H
#define SURREAL_H

#include "Surreal/SurrealBase.h"
#include "Surreal/SurrealAdd.h"
#include "Surreal/SurrealSubtract.h"
#include "Surreal/SurrealMultiply.h"
#include "Surreal/SurrealDivide.h"
#include "Surreal/SurrealPow.h"
#include "Surreal/SurrealLogic.h"
#include "Surreal/SurrealUnary.h"

//****************************************************************************80
//! \class Surreal
//! \brief Class defining surreal numbers
//! \details  This is the class defining a surreal number i.e. a number
//!           with value and derivative w.r.t. an arbitrary number of
//!           parameters.
//! \nick
//! \version $Rev$
//! \tparam realT The datatype of real numbers, i.e. float or double
//! \tparam N The number of derivatives you want to take
//****************************************************************************80
template< typename realT, int N = 1>
class Surreal : public SurrealBase< Surreal<realT,N>, N > {
public:
  typedef realT realT_;
//++++++++++++++++++++++++++++++++ PUBLIC STUFF ++++++++++++++++++++++++++++++++
//****************************************************************************80
//! \brief The default constructor, no arguments.
//! \details
//! \nick
//! \version $Rev$
//!
//****************************************************************************80
  Surreal() {
    
  }// End default Surreal constructor

//****************************************************************************80
//!
//! \brief The constructor of Surreal from a "real" number
//! \details
//! \nick
//! \version $Rev$
//!
//****************************************************************************80
  Surreal(const realT& a) {
    //---> Set value to input
    value_ = a;
    //---> Set value of derivative to zero
    for( unsigned i = 0; i < N; i++){
      deriv_[i] = 0.0;
    }
  }// End Surreal real value constructor

//****************************************************************************80
//!
//! \brief The "copy" constructor that copys from a type SurrealBase
//! \details Not exactly a copy constructor
//! \nick
//! \version $Rev$
//! \param[in] a The instance to copy from
//****************************************************************************80
  template< class ExprType>
  Surreal(const SurrealBase<ExprType, N>& a)
  {
    const ExprType& ObjectCopyFrom = a.CastToDerived();

    //---> Set value
    value_ = ObjectCopyFrom.Value();
    //---> Copy the derivative
    for(int i = 0; i < N; i++){
      deriv_[i] = ObjectCopyFrom.Deriv(i);
    }

  }// End Surreal Copy constructor

//****************************************************************************80
//!
//! \brief Value : Get the reference to the value...modifiable
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
  inline realT& Value() {
    //---> Return value to user
    return(value_);
  }// End Value

//****************************************************************************80
//!
//! \brief Value : Get the value const correct version to return to const
//! \details
//! \nick
//! \version $Rev$
//****************************************************************************80
  inline const realT& Value() const {
    //---> Return value to user
    return(value_);
  }// End Value

//****************************************************************************80
//!
//! \brief Deriv : Get the ith partial derivative, i starts from 0.
//!                Gives modifiable reference.
//! \details
//! \nick
//! \version $Rev$
//! \param[in] i
//****************************************************************************80
 inline realT& Deriv(const int& i) {
    //---> Return value to user
    return(deriv_[i]);
  }// End Deriv

//****************************************************************************80
//!
//! \brief Deriv : Get the ith partial derivative, i starts from 0.
//!                Const correct version to return to const.
//! \details
//! \nick
//! \version $Rev$
//! \param[in] i
//****************************************************************************80
 inline const realT& Deriv(const int& i) const {
    //---> Return value to user
    return(deriv_[i]);
  }// End Deriv

//****************************************************************************80
//!
//! \brief  Stream operator
//! \details
//! \nick
//! \version $Rev$
//! \param[in] os The ostream object to stream to
//! \param[in] a The surreal number to stream
//! \return os The ostream object we streamed to
//****************************************************************************80
  friend std::ostream& operator << (std::ostream& os,
                                    const Surreal<realT,N>& a)
  {
    //---> Write the value
    os << "(" << a.value_ << ", { ";
    for(int i = 0; i < N; i++) os << a.deriv_[i] << " ";
    os << "} )";
    return(os);
  }// End operator <<

//-------------------------- ASSIGNMENT OPERATORS ------------------------------
//****************************************************************************80
//!
//! \brief  Assignment operator from real number.
//! \details
//! \nick
//! \version $Rev$
//! \param[in] rhs - Real valued right hand side of equal sign
//! \return *this The pointer to current instance of a Surreal
//****************************************************************************80
  inline Surreal& operator = (const realT& v)
  {
    //---> Set derivative to zero because we assigned a real value
    for(unsigned i = 0; i < N; i++) { // Derivative initialization
      deriv_[i] = 0.0;
    }
    //---> Set the value to v
    value_ = v;
    //---> Return the pointer to this instance
    return(*this);
  } // End operator =

//****************************************************************************80
//!
//! \brief  Assignement operator from SurrealBase type
//! \details
//! \nick
//! \version $Rev$
//! \tparam ExprType - The type of expression that we want to assign from
//! \param[in]   rhs - Right hand side of equal sign of type SurrealBase
//! \return    *this - The pointer to current instance of a Surreal
//****************************************************************************80
  template<class ExprType>
  inline Surreal& operator = (const SurrealBase<ExprType, N>& rhs)
  {
    const ExprType& expression_tree = rhs.CastToDerived();

    for(unsigned i = 0; i < N; i++) {//Derivative evaluation
      deriv_[i] = expression_tree.Deriv(i);
    } // End derivative evaluation

    //---> Assign the value
    value_ = expression_tree.Value();

    //---> Return the pointer to this instance
    return(*this);
  } // End operator =

//------------------------ COMPOUND ASSIGNMENT OPERATORS -----------------------
//****************************************************************************80
//! \brief  Addition assignment operator for a SurrealBase argument
//! \details  For the addition assignment operator, the value of the SurrealBase
//!           argument is added to the value of the current instance of Surreal.
//!               this.value_ += rhs.Value;
//!           Likewise, the derivative of the SurrealBase argument is added to
//!           the derivatives of the current instance of Surreal.
//!               this.deriv_ += rhs.Deriv;
//! \nick
//! \tparam ExprType The type of expression we are going to apply
//! \param[in]   rhs - Right hand side
//! \return *this The pointer to current instance of a Surreal
//****************************************************************************80
  template<class ExprType>
  inline Surreal& operator += (const SurrealBase<ExprType, N>& rhs)
  {
    const ExprType& expression_tree = rhs.CastToDerived();

    //---> Add derivative of argument to current derivatives
    for(unsigned i = 0; i < N; i++) {//Derivative evaluation
      deriv_[i] += expression_tree.Deriv(i);
    } // End derivative evaluation

    //---> Add value of argument to current value
    value_ += expression_tree.Value();

    //---> Return the pointer to this instance
    return(*this);
  } // End operator +=

//****************************************************************************80
//! \brief  Addition assignment operator for a realT argument
//! \details  For the addition assignment operator, the value of the realT
//!           argument is added to the value of the current instance of Surreal.
//!               this.value_ += rhs;
//!           The derivative is unchanged.
//!               this.deriv_ = this.deriv_;
//! \nick
//! \tparam ExprType The type of expression we are going to apply
//! \param[in]   rhs - Right hand side
//! \return *this The pointer to current instance of a Surreal
//****************************************************************************80
  inline Surreal& operator += (const realT& rhs)
  {
    //---> Do nothing with current derivatives
    //...
    //---> Add value of argument to current value
    value_ += rhs;
    //---> Return the pointer to this instance
    return(*this);
  } // End operator +=

//****************************************************************************80
//! \brief  Subtraction assignment operator for a SurrealBase argument
//! \details  For the subtraction assignment operator, the value of the
//!           SurrealBase argument is subtracted from the value of the current
//!           instance of Surreal.
//!               this.value_ -= rhs.Value;
//!           Likewise, the derivative of the SurrealBase argument is
//!           subtracted from the derivatives of the current instance of
//!           Surreal.
//!               this.deriv_ -= rhs.Deriv;
//! \nick
//! \tparam ExprType The type of expression we are going to apply
//! \param[in]   rhs - Right hand side
//! \return *this The pointer to current instance of a Surreal
//****************************************************************************80
  template<class ExprType>
  inline Surreal& operator -= (const SurrealBase<ExprType, N>& rhs)
  {
    const ExprType& expression_tree = rhs.CastToDerived();

    //---> Subtract derivative of argument from current derivatives
    for(unsigned i = 0; i < N; i++) {//Derivative evaluation
      deriv_[i] -= expression_tree.Deriv(i);
    } // End derivative evaluation

    //---> Subtract value of argument from current value
    value_ -= expression_tree.Value();

    //---> Return the pointer to this instance
    return(*this);
  } // End operator -=

//****************************************************************************80
//! \brief  Subtraction assignment operator for a realT argument
//! \details  For the subtraction assignment operator, the value of the
//!           realT  argument is subtracted from the value of the current
//!           instance of Surreal.
//!               this.value_ -= rhs;
//!           The derivative is unchanged.
//!               this.deriv_ = this.deriv_;
//! \nick
//! \tparam ExprType The type of expression we are going to apply
//! \param[in]   rhs - Right hand side
//! \return *this The pointer to current instance of a Surreal
//****************************************************************************80
  inline Surreal& operator -= (const realT& rhs)
  {
    //---> Do nothing with current derivatives
    //...
    //---> Subtract value of argument from current value
    value_ -= rhs;

    //---> Return the pointer to this instance
    return(*this);
  } // End operator -=

//****************************************************************************80
//! \brief  Multiplication assignment operator for a SurrealBase argument
//! \details  For the multiplication assignment operator, the value of the
//!           current instance of Surreal is multiplied by the value of the
//!           SurrealBase argument.
//!               this.value_ *= rhs.Value;
//!           The new derivative of the current instance is calculated with the
//!           product rule.
//!               this.deriv_ = this.value_ *rhs.Deriv +
//!                             this.deriv_ *rhs.Value
//! \nick
//! \tparam ExprType The type of expression we are going to apply
//! \param[in]   rhs - Right hand side
//! \return *this The pointer to current instance of a Surreal
//****************************************************************************80
  template<class ExprType>
  inline Surreal& operator *= (const SurrealBase<ExprType, N>& rhs)
  {
    const ExprType& expression_tree = rhs.CastToDerived();

    //---> Use product rule to compute new derivatives
    for(unsigned i = 0; i < N; i++) {//Derivative evaluation
      deriv_[i] = value_ * expression_tree.Deriv(i) +
        deriv_[i]*expression_tree.Value();
    } // End derivative evaluation

    //---> Multiply value of argument with current value
    value_ *= expression_tree.Value();

    //---> Return the pointer to this instance
    return(*this);
  } // End operator *=

//****************************************************************************80
//! \brief  Multiplication assignment operator for a realT argument
//! \details  For the multiplication assignment operator, the value of the
//!           current instance of surreal is multiplied by the realT argument.
//!               this.value_ *= rhs;
//!           Likewise, the derivative of the current instance of surreal is
//!           multiplied by the realT argument
//!               this.deriv_ *= rhs;
//! \nick
//! \tparam ExprType The type of expression we are going to apply
//! \param[in]   rhs - Right hand side
//! \return *this The pointer to current instance of a Surreal
//****************************************************************************80
  inline Surreal& operator *= (const realT& rhs)
  {
    //---> Multiply current derivative by value of input
    for(unsigned i = 0; i < N; i++) {//Derivative evaluation
      deriv_[i] *= rhs;
    } // End derivative evaluation

    //---> Multiply value of argument with current value
    value_ *= rhs;

    //---> Return the pointer to this instance
    return(*this);
  } // End operator *=

//****************************************************************************80
//! \brief  Division assignment operator for a SurrealBase argument
//! \details  For the division assignment operator, the value of the
//!           current instance of Surreal is divided by the value of the
//!           SurrealBase argument.
//!               this.value_ /= rhs.Value;
//!           The new derivative of the current instance is calculated with the
//!           quotient rule.
//!               this.deriv_ = (this.deriv_ *rhs.Value -
//!                              this.value_ *rhs.Deriv_) / (rhs.Value)^2
//! \nick
//! \tparam ExprType The type of expression we are going to apply
//! \param[in]   rhs - Right hand side
//! \return *this The pointer to current instance of a Surreal
//****************************************************************************80
  template<class ExprType>
  inline Surreal& operator /= (const SurrealBase<ExprType, N>& rhs)
  {
    const ExprType& expression_tree = rhs.CastToDerived();

    //---> Use quotient rule to compute new derivatives
    for(unsigned i = 0; i < N; i++) {//Derivative evaluation
      deriv_[i] = (deriv_[i] * expression_tree.Value() -
                   value_    * expression_tree.Deriv(i)) /
        (expression_tree.Value()*expression_tree.Value());
    } // End derivative evaluation

    //---> Divide current value by value of argument
    value_ /= expression_tree.Value();

    //---> Return the pointer to this instance
    return(*this);
  } // End operator /=

//****************************************************************************80
//! \brief  Division assignment operator for a realT argument
//! \details  For the division assignment operator, the value of the
//!           current instance of Surreal is divided by the value of the
//!           realT argument.
//!               this.value_ /= rhs;
//!           Likewise, the derivative of the current instance of surreal is
//!           divided  by the realT argument
//!               this.deriv_ /= rhs;
//! \nick
//! \tparam ExprType The type of expression we are going to apply
//! \param[in]   rhs - Right hand side
//! \return *this The pointer to current instance of a Surreal
//****************************************************************************80
  inline Surreal& operator /= (const realT& rhs)
  {
    //---> Use quotient rule to compute new derivatives
    for(unsigned i = 0; i < N; i++) {//Derivative evaluation
      deriv_[i] /= rhs;
    } // End derivative evaluation

    //---> Divide current value by value of argument
    value_ /= rhs;

    //---> Return the pointer to this instance
    return(*this);
  } // End operator /=

private:
//++++++++++++++++++++++++++++++++ PRIVATE STUFF +++++++++++++++++++++++++++++++
  realT value_;    /*!< The value of the number. */
  realT deriv_[N]; /*!< The derivatives of the number. */

}; //End class Surreal

#endif
