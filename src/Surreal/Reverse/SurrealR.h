//-*-c++-*-
#ifndef SURREALR_H
#define SURREALR_H

#include "Surreal/Reverse/SurrealRBase.h"
#include "Surreal/Reverse/SurrealRAdd.h"
#include "Surreal/Reverse/SurrealRSubtract.h"
#include "Surreal/Reverse/SurrealRMultiply.h"
#include "Surreal/Reverse/SurrealRDivide.h"
#include "Surreal/Reverse/SurrealRPow.h"
#include "Surreal/Reverse/SurrealRLogic.h"
#include "Surreal/Reverse/SurrealRUnary.h"
#include <stdexcept>
#include <iostream>

//****************************************************************************80
//! \class SurrealR
//! \brief Class defining surreal numbers
//! \details  This is the class defining a surreal number i.e. a number
//!           with value and derivative w.r.t. an arbitrary number of
//!           parameters.  This class calculates the derivatives in the reverse
//!           mode.
//! \nick
//! \tparam realT The datatype of real numbers, i.e. float or double
//! \tparam N The number of derivatives you want to take
//****************************************************************************80
template< typename realT, int N = 1>
class SurrealR : public SurrealRBase< SurrealR<realT,N>, realT, N > {
public:
//++++++++++++++++++++++++++++++++ PUBLIC STUFF ++++++++++++++++++++++++++++++++
//****************************************************************************80
//! \brief The default constructor, no arguments.
//! \nick
//****************************************************************************80
  SurrealR() : count_(0) {
    std::fill_n(deriv_,N,0.0);
    std::fill_n(deriv_ix_,N,-1);
    std::fill_n(ix_,N,-1);
  }// End default SurrealR constructor

//****************************************************************************80
//! \brief The constructor of SurrealR from a "real" number
//! \nick
//****************************************************************************80
  SurrealR(const realT a) : count_(0) {
    this->value_ = a;
    std::fill_n(deriv_,N,0.0);
    std::fill_n(deriv_ix_,N,-1);
    std::fill_n(ix_,N,-1);
  }// End Surreal real value constructor

//****************************************************************************80
//! \brief The default copy constructor to copy from other Surreals
//! \nick
//! \param[in] rhs The instance to copy from
//****************************************************************************80
  SurrealR(SurrealR<realT,N> const & rhs) = default;

//****************************************************************************80
//! \brief The "copy" constructor that copys from a type SurrealBase lvalue
//! \nick
//! \param[in] a The instance to copy from
//****************************************************************************80
  template< class ExprType>
  SurrealR(SurrealRBase<ExprType, realT, N> const & a)
  {
    ExprType const & expression_tree = a.CastToDerived();

    //---> Assign the value
    this->value_ = expression_tree.GetValue();

    //--->reset state of derivative arrays to default state
    std::fill_n(deriv_,N,0.0);
    std::fill_n(deriv_ix_,N,-1);
    std::fill_n(ix_,N,-1);
    count_ = 0;

    //seed root of tree has adjoint of 1.0
    expression_tree.Diff(1.0, deriv_, deriv_ix_, ix_, count_);
  }// End Surreal Copy constructor

//****************************************************************************80
//!
//! \brief Deriv : Get the ith partial derivative, i starts from 0.
//!                Gives modifiable reference.
//! \details  This function modifies the 'sparsity' pattern of the underlying
//!          gradient arrays by adding an additional non-zero.
//! \nick
//! \version $Rev$
//! \param[in] i
//****************************************************************************80
 inline realT& Deriv(const int i) {
    int target_ix  = ix_[i];
    if(target_ix == -1){
      deriv_ix_[count_] = i;
      ix_[i]   = count_;
      target_ix = count_++;
    }
    return deriv_[target_ix];
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
 inline realT Deriv(const int i) const {
   return GetDeriv(i);
  }// End Deriv

//****************************************************************************80
//! \brief SetDeriv : Set the ith partial derivative to value
//! \details This function modifies the 'sparsity' pattern of the underlying
//!          gradient arrays by adding an additional non-zero.
//! \nick
//! \param[in] i   - index of derivative to set (ith partial)
//! \param[in] val - value to set to
//****************************************************************************80
  void SetDeriv(const int i, const realT val) {
    int target_ix  = ix_[i];
    if(target_ix == -1){
      deriv_ix_[count_] = i;
      ix_[i]   = count_;
      target_ix = count_++;
    }
    deriv_[target_ix] = val;
  }// End SetDeriv

//****************************************************************************80
//! \brief GetDeriv : Get the ith partial derivative
//! \nick
//! \param[in] i - index of derivative to get (ith partial)
//****************************************************************************80
 inline const realT GetDeriv(const int i) const {
    int target_ix  = ix_[i];
    if(target_ix == -1){ //not found - return 0.0
      return 0.0;
    }else{
      return deriv_[target_ix];
    }
 }// End GetDeriv

//****************************************************************************80
//! \brief  Stream operator
//! \nick
//! \param[in] os The ostream object to stream to
//! \param[in] a The surreal number to stream
//! \return os The ostream object we streamed to
//****************************************************************************80
  friend std::ostream& operator << (std::ostream& os,
                                    const SurrealR<realT,N>& a)
  {
    //---> Write the value
    os << "(" << a.value_ << ", { ";
    for(unsigned int i = 0; i < N; i++) os << a.GetDeriv(i) << " ";
    os << "}";

    return(os);
  }// End operator <<

//-------------------------- ASSIGNMENT OPERATORS ------------------------------
//****************************************************************************80
//! \brief  Assignment operator from real number.
//! \nick
//! \param[in] value - Real valued right hand side of equal sign
//! \return *this The pointer to current instance of a Surreal
//****************************************************************************80
  inline SurrealR& operator = (const realT value)
  {
    //--->reset state of derivative arrays to default state
    std::fill_n(deriv_,N,0.0);
    std::fill_n(deriv_ix_,N,-1);
    std::fill_n(ix_,N,-1);
    count_ = 0;

    //---> Assign the value
    this->value_ = value;

    //---> Return the pointer to this instance
    return(*this);
  } // End operator =

//****************************************************************************80
//!\brief  Assignement operator from lvalue reference to a Surreal
//! \nick
//! \tparam ExprType - The type of expression that we want to assign from
//! \param[in]   rhs - Right hand side of equal sign of type SurrealBase
//! \return    *this - The pointer to current instance of a Surreal
//****************************************************************************80
  inline SurrealR& operator = (SurrealR<realT,N> const & rhs) = default;

//****************************************************************************80
//!\brief  Assignment operator from lvalue reference to a SurrealBase type
//! \nick
//! \tparam ExprType - The type of expression that we want to assign from
//! \param[in]   rhs - Right hand side of equal sign of type SurrealBase
//! \return    *this - The pointer to current instance of a Surreal
//****************************************************************************80
  template<class ExprType>
  inline SurrealR& operator = (SurrealRBase<ExprType, realT, N> const & rhs)
  {
    ExprType const & expression_tree = rhs.CastToDerived();

    //---> Assign the value
    this->value_ = expression_tree.GetValue();

    //--->reset state of derivative arrays to default state
    std::fill_n(deriv_,N,0.0);
    std::fill_n(deriv_ix_,N,-1);
    std::fill_n(ix_,N,-1);
    count_ = 0;

    //seed root of tree has adjoint of 1.0
    expression_tree.Diff(1.0, deriv_, deriv_ix_, ix_, count_);

    //---> Return the pointer to this instance
    return(*this);
  } // End operator =

//------------------------ COMPOUND ASSIGNMENT OPERATORS -----------------------
//****************************************************************************80
//! \brief  Addition assignment operator for a lvalue SurrealBase argument
//! \nick
//! \tparam ExprType The type of expression we are going to apply
//! \param[in]   rhs - Right hand side
//! \return *this The pointer to current instance of a Surreal
//****************************************************************************80
  template<class ExprType>
  inline SurrealR& operator += (SurrealRBase<ExprType, realT, N> const & rhs)
  {
    ExprType const & expression_tree = rhs.CastToDerived();
    //---> Calculate contribution of self contribution to derivative
    //--> lhs_adjoint is 1.0 - self contribution is identity so do nothing

    //---> Calculate contribution from right hand side
    expression_tree.Diff(1.0, deriv_, deriv_ix_, ix_, count_);

    //---> Add value of argument to current value
    this->value_ += expression_tree.GetValue();

    //---> Return the pointer to this instance
    return(*this);
  } // End operator +=

//****************************************************************************80
//! \brief  Addition assignment operator for a realT argument
//! \nick
//! \param[in]   rhs_value - Right hand side value
//! \return *this The pointer to current instance of a Surreal
//****************************************************************************80
  inline SurrealR& operator += (const realT rhs_value)
  {
    //---> Calculate contribution of self contribution to derivative
    //--> lhs_adjoint is 1.0 - self contribution is identity so do nothing

    //---> Multiply value of argument with current value
    this->value_ += rhs_value;

    //---> No contribution to derivative from right hand side (constant node)
    //...

    //---> Return the pointer to this instance
    return(*this);
  } // End operator +=

//****************************************************************************80
//! \brief  Subtraction assignment operator for a lvalue SurrealBase argument
//! \nick
//! \tparam ExprType The type of expression we are going to apply
//! \param[in]   rhs - Right hand side
//! \return *this The pointer to current instance of a Surreal
//****************************************************************************80
  template<class ExprType>
  inline SurrealR& operator -= (SurrealRBase<ExprType, realT, N> const & rhs)
  {
    ExprType const & expression_tree = rhs.CastToDerived();
    //---> Calculate contribution of self contribution to derivative
    //--> lhs_adjoint is 1.0 - self contribution is identity so do nothing

    //---> Calculate contribution from right hand side
    expression_tree.Diff(-1.0, deriv_, deriv_ix_, ix_, count_);

    //---> Add value of argument to current value
    this->value_ -= expression_tree.GetValue();

    //---> Return the pointer to this instance
    return(*this);
  } // End operator -=

//****************************************************************************80
//! \brief  Subtraction assignment operator for a realT argument
//! \nick
//! \param[in]   rhs_value - Right hand side value
//! \return *this The pointer to current instance of a Surreal
//****************************************************************************80
  inline SurrealR& operator -= (const realT rhs_value)
  {
    //---> Calculate contribution of self contribution to derivative
    //--> lhs_adjoint is 1.0 - self contribution is identity so do nothing

    //---> Multiply value of argument with current value
    this->value_ -= rhs_value;

    //---> No contribution to derivative from right hand side (constant node)
    //...

    //---> Return the pointer to this instance
    return(*this);
  } // End operator +=


//****************************************************************************80
//! \brief  Multiplication assignment operator for a lvalue SurrealBase argument
//! \details
//! \nick
//! \tparam ExprType The type of expression we are going to apply
//! \param[in]   rhs - Right hand side
//! \return *this The pointer to current instance of a Surreal
//****************************************************************************80
  template<class ExprType>
  inline SurrealR& operator *= (SurrealRBase<ExprType, realT, N> const & rhs)
  {
    ExprType const & expression_tree = rhs.CastToDerived();

    //---> Calculate the rhs value
    realT rhs_value = expression_tree.GetValue();

    //---> Calculate contribution of self contribution to derivative
    realT lhs_adjoint = rhs_value;
    for(unsigned i = 0; i < count_; i++) {//Derivative evaluation
      deriv_[i] *= lhs_adjoint;
    } // End derivative evaluation

    //---> Calculate contribution from right hand side
    expression_tree.Diff(this->value_, deriv_, deriv_ix_, ix_, count_);

    //---> Multiply value of argument with current value
    this->value_ *= rhs_value;

    //---> Return the pointer to this instance
    return(*this);
  } // End operator *=

//****************************************************************************80
//! \brief  Multiplication assignment operator for a realT argument
//! \nick
//! \param[in]   rhs_value - Right hand side value
//! \return *this The pointer to current instance of a Surreal
//****************************************************************************80
  inline SurrealR& operator *= (const realT rhs_value)
  {
    //---> Calculate contribution of self contribution to derivative
    realT lhs_adjoint = rhs_value;
    for(unsigned i = 0; i < count_; i++) {//Derivative evaluation
      deriv_[i] *= lhs_adjoint;
    } // End derivative evaluation

    //---> No contribution from right hand side (constant node)
    //...

    //---> Multiply value of argument with current value
    this->value_ *= rhs_value;

    //---> Return the pointer to this instance
    return(*this);
  } // End operator *=


//****************************************************************************80
//! \brief  Division assignment operator for a lvalue SurrealBase argument
//! \nick
//! \tparam ExprType The type of expression we are going to apply
//! \param[in]   rhs - Right hand side
//! \return *this The pointer to current instance of a Surreal
//****************************************************************************80
  template<class ExprType>
  inline SurrealR& operator /= (SurrealRBase<ExprType, realT, N> const & rhs)
  {
    ExprType const & expression_tree = rhs.CastToDerived();

    //---> Calculate the rhs value
    realT rhs_value = expression_tree.GetValue();

    //---> Calculate contribution of self contribution to derivative
    realT lhs_adjoint = 1.0 / rhs_value;
    for(unsigned i = 0; i < count_; i++) {//Derivative evaluation
      deriv_[i] *= lhs_adjoint;
    } // End derivative evaluation

    //---> Calculate contribution from right hand side
    expression_tree.Diff(-this->value_ / (rhs_value * rhs_value),
                         deriv_, deriv_ix_, ix_, count_);

    //---> Multiply value of argument with current value
    this->value_ /= rhs_value;

    //---> Return the pointer to this instance
    return(*this);
  } // End operator /=


//****************************************************************************80
//! \brief  Division assignment operator for a realT argument
//! \nick
//! \param[in]   rhs_value - Right hand side value
//! \return *this The pointer to current instance of a Surreal
//****************************************************************************80
  inline SurrealR& operator /= (const realT rhs_value)
  {
    //---> Calculate contribution of self contribution to derivative
    realT lhs_adjoint = 1.0 / rhs_value;
    for(unsigned i = 0; i < count_; i++) {//Derivative evaluation
      deriv_[i] *= lhs_adjoint;
    } // End derivative evaluation

    //---> No contribution from right hand side (constant node)
    //...

    //---> Multiply value of argument with current value
    this->value_ /= rhs_value;

    //---> Return the pointer to this instance
    return(*this);
  } // End operator /=

//****************************************************************************80
//! \brief Conversion operator, allows calls
//!                   realT x = static_cast<realT> Surreal y;
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[in] rhs Surreal that you want to convert to a realT
//****************************************************************************80
  inline operator realT() {
    return this->value_;
  } // End convert to realT

//****************************************************************************80
//! \brief Conversion operator, allows calls
//!                   realT x = static_cast<realT> Surreal y; (const correct)
//! \details
//! \nick
//! \version $Rev$
//! \date $Date$
//! \param[in] rhs Surreal that you want to convert to a realT
//****************************************************************************80
  inline operator realT() const {
    return this->value_;
  } // End convert to realT

//****************************************************************************80
//! \brief  Diff - differentiation function
//! \details USER - DO NOT CALL
//! \nick
//! \param[in]
//****************************************************************************80
  void Diff(realT adjoint, realT (&deriv)[N], int (&deriv_ix)[N], int (&ix)[N],
            unsigned int & count) const {
#ifdef DEV_DEBUG
    //check if this leaf node is referring to itself - if it is, throw an error
    if(deriv_ == deriv){throw std::invalid_argument("A surreal object cannot "
                                                    "be on both the left and "
                                                    "right hand sides of the "
                                                    "equation.");}
#endif

    //LEAF NODE CONTRIBUTION
    for(unsigned i = 0; i < count_; i++) {
      int primary_ix = deriv_ix_[i];
      int target_ix  = ix[primary_ix];
      if(target_ix == -1){
        deriv_ix[count] = primary_ix;
        ix[primary_ix]  = count;
        target_ix = count++;
      }
      deriv[target_ix] += adjoint*deriv_[i];
    } // End derivative evaluation
  }// End Deriv

private:
//++++++++++++++++++++++++++++++++ PRIVATE STUFF +++++++++++++++++++++++++++++++
  //--->'sparse' storage of derivatives
  realT deriv_[N];  /*!< The value of derivatives wrt deriv_ix */
  int deriv_ix_[N]; /*!< index of var that derivatives are taken wrt */
  int ix_[N];       /*!< mapping from deriv_ix to position in deriv */
  unsigned int count_ = 0;
}; //End class SurrealR

#endif
