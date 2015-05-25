// -*-c++-*-
#ifndef BLOCK_MATRIX_H
#define BLOCK_MATRIX_H
//****************************************************************************80
//! 
//! \brief This is the header file of the matrix base class SquareMatrix.    
//! \details The idea here is that functions and data members that are 
//!  common to square matrices live in the base class. The class extends an
//! Array2D class in order to provide square matrix functionality rather than
//! accepting an Array2D and performing operations on it. 
//! \ryan
//!  $Rev: 6 $
//!  $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//****************************************************************************80
#include "Array1D.h"
#include "Array2D.h"
#include "my_incl.h"
#include "consts.h"

template < class realT>
class SquareMatrix : public Array2D<realT>  {

protected:
//+++++++++++++++++++++++++++++++ PROTECTED STUFF ++++++++++++++++++++++++++++++
  /*---> The following variables are intended to be accessed from the derived 
    classes. */
  
  intT nSize; /*!< # of rows/columns for input matrix*/
  

public:
//****************************************************************************80
//!
//! \brief SquareMatrix: The constructor accepts a seedSize and guarentees that
//! the underlying Array2D is square
//! \details
//! \ryan
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] seedSize the input size
//****************************************************************************80
  SquareMatrix(const intT seedSize) : Array2D<realT>(seedSize,seedSize) { 
    nSize = seedSize;

  }// End Square Matrix Constructor

//****************************************************************************80
//!
//! \brief ~SquareMatrix : The destructor but we declare it anyway
//! \details
//! \ryan
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! 
//****************************************************************************80
  ~SquareMatrix() {

  } // End ~SquareMatrix
  
//****************************************************************************80
//!
//! \brief invert : inverts this square matrix that contains any realT 
//!    including Surreal variables
//! \details : Pivoted LU decomposition along with backward/forward substitution
//! is employed to invert a small square matrix. The class is templated and can
//! be used with double/float variables or Surreal variables. If Surreal 
//! variables are employed, derivatives of the square matrix inverse operation
//! wrt the surreal variables are computed with operator overloading. 
//! \ryan
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//****************************************************************************80
  void invert()
  {
    // Temporary storage for matrix inverse
    Array2D<realT> matInInv(nSize,nSize);  
    // integer pivot Array1D, swaps rows in the matrix for optimal LU decomp
    // based on the abs of the largest element value for a given column
    Array1D<intT> pivot(nSize);
    // the inverse of matrix A is actually computed by solving Ax = b for x
    // the ith element of b is set to 1.0 and the rest of b is set to 0.0
    // which allows for the solution x of Ax = b to corrrespond to the ith
    // column of inv(A) 
    Array1D<realT> b(nSize);
    // each column of inv(A) is stored in x
    Array1D<realT> x(nSize);

    // function to perform pivoted LU decomposition on this square matrix.
    // This square matrix is overwritten with the LU decomposition operation
    // the integer pivot array is passed by reference and modified within
    // the function.
   
    LUdecomp(pivot);
    //lu(pivot);
    for (intT i = 0; i < nSize; i++) {
      // set up the b array 
      for (intT j = 0; j < nSize; j++) b(j) = 0.0;
      b(i) = 1.0;
      // this is the backward/forward substutute funciton 
      // this function calculates the ith column of inv(A) and accepts 
      // const pivot and b arrays
      _backForwardSub(pivot,b,x);
      for (intT j = 0; j < nSize; j++) matInInv(j,i) = x(j);
    } // End For

    // this square matrix is overwritten with the inv(A)
    for (intT i = 0; i < nSize; i++) {
      for (intT j = 0; j < nSize; j++) {
        Array2D<realT>::operator()(i,j) = matInInv(i,j);
      } 
    }
      
  }// End invert

//****************************************************************************80
//!
//! \brief squareSolve:solve Ax=b where A is this square matrix that contains 
//!    any realT including Surreal variables
//! \details : Pivoted LU decomposition along with backward/forward substitution
//! is employed to solve a small square matrix. The class is templated and can
//! be used with double/float variables or Surreal variables. If Surreal 
//! variables are employed, derivatives of the square matrix inverse operation
//! wrt the surreal variables are computed with operator overloading. 
//! \ryan
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] b The rhs
//! \param[out] x The solution
//****************************************************************************80
  void squareSolve(const Array1D<realT>& b, Array1D<realT> &x)
  {
    // integer pivot Array1D, swaps rows in the matrix for optimal LU decomp
    // based on the abs of the largest element value for a given column
    Array1D<intT> pivot(nSize);

    // function to perform pivoted LU decomposition on this square matrix.
    // This square matrix is overwritten with the LU decomposition operation
    // the integer pivot array is passed by reference and modified within
    // the function.
    LUdecomp(pivot);
    // this is the backward/forward substutute funciton 
    // this function calculates b*inv(A) and accepts 
    // const pivot and b arrays
    _backForwardSub(pivot,b,x);
 
  }// End squareSolve

//****************************************************************************80
//!
//! \brief LUdecomp:calculate the LU decomposition of this square matrix  
//! \details
//! \ryan
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] pivot array
//****************************************************************************80
  void LUdecomp(Array1D<intT>& pivot) { 
    
    /*---> Maxpivotrow is the row of a column of this square matrix that 
      contains the largest magnitude value as compared to all of the values 
      in a column of this square matrix */
    intT maxpivotrow;
    /*---> Initial first element of a column of this square matrix before 
      pivoting procedure is commenced. */
    intT initialfirstrow;
    /*---> maxpivotvalue is the largest magnitude value as compared to all the 
      values in a column of this square matrix. */
    realT maxpivotvalue;
    Array1D<realT> temp(nSize);
    
    // the pivot array is initialized to the original orientation of the 
    // square matrix
    for (intT i = 0; i < nSize; i++) pivot(i) = i;
    
    for (intT j = 0; j < nSize-1; j++) {   /* loop through all the columns
					      of this square matrix except 
					      the last one */
      
      maxpivotvalue = 0.0;  //---> Maxpivot value is initialized to 0.0
      maxpivotrow = j; //---> Initialize to the current column.
      
      /*---> For a given column of this square matrix find the largest magnitude
	value. */ 
      /*---> The (i,j) element of this square matrix is accessed using the base 
	class Array2D<realT>::operator()(i,j) */
      for (intT i = j; i < nSize; i++) {  //pivot loop
	if (abs(Array2D<realT>::operator()(pivot(i),j))>
            abs(maxpivotvalue)) {
          maxpivotvalue = abs(Array2D<realT>::operator()(pivot(i),j));
          maxpivotrow = i;
        }
      }  // end pivot loop
      /*---> the initial first row and the maxpivot row values are swapped
	in the pivot array */
      initialfirstrow = pivot(j);
      pivot(j) = pivot(maxpivotrow);
      pivot(maxpivotrow) = initialfirstrow;
        
      for(intT i = j + 1; i < nSize; i++){
	Array2D<realT>::operator()(pivot(i),j) /=
	  Array2D<realT>::operator()(pivot(j),j);
	for(intT k = j + 1; k < nSize; k++){
	  Array2D<realT>::operator()(pivot(i),k) -= 
	    Array2D<realT>::operator()(pivot(i),j)*
	    Array2D<realT>::operator()(pivot(j),k);
	}
      }
    }  // end column loop
   
  }// End LUdecomp

//****************************************************************************80
//!
//! \brief _backForwardSub :return by reference each column vector of this
//!                         square matrix inverse
//! \details
//! \ryan
//! \version $Rev: 6 $
//! \date $Date: 2013-10-28 10:08:01 -0700 (Mon, 28 Oct 2013) $
//! \param[in] pivot Pivot array
//! \param[in] b Rhs array
//! \param[out] x Solution of Ax=b
//****************************************************************************80
  void _backForwardSub(const Array1D<intT>& pivot, const Array1D<realT>& b, 
                       Array1D<realT>& x) { 
    Array1D<realT> y(nSize);
    // perform backward substitution 
    for (intT i = 0; i < nSize; i++) y(i) = b(pivot(i));

    for (intT j = 0; j < nSize-1; j++) {
      for (intT i = j+1; i < nSize; i++) {
        y(i) = y(i) - Array2D<realT>::operator()(pivot(i),j)*y(j);
      }
    }
    // perform forward substitution
    for (intT i = 0; i < nSize; i++) {
      x(i) = y(i)/Array2D<realT>::operator()(pivot(i),i);
    }

    for (intT j = nSize-1; j > 0; j--) {
      for (intT i = j-1; i > -1; i--) {
        x(i) = x(i)- Array2D<realT>::operator()(pivot(i),j)*x(j)/
                Array2D<realT>::operator()(pivot(i),i);

      }
    }

  } // End _backForwardSub

//****************************************************************************80
//!  
//! \brief backup lu function
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! 
//****************************************************************************80
  void lu(Array1D<intT>& pivot)
  {
    
    for(intT c = 0; c < nSize; c++){pivot(c) = c;}

    for(intT c = 0; c < nSize - 1; c++) { // Column loop 
      //---> For the row we're on find the colume with the largest value
      intT max_loc = c;
      realT max_val = abs( Array2D<realT>::operator()(c,c) );
      
      for(intT r = c; r < nSize; r++) { // Pivot 
    	if( abs(Array2D<realT>::operator()(r, c)) > max_val ) {
    	  max_val = abs(Array2D<realT>::operator()(r, c));
    	  max_loc = r;
    	}
      } // End Pivot
      
      if( max_loc > c ) { // Move row 
    	intT itemp = pivot(c);
    	pivot(c) = pivot(max_loc);
    	pivot(max_loc) = itemp;
	
   
	//---> Do the pivoting
	for(intT j = 0; j < nSize; j++){ // Loop columns of max_loc
	  realT temp_element = Array2D<realT>::operator()(c, j);
	  Array2D<realT>::operator()(c, j) = 
	    Array2D<realT>::operator()(max_loc, j);
	  
	  Array2D<realT>::operator()(max_loc, j) = temp_element;
	} // End loop columns of max_loc
      } // End move row
      
      //---> Factorize
      for(intT r = c + 1; r < nSize; r++){ 
    	Array2D<realT>::operator()(r,c) /= Array2D<realT>::operator()(c,c);
    	for(intT k = c + 1; k < nSize; k++){
    	  Array2D<realT>::operator()(r,k) -= 
    	    Array2D<realT>::operator()(r,c)*Array2D<realT>::operator()(c,k); 
    	}
      }

    }// End outer loop 
    
  } // End lu

//****************************************************************************80
//!
//! \brief
//! \details
//! \nick 
//! \version $Rev$
//! \date $Date$
//! 
//****************************************************************************80
  void lu_solve(const Array1D<intT>& pivot, const Array1D<realT>& b, 
	   Array1D<realT>& x)
  {

    for(intT i = 0; i < nSize; i++){
      x(i) = b(pivot(i));
    }

    for(intT i = 1; i < nSize; i++){
      for(intT j = 0; j < i; j++){ 
	x(i) -= Array2D<realT>::operator()(i,j)*x(j);
      }
    }
     
    x(nSize-1) /= Array2D<realT>::operator()(nSize-1,nSize-1);
    
    for(intT i = nSize - 2; i > -1; i--){
      for(intT j = nSize - 1; j > i; j--){
	x(i) -= Array2D<realT>::operator()(i,j)*x(j);
      }
      x(i) /= Array2D<realT>::operator()(i,i);
    }
    

  } // End lu_solve

}; // End class SquareMatrix
#endif
