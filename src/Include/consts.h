// -*-c++-*-
#ifndef MY_KINDDEFS_H
#define MY_KINDDEFS_H
//****************************************************************************80
//! \file consts.h 
//! \brief consts.h : Holds Various numerical constants that can be useful
//! \nick
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $ 
//****************************************************************************80

//---> Basic includes... we'll have these all the time. */
#include <math.h>    /*!< Always include the math library. */  
#include <iostream> /*!< Always include the output library. */
#include <complex> /*!< Include the native C++ complex template we'll use 
		     this. */
#include <iomanip> /*!< Include the ability to manipulate io. */	   

//!
//--> Type definitions		    
#ifdef SINGLE_PRECISION
typedef float DOUBLE;
#else
typedef double DOUBLE;
#endif

typedef std::complex <DOUBLE> COMPLEX; /*!< Definition of a double complex */ 
//---> Compile time constants that are useful
 static const double pi = 
   double(3.141592653589793238462643383279502884197); //!< Value of PI */
static const double half = double (.5); //!< Value of 1/2 
static const double onethird = double(double(1.0)/double(3.0)); /*!< Value of 
								  1.0/3.0 */
static const double twothird = double(double(2.0)/double(3.0)); /*!< Value of 
								  2.0/3.0 */

#endif
