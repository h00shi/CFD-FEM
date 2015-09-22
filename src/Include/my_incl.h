// -*-c++-*-
#ifndef MY_INCL_H
#define MY_INCL_H
//****************************************************************************80
//! \file my_incl.h 
//! \brief my_incl.h : This is the basic include that should be used for every
//!                    source file in the code.  It ensures screen io, math 
//!                    and complex numbers are always availible.  Why C++ 
//!                    doesn't do this is beyond me.   
//!
//! \nick
//! \version $Rev: 5 $
//! \date $Date: 2013-10-21 14:35:02 -0700 (Mon, 21 Oct 2013) $
//****************************************************************************80

//---> Basic includes... we'll have these all the time. */
#ifndef INC_MATH
#define INC_MATH
#include <math.h>    /*!< Always include the math library. */ 
#endif

#ifndef INC_COMPLEX
#define INC_COMPLEX
#include <complex> /*!< Include the native C++ complex template we'll use 
		     this. */
#endif

#ifndef INC_STDIO
#define INC_STDIO
#include <stdio.h> /*!< Standard IO header file...may need for IO. */
#endif

#ifndef INC_IOSTREAM
#define INC_IOSTREAM
#include <iostream> /*!< IOSTREAM header file gives access to ofstream object. */
#endif

#ifndef INC_IOMANIP
#define INC_IOMANIP
#include <iomanip> /*!< Header file for manipulating IO. */
#endif

#ifndef INC_SSTREAM
#define INC_SSTREAM
#include <sstream> /*!< Yet another IO related header file...why so many C++? */
#endif // END INC_SSTREAM

//---> 
#ifndef INC_STRING
#define INC_STRING
#include <string> 
#endif

//---> Namespaces
#ifndef INC_STD
#define INC_STD
//#include <stdlib> /*!< Header file for stdlib */
#endif // END INC_STD

//---> MPI
#ifndef INC_MPI
#define INC_MPI

#ifdef HAVE_MPI
#include <mpi.h> /*<! Include MPI header file if we are compiling with it. */
#endif

#endif // END INC_MPI

//----> Finally define precision of integer and real 
typedef int intT;
typedef double realT;

#include "consts.h"
#include "MathUtilities/complexify.h"
#endif //END MY_INCL_H
