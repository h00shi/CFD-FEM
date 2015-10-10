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

#ifndef INC_STRING
#define INC_STRING
#include <string> 
#endif

#ifndef INC_IOSTREAM
#define INC_IOSTREAM
#include <iostream>
#endif


//----> Finally define precision of integer and real 
#include "precision.h"
#endif //END MY_INCL_H
