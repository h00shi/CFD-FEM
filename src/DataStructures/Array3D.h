//-*-c++-*-

#ifndef ARRAY3D_H
#define ARRAY3D_H
//****************************************************************************80
//! \class Array3D
//! \brief  Array3D :  A 3-D array class that we can use for small stuff only
//! \details The ideas in this class are due to
//! A. Katz but have been re-implemented by N. Burgess
//! \nick
//! \version $Rev: 5 $
//! \tparam dataT Template parameter for type of data
//****************************************************************************80
#include "my_incl.h"
#include "SystemUtils/SystemModule.h"
#include <sstream>
#include "DataStructures/Array1D.h"
template <typename dataT>
class Array3D : public Array1D<dataT>
{
  //+++++++++++++++++++++++++++++++ PUBLIC STUFF +++++++++++++++++++++++++++++++
public:
//****************************************************************************80
//!
//! \brief Array3D : Array class default constructor
//! \details
//! \nick
//! \version $Rev: 5 $
//!
//****************************************************************************80
  Array3D(){

    size1_ = 0;
    size2_ = 0;
    size3_ = 0;

  } // End Array3D

//****************************************************************************80
//!
//! \brief Array3D : Array class constructor with size specification
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] n1 The size of the first array dimension
//! \param[in] n2 The size of the second array dimension
//! \param[in] n3 The size of the thrid array dimension
//****************************************************************************80
  Array3D(intT n1, intT n2, intT n3) : Array3D()
  {
    this->initialize(n1,n2,n3);
  } // End Array3D

//****************************************************************************80
//!
//! \brief Array3D : Array3D class move constructor
//! \details
//! \nick
//! \param[in] from_array - rvalue reference that we are "moving" data from
//****************************************************************************80
  Array3D(Array3D<dataT>&& from_array) :
    Array1D<dataT>(std::move(from_array))
  {
    //pilfer other's resource
    size1_ = from_array.size1_;
    size2_ = from_array.size2_;
    size3_ = from_array.size3_;


    //reset from_array
    from_array.size1_ = 0;
    from_array.size2_ = 0;
    from_array.size3_ = 0;

  }

//****************************************************************************80
//!
//! \brief operator= : Array3D class move assignment operator
//! \details
//! \nick
//! \param[in] from_array - rvalue reference that we are "moving" data from
//! \return    this  - returns this object with data moved into it
//****************************************************************************80
  Array3D<dataT>& operator=(Array3D<dataT>&& from_array)
  {
    Array1D<dataT>::operator=(std::move(from_array));

    //--->pilfer from_array's resource
    size1_ = from_array.size1_;
    size2_ = from_array.size2_;
    size3_ = from_array.size3_;

    //--->reset from_array
    from_array.size1_ = 0;
    from_array.size2_ = 0;
    from_array.size3_ = 0;

    return *this;
  }



//****************************************************************************80
//!
//! \brief ~Array3D : Destructor for class Array3D
//! \details
//! \nick
//! \version $Rev: 5 $
//!
//****************************************************************************80
  ~Array3D(){
    //---> Reset the size variable
    size1_ = 0;
    size2_ = 0;
    size3_ = 0;
  } // End ~Array3D

//****************************************************************************80
//!
//! \brief initialize : A function to allocate the size of the array in the case
//!                    that the size is not known at instantiation.
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] n1 The size of the first array dimension
//! \param[in] n2 The size of the second array dimension
//! \param[in] n3 The size of the thrid array dimension
//****************************************************************************80
  void initialize(intT n1, intT n2, intT n3){
    //---> Set value
    size1_ = n1;
    size2_ = n2;
    size3_ = n3;
    Array1D<dataT>::initialize(n1*n2*n3);

  } // End initialize

//****************************************************************************80
//!
//! \brief () : Parenthesis operator for acessing array, for const. correctness
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] i First index of the array you want to acess with parenthesis
//! \param[in] j Second index of the array you want to access with parenthesis
//! \param[in] k Third index of the arrary you want to access with parenthesis
//****************************************************************************80
  const dataT& operator () (intT i, intT j, intT k) const
  {
#ifdef DEV_DEBUG
    this->CheckBounds(i,j,k);
#endif

    //---> The return of this operator is the ith reference of data pointer
    return Array1D<dataT>::operator()(size3_*(i*size2_ + j) + k);
  }

//****************************************************************************80
//!
//! \brief () : Parenthesis operator for acessing array
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] i First index of the array you want to acess with parenthesis
//! \param[in] j Second index of the array you want to access with parenthesis
//! \param[in] k Third index of the arrary you want to access with parenthesis
//****************************************************************************80
  dataT& operator () (intT i, intT j, intT k)
  {
#ifdef DEV_DEBUG
    this->CheckBounds(i,j,k);
#endif

    //---> The return of this operator is the ith reference of data pointer
    //calls constant version of operator (i,j,k)
    return const_cast<dataT&>(static_cast<const Array3D&>(*this).
                                 operator()(i,j,k));
  } // End

//****************************************************************************80
//!
//! \brief get_dims  : For the specified index gets the size of that dimension
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] dim The dimension you want the size of
//****************************************************************************80
  intT get_size(intT dim)
  {

    //---> Check the value of dim to return the correct size
    switch (dim) { // check_dim
    case 0:
      return(size1_);
      break;
    case 1:
      return(size2_);
      break;
    case 2:
      return(size3_);
      break;
    default:
      return(-99);
    } // End check_dim

  }// End get_size

//****************************************************************************80
//!
//! \brief operator << : A facility to print the Array3D class the screen
//! \details
//! \nick
//! \version $Rev$
//! \param[in] os The ostream operator we stream to
//! \param[in] a The array we wish to stream
//! \return os Returns the modified value of os
//****************************************************************************80
  friend std::ostream& operator << (std::ostream& os, const Array3D<dataT>& a)
  {
    //---> Write the data to ostream object
    for(intT i = 0; i < a.size1_; i++){
      os << i << ":" << std::endl;
      for(intT j = 0; j < a.size2_; j++) {
        os << "\t" << j << ": ";
        for(intT k = 0; k < a.size3_; k++){
          os << a(i,j,k) << " ";
        }
        os << std::endl;
      }
    }
    return(os);
  } // End operator <<

//****************************************************************************80
//!
//! \brief get_ptr : Allow acess to pointer at specified address
//! \details  Constant correct version
//! \jun
//! \version
//! \param[in] i The first address you want to access
//! \param[in] j The second address you want to access
//! \param[in] k The third address you want to access
//****************************************************************************80
  dataT const * get_ptr(const intT& i, const intT& j,
                       const intT& k) const
  {
#ifdef DEV_DEBUG
    this->CheckBounds(i,j,k);

#endif

    return Array1D<dataT>::get_ptr( size3_*(i*size2_ + j) + k );
  }// End get_ptr

//****************************************************************************80
//!
//! \brief get_ptr : Allow acess to pointer at specified address
//!
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] i The first address you want to access
//! \param[in] j The second address you want to access
//! \param[in] k The third address you want to access
//****************************************************************************80
  dataT* get_ptr(const intT& i, const intT& j,
                       const intT& k)
  {
 #ifdef DEV_DEBUG
    this->CheckBounds(i,j,k);

#endif

    return Array1D<dataT>::get_ptr(+ size3_*(i*size2_ + j) + k );
  }// End get_ptr

//****************************************************************************80
//! \brief MemoryDiagnostic : Prints the size and memory information to user.
//!        
//! \details This is a very useful feature that helps the user/developer know
//!          how big the variable is and how much memory it consumes 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] var_name A string containing the variable name
//****************************************************************************80
  std::string MemoryDiagnostic(std::string const & var_name) {
   
     std::ostringstream stream;
     //---> Use operators to write to ostringstream...this is nice an easy
     stream << var_name << "(" << this->get_size(0) << ", "
            << this->get_size(1) << ", " 
	    << this->get_size(2) <<  "):\t "
	    << this->get_mem() << " MB" << std::endl;
     
     //---> However ostringstream is not copyable...so return the string
     return stream.str();

  }

private:
  //+++++++++++++++++++++++++++++++ PRIVATE STUFF ++++++++++++++++++++++++++++++

  intT size1_; /*!< Size of dimension 1 */
  intT size2_; /*!< Size of dimension 2 */
  intT size3_; /*!< Size of dimension 3 */

//****************************************************************************80
//!
//! \brief Array3D : Copy constructor...private so it can't be called, which
//!                  blocks copy construction of this class.
//! \details A copy constructor is a constructor that takes a argument of the
//!      type defined by the class.  In this case Array3D is the class.
//!      It then seeks to define operators for constructing a new instance
//!      using the old one.  We are purposely blocking this capability by
//!      intersting this function as a private functions with NO CODE.
//!
//!      NOTE: Just because this function does nothing does not make it
//!      acceptable to omit the documentation and clear detailed commentary.
//!      Note the detail above and follow it as a rule.
//! \nick
//! \version $Rev: 5 $
//****************************************************************************80
  Array3D(const Array3D<dataT>&) = delete; //blocked

//****************************************************************************80
//!
//! \brief = : The assignment operator... private so it can't be called, which
//!            blocks assignment of one array to another.
//! \details This is the assignement operator.  We are declaring it private to
//!      block the capability Array3D x = Array3D y type of behavior.
//!      We are purposely blocking this capability by
//!      intersting this function as a private operator with NO CODE.
//!
//!      NOTE: Just because this function does nothing does not make it
//!      acceptable to omit the documentation and clear detailed commentary.
//!      Note the detail above and follow it as a rule.
//! \nick
//! \version $Rev: 5 $
//****************************************************************************80
  Array3D& operator = (const Array3D<dataT>&) = delete; //blocked

//****************************************************************************80
//!
//! \brief CheckBounds: Debug feature - checks the bounds of the array
//! \details
//! \nick
//! \version
//! \param[in] i - i index
//****************************************************************************80
  void CheckBounds(intT i, intT j, intT k) const {
    if( i >= size1_) {
      std::cerr << "ERROR: In Array3D.h - Over bounds on 1st index. "
		<< "Acessing Array3D::("<< i <<","<< j <<","<< k 
		<< "). Size of Array3D "
		<< "is " << size1_ << ","<< size2_ << "," << size3_ << std::endl;
      SystemModule::my_exit();
    }
    if( j >= size2_) {
      std::cerr << "ERROR: In Array3D.h - Over bounds on 2nd index. "
		<< "Acessing Array3D::("<< i <<","<< j <<","<< k 
		<< "). Size of Array3D "
		<< "is " << size1_ << ","<< size2_ << "," << size3_ << std::endl;
      SystemModule::my_exit();  
    }
    if( k >= size3_) {
      std::cerr << "ERROR: In Array3D.h - Over bounds on 3rd index. "
		<< "Acessing Array3D::("<< i <<","<< j <<","<< k 
		<< "). Size of Array3D "
		<< "is " << size1_ << ","<< size2_ << "," << size3_ << std::endl;
      SystemModule::my_exit();
    } 
  }// End CheckBounds
  
}; // End class Array3D

#endif
