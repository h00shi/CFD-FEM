//-*-c++-*-

#ifndef ARRAY2D_H
#define ARRAY2D_H

//****************************************************************************80
//! \class Array2D
//! \brief  Array2D :  A 2-D array class that we can use for small stuff only
//! \details The ideas in this class are due to
//! A. Katz but have been re-implemented by N. Burgess
//! \nick
//! \version $Rev: 5 $
//! \tparam dataT Template parameter for type of data
//****************************************************************************80
#include "my_incl.h"
#include "SystemModule.h"

template <typename dataT>
class Array2D
{
  //+++++++++++++++++++++++++++++++ PUBLIC STUFF +++++++++++++++++++++++++++++++
public:
//****************************************************************************80
//!
//! \brief Array2D : Array class default constructor
//! \details
//! \nick
//! \version $Rev: 5 $
//!
//****************************************************************************80
  Array2D(){
    mem = 0.0;
    size1 = 0;
    size2 = 0;
    data  = NULL;
  } // End Array2D

//****************************************************************************80
//!
//! \brief Array2D : Array class constructor with size specification
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] n1 The size of the first array dimension
//! \param[in] n2 The size of the second array dimension
//****************************************************************************80
  Array2D(intT n1, intT n2){
    //---> Set array to null
    data = NULL;

    //---> Set default values of memory and sizes
    mem = 0.0;
    size1 = 0;
    size2 = 0;

    //---> Allocate the data if size request is bigger than 1
    if( n1*n2 > 0) { // check_size
      //---> Set the sizes
      size1 = n1;
      size2 = n2;
      //---> Allocate
      mem = SystemModule::alloc_mem< dataT, int, double>(data, size1*size2);
      //---> Loop over data and set the value to zero;
      for( intT i = 0; i < size1*size2; i++) {// init_loop
        data[i] = (dataT) 0;
      }// end init_loop
    }// End check_destructor
  } // End Array2D

//****************************************************************************80
//!
//! \brief Array2D : Array2D class move constructor
//! \details
//! \nick
//! \param[in] from_array - rvalue reference that we are "moving" data from
//****************************************************************************80
  Array2D(Array2D<dataT>&& from_array)
  {
    //pilfer other's resource
    size1 = from_array.size1;
    size2 = from_array.size2;
    data  = from_array.data;
    mem   = from_array.mem;

    //reset from_array
    from_array.size1 = 0;
    from_array.size2 = 0;
    from_array.mem   = 0.0;
    from_array.data  = NULL;
  }

//****************************************************************************80
//!
//! \brief operator= : Array2D class move assignment operator
//! \details
//! \nick
//! \param[in] from_array - rvalue reference that we are "moving" data from
//! \return    this  - returns this object with data moved into it
//****************************************************************************80
  Array2D<dataT>& operator=(Array2D<dataT>&& from_array)
  {
    //--->release the current object's resources
    //--> Delete pointer to data
    if (data != NULL)  delete[] data;

    //--->pilfer from_array's resource
    size1 = from_array.size1;
    size2 = from_array.size2;
    data  = from_array.data;
    mem   = from_array.mem;

    //--->reset from_array
    from_array.size1 = 0;
    from_array.size2 = 0;
    from_array.mem   = 0.0;
    from_array.data  = NULL;

    return *this;
  }

//****************************************************************************80
//!
//! \brief ~Array2D : Destructor for class Array2D
//! \details
//! \nick
//! \version $Rev: 5 $
//!
//****************************************************************************80
  virtual ~Array2D(){
    //---> Delete pointer to data
    if(data != NULL)  delete[] data;
    //---> Reset data pointer to NULL
    data = NULL;
    //---> Reset the size variable
    size1 = 0;
    size2 = 0;
    mem   = 0;
  } // End ~Array2D

//****************************************************************************80
//!
//! \brief initialize : A function to allocate the size of the array in the case
//!                    that the size is not known at instantiation.
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] n1 The size of the first array dimension
//! \param[in] n2 The size of the second array dimension
//****************************************************************************80
  void initialize(intT n1, intT n2){
    /*-->only perform initialization if the data is null
    and input size request is positive*/
    if ( data == NULL && n1*n2 > 0) { // check_size
      //---> Set the sizes
      size1 = n1;
      size2 = n2;
      //---> Allocate
      mem = SystemModule::alloc_mem< dataT, int, double>(data, size1*size2);
      //---> Loop over data and set the value to zero;
      for( intT i = 0; i < size1*size2; i++) {// init_loop
        data[i] = (dataT) 0;
      }// end init_loop
    } // End check_size
  } // End initialize

//****************************************************************************80
//!
//! \brief set_value : Sets the whole array to a value
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] val The value you want to set the whole array to
//****************************************************************************80
  void set_value(const dataT& val)
  {
    for(intT i = 0; i < size1*size2; i++) {// set_loop
      data[i] = val;
    }// End set_loop
  }// End set_value

//****************************************************************************80
//!
//! \brief () : Parenthesis operator for acessing array, for const. correctness
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] i First index of the array you want to access with parenthesis
//! \param[in] j Second index of the array you want to access with parenthesis
//****************************************************************************80
  const dataT& operator () (intT i, intT j) const
  {
#ifdef DEV_DEBUG
    this->CheckBounds(i,j);
#endif

    //---> The return of this operator is the ith reference of data pointer
    return (data[i*size2 + j]);
  }

//****************************************************************************80
//!
//! \brief () : Parenthesis operator for acessing array
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] i First index of the array you want to access with parenthesis
//! \param[in] j Second index of the array you want to access with parenthesis
//****************************************************************************80
  dataT& operator () (intT i, intT j)
  {
#ifdef DEV_DEBUG
    this->CheckBounds(i,j);
#endif

    //---> The return of this operator is the ith reference of data pointer
    return (data[i*size2 + j]);
  } // End

//****************************************************************************80
//!
//! \brief get_dims  : For the specified index gets the size of that dimension
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] dim The dimension you want the size of
//****************************************************************************80
  intT get_size(intT dim) const
  {
    //---> Check the value of dim to return the correct size
    switch (dim) { // check_dim
    case 0:
      return(size1);
      break;
    case 1:
      return(size2);
      break;
    default:
      return(-99);
    } // End check_dim

  }// End get_size

//****************************************************************************80
//!
//! \brief get_mem : Diagnostic routine to query how much memory the class is
//!                  using for the data part of the class
//! \details
//! \nick
//! \version $Rev: 5 $
//!
//****************************************************************************80
  double get_mem( ) const
  {
    /*---> Return the amount of memory used to store pointer data* to user.
      Remember we stored this in variable mem at allocation */
    return(mem);
  } // End get_mem

//****************************************************************************80
//!
//! \brief operator << : A facility to print the Array2D class the screen
//! \details
//! \nick
//! \version $Rev$
//! \param[in] os The ostream operator we stream to
//! \param[in] a The array we wish to stream
//! \return os Returns the modified value of os
//****************************************************************************80
  friend std::ostream& operator << (std::ostream& os, const Array2D<dataT>& a)
  {
    //---> Write the data to ostream object
    for (intT i = 0; i < a.size1; i++) {
      os << i << ": ";
      for (intT j = 0; j < a.size2; j++){
        os << a(i,j) << " ";
      }
      os << std::endl;
    }

    return(os);
  } // End operator <<

//****************************************************************************80
//!
//! \brief get_ptr : Allow access to pointer at specified address
//! \details Constant correct version
//! \nick
//! \version
//! \param[in] i The first address you want to access
//! \param[in] j The second address you want to access
//****************************************************************************80
  dataT const * get_ptr(const intT& i, const intT& j) const
  {
#ifdef DEV_DEBUG
    this->CheckBounds(i,j);
#endif
    return (data + i*size2 + j);

  }// End get_ptr

//****************************************************************************80
//! \brief get_ptr : Allow acess to pointer at specified address
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] i The first address you want to access
//! \param[in] j The second address you want to access
//****************************************************************************80
  dataT* get_ptr(const intT& i, const intT& j)
  {
#ifdef DEV_DEBUG
    this->CheckBounds(i,j);
#endif
    return(data + i*size2 + j);
  }// End get_ptr

//****************************************************************************80
//! \brief begin : Returns pointer to start of the array
//! \nick
//****************************************************************************80
  dataT* begin()
  {
    return (data);
  }// End begin

//****************************************************************************80
//! \brief begin : Returns pointer to start of a row of the array
//! \nick
//****************************************************************************80
  dataT const * begin() const
  {
    return (data);
  }// End begin

//****************************************************************************80
//! \brief end : Returns pointer to end of a row of the array
//! \nick
//****************************************************************************80
  dataT* end()
  {
    return (data + size1*size2);
  }// End end

//****************************************************************************80
//! \brief end : Returns pointer to end of a row of the array
//! \nick
//****************************************************************************80
  dataT const * end() const
  {
    return (data + size1*size2);
  }// End end

private:
  //+++++++++++++++++++++++++++++++ PRIVATE STUFF ++++++++++++++++++++++++++++++

  intT size1; /*!< Size of dimension 1 */
  intT size2; /*!< Size of dimension 2 */
  dataT* data; /*!< Data pointer for 1-D array */
  double mem; /*!< Amount of memory in megabytes specified for the array */

//****************************************************************************80
//!
//! \brief Array2D : Copy constructor...private so it can't be called, which
//!                  blocks copy construction of this class.
//! \details A copy constructor is a constructor that takes a argument of the
//!      type defined by the class.  In this case Array2D is the class.
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
  Array2D(const Array2D<dataT>&) = delete; //blocked

//****************************************************************************80
//!
//! \brief = : The assignment operator... private so it can't be called, which
//!            blocks assignment of one array to another.
//! \details This is the assignement operator.  We are declaring it private to
//!      block the capability Array2D x = Array2D y type of behavior.
//!      We are purposely blocking this capability by
//!      intersting this function as a private operator with NO CODE.
//!
//!      NOTE: Just because this function does nothing does not make it
//!      acceptable to omit the documentation and clear detailed commentary.
//!      Note the detail above and follow it as a rule.
//! \nick
//! \version $Rev: 5 $
//****************************************************************************80
  Array2D& operator = (const Array2D<dataT>&) = delete; //blocked

//****************************************************************************80
//!
//! \brief CheckBounds: Debug feature - checks the bounds of the array
//! \details
//! \nick
//! \version
//! \param[in] i - i index
//****************************************************************************80
  void CheckBounds(intT i, intT j) const {
    if( i >= size1) {
      std::cerr << "ERROR: In Array2D.h - Over bounds on 1st index. "
		<< "Acessing Array2D::("<< i <<","<< j <<"). Size of Array2D "
		<< "is " << size1 << ","<< size2 << std::endl;
      SystemModule::my_exit();
    }
    if( j >= size2) {
      std::cerr << "ERROR: In Array2D.h - Overbounds on 2nd index. "
		<< "Acessing Array2D::("<< i <<","<< j <<"). Size of Array2D "
		<< "is " << size1 << ","<< size2 << std::endl;
      SystemModule::my_exit();
    }
      
  }// End CheckBounds

}; // End class Array2D

#endif
