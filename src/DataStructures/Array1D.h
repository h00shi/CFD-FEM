//-*-c++-*-

#ifndef ARRAY1D_H
#define ARRAY1D_H
//****************************************************************************80
//! \class Array1D
//! \brief  Array1D :  A 1-D array class that we can use for small stuff only
//! \details The ideas in this class are due to
//! A. Katz but have been re-implemented by N. Burgess
//! \nick
//! \version $Rev: 5 $
//! \tparam dataT Template parameter for type of data
//****************************************************************************80
#include "my_incl.h"
#include "SystemModule.h"

template <typename dataT>
class Array1D
{
  //+++++++++++++++++++++++++++++++ PUBLIC STUFF +++++++++++++++++++++++++++++++
public:
//****************************************************************************80
//!
//! \brief Array1D : Array class default constructor
//! \details
//! \nick
//! \version $Rev: 5 $
//!
//****************************************************************************80
  Array1D(){
    mem = 0.0;
    size1 = 0;
    data  = NULL;
  } // End Array1D

//****************************************************************************80
//!
//! \brief Array1D : Array class constructor with size specification
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] n The size of the array you are declaring
//****************************************************************************80
  Array1D(intT n){
    //---> Set array to null
    data = NULL;

    //---> Set default values of memory and size
    mem = 0.0;
    size1 = n;

    //---> Allocate the data in size1 is bigger than 1
    if( n > 0) { // check_size
      //---> Set the size
      size1 = n;
      //---> Allocate
      mem = SystemModule::alloc_mem< dataT, intT, double>(data, size1);
      //---> Loop over data and set the value to zero;
      for( intT i = 0; i < size1; i++) {// init_loop
        data[i] = (dataT) 0;
      }// end init_loop

    }// End check_destructor
  } // End Array1D

//****************************************************************************80
//!
//! \brief Array1D : Array class move constructor
//! \details
//! \nick
//! \param[in] other - rvalue reference that we are "moving" data from
//****************************************************************************80
  Array1D(Array1D<dataT>&& from_array)
  {


    //pilfer other's resource
    size1 = from_array.size1;
    data  = from_array.data;
    mem   = from_array.mem;

    //reset other
    from_array.size1 = 0;
    from_array.mem   = 0.0;
    from_array.data  = NULL;
  }

//****************************************************************************80
//!
//! \brief operator= : Array class move assignment operator
//! \details
//! \nick
//! \param[in] other - rvalue reference that we are "moving" data from
//! \return    this  - returns this object with data moved into it
//****************************************************************************80
  Array1D<dataT>& operator=(Array1D<dataT>&& from_array)
  {
    //--->release the current object's resources
    //--> Delete pointer to data
    if (data != NULL)  delete[] data;

    //--->pilfer from_array's resource
    size1 = from_array.size1;
    data  = from_array.data;
    mem   = from_array.mem;

    //--->reset from_array
    from_array.size1 = 0;
    from_array.mem   = 0.0;
    from_array.data  = NULL;

    return *this;
  }

//****************************************************************************80
//!
//! \brief ~Array1D : Destructor for class Array1D
//! \details
//! \nick
//! \version $Rev: 5 $
//!
//****************************************************************************80
  ~Array1D(){
    //---> Delete pointer to data
    if (data != NULL)  delete[] data;
    //---> Reset data pointer to NULL
    data = NULL;
    //---> Reset the size variable
    size1 = 0;
  } // End ~Array1D

//****************************************************************************80
//!
//! \brief initialize : A function to allocate the size of the array in the case
//!                    that the size is not known at instantiation.
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] n The size of the array to be initialized
//****************************************************************************80
  void initialize(intT n){
    /*-->only perform initialization if the data is null
    and input size request is positive*/
    if ( data == NULL && n > 0) {
      //---> Set value
      size1 = n;
      //---> Allocate
      mem = SystemModule::alloc_mem< dataT, intT, double>(data, size1);
      //---> Loop over data and set the value to zero;
      for( intT i = 0; i < size1; i++) {// init_loop
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
    for(intT i = 0; i < size1; i++) {// set_loop
      data[i] = val;
    }// End set_loop
  }// End set_value

//****************************************************************************80
//!
//! \brief () : Parenthesis operator for acessing array, for const. correctness
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] i The index of the array you want to acess with parenthesis
//****************************************************************************80
const dataT& operator () (intT i) const
  {

#ifdef DEV_DEBUG
    this->CheckBounds(i);
#endif

    //---> The return of this operator is the ith reference of data pointer
    return (data[i]);
  } // End ()

//****************************************************************************80
//!
//! \brief () : Parenthesis operator for acessing array
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] i The index of the array you want to acess with parenthesis
//****************************************************************************80
  dataT& operator () (intT i)
  {
    //calls constant version of operator (i)
    return const_cast<dataT&>(static_cast<const Array1D&>(*this).
                              operator()(i));
  } // End ()

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
    //---> It's only a 1-D array so just give back the size
    switch (dim) {
    case 0: 
      return(size1);
      break;
    default :
      return(size1);
    }
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
  double get_mem() const
  {
    /*---> Return the amount of memory used to store pointer data* to user.
      Remember we stored this in variable mem at allocation */
    return(mem);
  } // End get_mem

//****************************************************************************80
//!
//! \brief operator << : A facility to print the Array1D class the screen
//! \details
//! \nick
//! \version $Rev$
//! \param[in] os The ostream operator we stream to
//! \param[in] a The array we wish to stream
//! \return os Returns the modified value of os
//****************************************************************************80
  friend std::ostream& operator << (std::ostream& os, const Array1D<dataT>& a)
  {
    //---> Write the data to ostream object
    for (intT i = 0; i < a.size1; i++) {
      os << i << ": " << a(i) << std::endl;;
    }
    return(os);
  } // End operator <<

//****************************************************************************80
//!
//! \brief get_ptr : Allow access to pointer at specified address
//! \details  Constant correct version
//! \nick
//! \version
//! \param[in] i The first address you want to access
//****************************************************************************80
  dataT const * get_ptr(const intT& i) const
  {
#ifdef DEV_DEBUG
    this->CheckBounds(i);
#endif

    return(data + i);
  }// End get_ptr

//****************************************************************************80
//!
//! \brief get_ptr : Allow acess to pointer at specified address
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] i The address you want to access
//****************************************************************************80
  dataT* get_ptr(const intT& i)
  {
#ifdef DEV_DEBUG
    this->CheckBounds(i);
#endif
    //---> Just return the point + i 
    return(data + i);
  }// End get_ptr

//****************************************************************************80
//! \brief begin : Returns pointer to start of array
//! \nick
//****************************************************************************80
  dataT* begin()
  {
    return(data);
  }// End begin

//****************************************************************************80
//!
//! \brief begin : Returns pointer to start of array
//! \nick
//****************************************************************************80
  dataT const * begin() const
  {
    return(data);
  }// End begin

//****************************************************************************80
//!
//! \brief end : Returns pointer to end of array
//! \nick
//****************************************************************************80
  dataT* end()
  {
    return(data + size1);
  }// End end

//****************************************************************************80
//!
//! \brief end : Returns pointer to end of array
//! \nick
//****************************************************************************80
  dataT const * end() const
  {
    return(data + size1);
  }// End end

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
     stream << var_name << "(" << this->get_size(0) << "):\t "
     	   << this->get_mem() << " MB" << std::endl;

     //---> However ostringstream is not copyable...so return the string
     return stream.str();

  }

private:
  //+++++++++++++++++++++++++++++++ PRIVATE STUFF ++++++++++++++++++++++++++++++

  intT size1; /*!< Size of dimension 1 */
  dataT* data; /*!< Data pointer for 1-D array */
  double mem; /*!< Amount of memory in megabytes specified for the array */

//****************************************************************************80
//!
//! \brief Array1D : Copy constructor...private so it can't be called, which
//!                  blocks copy construction of this class.
//! \details A copy constructor is a constructor that takes a argument of the
//!      type defined by the class.  In this case Array1D is the class.
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
  Array1D(const Array1D<dataT>&) = delete; //blocked

//****************************************************************************80
//!
//! \brief = : The assignment operator... private so it can't be called, which
//!            blocks assignment of one array to another.
//! \details This is the assignement operator.  We are declaring it private to
//!      block the capability Array1D x = Array1D y type of behavior.
//!      We are purposely blocking this capability by
//!      intersting this function as a private operator with NO CODE.
//!
//!      NOTE: Just because this function does nothing does not make it
//!      acceptable to omit the documentation and clear detailed commentary.
//!      Note the detail above and follow it as a rule.
//! \nick
//! \version $Rev: 5 $
//****************************************************************************80
  Array1D& operator = (const Array1D<dataT>&) = delete; //blocked

  
//****************************************************************************80
//!
//! \brief CheckBounds: Debug feature - checks if the input ith IX is allowable
//! \details
//! \nick
//! \version
//! \param[in] i - i index
//****************************************************************************80
  void CheckBounds(intT i) const{
    
    if(i >= size1){
      std::cerr << "ERROR: In Array1D.h - Over bounds on 1st index. "
                << "Accessing Array1D::("<< i <<"). Size of Array1D is " 
		<< size1 << std::endl;
      SystemModule::my_exit();
    }
  } // End CheckBounds

}; // End class Array1D

#endif
