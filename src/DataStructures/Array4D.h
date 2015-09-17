//-*-c++-*-
#ifndef ARRAY4D_H
#define ARRAY4D_H
//****************************************************************************80
//! \class Array4D
//! \brief  Array4D :  A 4-D array class that we can use for small stuff only
//! \details The ideas in this class are due to
//! A. Katz but have been re-implemented by N. Burgess
//! \nick
//! \version $Rev: 5 $
//! \tparam dataT Template parameter for type of data
//****************************************************************************80
#include "my_incl.h"
#include "SystemUtils/SystemModule.h"

template <typename dataT>
class Array4D
{
  //+++++++++++++++++++++++++++++++ PUBLIC STUFF +++++++++++++++++++++++++++++++
public:
//****************************************************************************80
//!
//! \brief Array4D : Array class default constructor
//! \details
//! \nick
//! \version $Rev: 5 $
//!
//****************************************************************************80
  Array4D(){
    size1 = 0;
    size2 = 0;
    size3 = 0;
    size4 = 0;
    data  = NULL;
  } // End Array4D

//****************************************************************************80
//!
//! \brief Array4D : Array class constructor with size specification
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] n1 The size of the first array dimension
//! \param[in] n2 The size of the second array dimension
//! \param[in] n3 The size of the third array dimension
//! \param[in] n4 The size of the fourth array dimension
//****************************************************************************80
  Array4D(intT n1, intT n2, intT n3, intT n4){
    //---> Set the zie
    size1 = n1;
    size2 = n2;
    size3 = n3;
    size4 = n4;

    //---> Set array to null
    data = NULL;

    //---> Allocate the data in size1 is bigger than 1
    if( size1*size2*size3*size4 > 0) { // check_size
      mem = 
	SystemModule::alloc_mem< dataT, int, double>(data,
						     size1*size2*size3*size4);

      //---> Loop over data and set the value to zero;
      for( intT i = 0; i < size1*size2*size3*size4; i++) {// init_loop
        data[i] = (dataT) -999999999999999999999;
      }// end init_loop

    }// End check_destructor

  } // End Array4D

//****************************************************************************80
//!
//! \brief Array4D : Array4D class move constructor
//! \details
//! \nick
//! \param[in] from_array - rvalue reference that we are "moving" data from
//****************************************************************************80
  Array4D(Array4D<dataT>&& from_array)
  {
    //pilfer other's resource
    size1 = from_array.size1;
    size2 = from_array.size2;
    size3 = from_array.size3;
    size4 = from_array.size4;
    data  = from_array.data;
    mem   = from_array.mem;

    //reset from_array
    from_array.size1 = 0;
    from_array.size2 = 0;
    from_array.size3 = 0;
    from_array.size4 = 0;
    from_array.mem   = 0.0;
    from_array.data  = NULL;
  }

//****************************************************************************80
//!
//! \brief operator= : Array4D class move assignment operator
//! \details
//! \nick
//! \param[in] from_array - rvalue reference that we are "moving" data from
//! \return    this  - returns this object with data moved into it
//****************************************************************************80
  Array4D<dataT>& operator=(Array4D<dataT>&& from_array)
  {
    //--->release the current object's resources
    //--> Delete pointer to data
    if (data != NULL)  delete[] data;

    //--->pilfer from_array's resource
    size1 = from_array.size1;
    size2 = from_array.size2;
    size3 = from_array.size3;
    size4 = from_array.size4;
    data  = from_array.data;
    mem   = from_array.mem;

    //--->reset from_array
    from_array.size1 = 0;
    from_array.size2 = 0;
    from_array.size3 = 0;
    from_array.size4 = 0;
    from_array.mem   = 0.0;
    from_array.data  = NULL;

    return *this;
  }

//****************************************************************************80
//!
//! \brief ~Array4D : Destructor for class Array4D
//! \details
//! \nick
//! \version $Rev: 5 $
//!
//****************************************************************************80
  ~Array4D(){
    //---> Delete pointer to data
    if( size1 > 0) {
      if (data != NULL)  delete[] data;
    }
    //---> Reset data pointer to NULL
    data = NULL;
    //---> Reset the size variable
    size1 = 0;
    size2 = 0;
    size3 = 0;
    size4 = 0;
  } // End ~Array4D

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
//! \param[in] n4 The size of the thrid array dimension
//****************************************************************************80
  void initialize(intT n1, intT n2, intT n3, intT n4){
    //---> Set value
    size1 = n1;
    size2 = n2;
    size3 = n3;
    size4 = n4;
    if ( data == NULL && size1 > 0) { // check_size
      mem = SystemModule::alloc_mem< dataT, int, double>(data,
                                                          size1*size2*size3*size4);

      //---> Loop over data and set the value to zero;
      for( intT i = 0; i < size1*size2*size3*size4; i++) {// init_loop
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
    for(intT i = 0; i < size1*size2*size3*size4; i++) {// set_loop
      data[i] = val;
    }// End set_loop
  }// End set_value

//****************************************************************************80
//!
//! \brief () : Parenthesis operator for acessing array
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] i First index of the array you want to acess with parenthesis
//! \param[in] j Second index of the array you want to access with parenthesis
//! \param[in] k Third index of the arrary you want to access with parenthesis
//! \param[in] l Third index of the arrary you want to access with parenthesis
//****************************************************************************80
  dataT& operator () (intT i, intT j, intT k, intT l)
  {
#ifdef DEV_DEBUG
    this->CheckBounds(i,j,k,l);
#endif
    //---> The return of this operator is the ith reference of data pointer
    return (data[ size4*(size3*(size2*i+j)+k) + l ]);
  } // End

//****************************************************************************80
//!
//! \brief () : Parenthesis operator for acessing array, for const. correctness
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] i First index of the array you want to acess with parenthesis
//! \param[in] j Second index of the array you want to access with parenthesis
//! \param[in] k Third index of the arrary you want to access with parenthesis
//! \param[in] l Fourth index of the arrary you want to access with parenthesis
//****************************************************************************80
  const dataT& operator () (intT i, intT j, intT k, intT l) const
  {
#ifdef DEV_DEBUG
    this->CheckBounds(i,j,k,l);
#endif
    //---> The return of this operator is the ith reference of data pointer
    return (data[ size4*(i*size3 + j*size2 + k) + l ]);
  }

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
      return(size1);
      break;
    case 1:
      return(size2);
      break;
    case 2:
      return(size3);
      break;
    case 3:
      return(size4);
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
  double get_mem( )
  {
    /*---> Return the amount of memory used to store pointer data* to user.
      Remember we stored this in variable mem at allocation */
    return(mem);
  } // End get_mem

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
//! \param[in] l The fourth index you want to access
//****************************************************************************80
  dataT* get_ptr(const intT& i, const intT& j,
                       const intT& k, const intT& l)
  {
#ifdef DEV_DEBUG
    this->CheckBounds(i,j,k,l);
#endif
    return( data + size4*(size3*(size2*i + j) + k) + l );
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
//! \param[in] l The fourth index you want to access
//****************************************************************************80
  dataT const * get_ptr(const intT& i, const intT& j,
                       const intT& k, const intT& l) const
  {
#ifdef DEV_DEBUG
    this->CheckBounds(i,j,k,l);
#endif
    return( data + size4*(size3*(size2*i + j) + k) + l );
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
    return (data + size1*size2*size3*size4);
  }// End end

//****************************************************************************80
//! \brief end : Returns pointer to end of a row of the array
//! \nick
//****************************************************************************80
  dataT const * end() const
  {
    return (data + size1*size2*size3*size4);
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
     stream << var_name << "(" << this->get_size(0) << ", "
            << this->get_size(1) << ", " 
	    << this->get_size(2) << ", " 
	    << this->get_size(3) <<  "):\t "
	    << this->get_mem() << " MB" << std::endl;
     
     //---> However ostringstream is not copyable...so return the string
     return stream.str();

  }
private:
  //+++++++++++++++++++++++++++++++ PRIVATE STUFF ++++++++++++++++++++++++++++++

  intT size1; /*!< Size of dimension 1 */
  intT size2; /*!< Size of dimension 2 */
  intT size3; /*!< Size of dimension 3 */
  intT size4; /*!< Size of dimension 4 */
  dataT* data; /*!< Data pointer for 1-D array */
  double mem; /*!< Amount of memory in megabytes specified for the array */

//****************************************************************************80
//!
//! \brief Array4D : Copy constructor...private so it can't be called, which
//!                  blocks copy construction of this class.
//! \details A copy constructor is a constructor that takes a argument of the
//!      type defined by the class.  In this case Array4D is the class.
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
  Array4D(const Array4D<dataT>&) = delete; //blocked

//****************************************************************************80
//!
//! \brief = : The assignment operator... private so it can't be called, which
//!            blocks assignment of one array to another.
//! \details This is the assignement operator.  We are declaring it private to
//!      block the capability Array4D x = Array4D y type of behavior.
//!      We are purposely blocking this capability by
//!      intersting this function as a private operator with NO CODE.
//!
//!      NOTE: Just because this function does nothing does not make it
//!      acceptable to omit the documentation and clear detailed commentary.
//!      Note the detail above and follow it as a rule.
//! \nick
//! \version $Rev: 5 $
//****************************************************************************80
  Array4D& operator = (const Array4D<dataT>&) = delete; //blocked

//****************************************************************************80
//!
//! \brief CheckBounds: Debug feature - checks the bounds of the array
//! \details
//! \nick
//! \version
//! \param[in] i - i index
//****************************************************************************80
  void CheckBounds(intT i, intT j, intT k, intT l) const {
    if( i >= size1) {
      std::cerr << "ERROR: In Array4D.h - Over bounds on 1st index. "
		<< "Acessing Array4D::("<< i <<","<< j <<","<< k << "," << l  
		<< "). Size of Array4D "
		<< "is " << size1 << ","<< size2 << "," << size3 << "," 
		<< size4 << std::endl;
      SystemModule::my_exit();
    }
    if( j >= size2) {
      std::cerr << "ERROR: In Array4D.h - Over bounds on 2nd index. "
	      	<< "Acessing Array4D::("<< i <<","<< j <<","<< k << "," << l  
		<< "). Size of Array4D "
		<< "is " << size1 << ","<< size2 << "," << size3 << "," 
		<< size4 << std::endl;
      SystemModule::my_exit();
    }
    if( k >= size3) {
      std::cerr << "ERROR: In Array4D.h - Over bounds on 3rd index. "
		<< "Acessing Array4D::("<< i <<","<< j <<","<< k << "," << l  
		<< "). Size of Array4D "
		<< "is " << size1 << ","<< size2 << "," << size3 << "," 
		<< size4 << std::endl;
      SystemModule::my_exit();
    } 
    if( l >= size4) {
      std::cerr << "ERROR: In Array4D.h - Over bounds on 4th index. "
		<< "Acessing Array4D::("<< i <<","<< j <<","<< k << "," << l  
		<< "). Size of Array4D "
		<< "is " << size1 << ","<< size2 << "," << size3 << "," 
		<< size4 << std::endl;
      SystemModule::my_exit(); 
    } 
  }// End CheckBounds
  


}; // End class Array4D

#endif
