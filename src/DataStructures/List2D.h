//-*-c++-*-
#ifndef LIST2D_H
#define LIST2D_H

//****************************************************************************80
//! \class List2D
//! \brief  List2D : A 2-D ordered list class.  Think of it like 
//!                  can non-uniform 2D array.
//! \nick
//! \version $Rev: 5 $
//! \tparam dataT Template parameter for type of data
//****************************************************************************80
#include "my_incl.h"
#include "SystemModule.h"
#include "Array1D.h"

template <class dataT>
class List2D
{
  //+++++++++++++++++++++++++++++++ PUBLIC STUFF +++++++++++++++++++++++++++++++
public:
//****************************************************************************80
//!
//! \brief List2D : Array class default constructor
//! \details
//! \nick
//! \version $Rev: 5 $
//!
//****************************************************************************80
  List2D(){
    mem = 0.0;
    size1 = 0;
    size = 0;
    index = NULL;
    data  = NULL;
  } // End List2D

//****************************************************************************80
//!
//! \brief List2D : Array class constructor with size specification
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] n1 The size of the leading list dimension
//! \param[in] n Total size of the all the data
//****************************************************************************80
  List2D(intT nrow, intT n){
    //---> Set array to null
    index = NULL;
    data  = NULL;

    //---> Set default values of memory and sizes
    mem   = 0.0;
    size1 = nrow;
    size  = n;

    //---> Allocate index if requested number of rows is positive
    if(nrow > 0) { //check row request
      //---> Set the size
      size1 = nrow;
      //---> Allocate
      mem += SystemModule::alloc_mem< int, int, double>(index, size1 + 1);
      //---> Set up columns of each row to be zero
      index[0] = 0;
      for(intT i = 0; i < nrow; i++){
        index[i + 1] = index[i] + 0;
      }
    }//end check row request

    //---> Allocate data if requested number of total size is positive
    if(n > 0) { // check input total size
      //---> Set the size
      size = n;
      //---> Allocate
      mem += SystemModule::alloc_mem< dataT, int, double>(data, size);
      //---> Loop over data and set the value to zero;
      for( intT i = 0; i < size; i++) {// init_loop
        data[i] = (dataT) 0;
      }// end init_loop
    }//end check total size

  } // End List2D

//****************************************************************************80
//!
//! \brief List2D : List2D class move constructor
//! \details
//! \nick
//! \param[in] other - rvalue reference that we are "moving" data from
//****************************************************************************80
  List2D(List2D<dataT>&& other)
  {
    //--->pilfer other's resource
    size1 = other.size1;
    size  = other.size;
    index = other.index;
    data  = other.data;
    mem   = other.mem;

    //--->reset other
    other.size1 = 0;
    other.size  = 0;
    other.mem   = 0.0;
    other.index = NULL;
    other.data  = NULL;
  }

//****************************************************************************80
//!
//! \brief operator= : List2D class move assignment operator
//! \details
//! \nick
//! \param[in] other - rvalue reference that we are "moving" data from
//! \return    this  - returns this object with data moved into it
//****************************************************************************80
  List2D<dataT>& operator=(List2D<dataT>&& other)
  {
    //--->release the current object's resources
    //--> Delete pointer to data
    if (index !=NULL)  delete[] index;
    if (data != NULL)  delete[] data;

    //--->pilfer other's resource
    size1 = other.size1;
    size  = other.size;
    index = other.index;
    data  = other.data;
    mem   = other.mem;

    //--->reset other
    other.size1 = 0;
    other.size  = 0;
    other.mem   = 0.0;
    other.index = NULL;
    other.data  = NULL;

    return *this;
  }

//****************************************************************************80
//!
//! \brief ~List2D : Destructor for class List2D
//! \details
//! \nick
//! \version $Rev: 5 $
//!
//****************************************************************************80
  ~List2D(){
    //---> Delete pointer to data if it exists
    if(index !=NULL) delete[] index;
    if(data != NULL) delete[] data;

    //---> Reset data pointer to NULL
    index = NULL;
    data = NULL;
    //---> Reset the size variable
    size1 = 0;
    size = 0;
    mem = 0.0;
  } // End ~List2D
//****************************************************************************80
//!
//! \brief CopyPattern : Take an existing List2D and copy it's pattern
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] other The other list 2D
//****************************************************************************80
  template<class otherT>
  void intialize_copy_pattern(const List2D<otherT>& other)
  {
    if (this->index == NULL && this->size1 == 0){
      //---> Call the initialize function
      this->initialize(other.get_lead_size(), other.get_total_size());
      //---> Set the number of columns to be the same
      for(intT i = 0; i < size1; i++) {this->set_ncol(i,other.get_ncol(i));}
    }
    
    else {
      std::cerr << "ERROR: In List2D.h - "
                << "Attemping to initalize and copy pattern for a List2D "
		<< "that has already been intialized.  List2D data follows."
		<< std::endl 
                << "This List2D, nrow: " << size1 << " nnz: " << size 
		<< std::endl
		<< "Attemping to pattern from List2D, nrow: " << size1 
		<< " nnz: " << size << std::endl;
      SystemModule::my_exit();
      
    }

  }// End CopyPattern;

//****************************************************************************80
//!
//! \brief initialize : A function to allocate the size of the list in the case
//!                    that the size is not known at instantiation.
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] n1 The size of the leading list dimension
//! \param[in] n Total size of the all the data
//****************************************************************************80
  void initialize(intT nrow, intT n){
    /*-->only perform index initialization if the index is null
    and input row request is positive*/
    if (index == NULL && nrow > 0){
      //---> Set the size
      size1 = nrow;
      //---> Allocate
      mem += SystemModule::alloc_mem< int, int, double>(index, size1 + 1);
      //---> Set up columns of each row to be zero
      index[0] = 0;
      for(intT i = 0; i < size1; i++){
        index[i + 1] = index[i] + 0;
      }
    }// end check input sizes

    /*-->only perform data initialization if the data is null
    and input total size request is positive*/
    if (data  == NULL && n > 0) {
      //---> Set the size
      size = n;
      //---> Allocate
      mem += SystemModule::alloc_mem< dataT, int, double>(data, size);
      //---> Loop over data and set the value to zero;
      for( intT i = 0; i < size; i++) {// init_loop
        data[i] = (dataT) 0;
      }// end init_loop
    }// end check input sizes
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
    for(intT i = 0; i < size; i++) {// set_loop
      data[i] = val;
    }// End set_loop
  }// End set_value

//****************************************************************************80
//! \brief setup_ncol : Set number of columns and create index array for row
//! \details NOTE: must be called in lexigraphical order!
//!                JKU - i think this function is scary.
//! \nick
//! \version $Rev$
//! \param[in] row The row number having it's number of columns set
//****************************************************************************80
  void set_ncol(intT row, intT nc)
  {
#ifdef DEV_DEBUG
    this->CheckSizeIndex(row);
#endif

    index[row + 1] = index[row] + nc;
  } // end set_ncol

//****************************************************************************80
//! \brief set_ncol : Sets number of columns and creates index array for all
//!                   rows.
//! \details This function modifies the index array such that the number of
//!          columns of row i is equal to ncol(i).
//!          Array1D ncol should be the same size of size1.
//! \nick
//! \version $Rev$
//! \param[in] ncol An Array1D of requested column numbers (signed int)
//****************************************************************************80
  void set_ncol(Array1D<int> const & ncol)
  {
#ifdef DEV_DEBUG
    //debug check - checking if size(Array1D) == nrow
    if(ncol.get_size(0) != size1){
      std::cerr << "In List2D.h - "
                << "ERROR: Mismatch of num rows and num of col specification."
                << "\n Number of rows: "<< size1
                << "\n Number of ncol specs: " << ncol.get_size(0)
                << std::endl;
      SystemModule::my_exit();
    }
#endif

    if(size1 > 0) {//check if we have columns
      //setup index array based on ncol input
      index[0] = 0;
      for(intT irow = 0; irow < size1; irow++){
        //---> Index;
        index[irow+1] = index[irow] + ncol(irow);
      }

#ifdef DEV_DEBUG
      //debug check - checking if total num of entries <= entire size of List
      intT targetSize = index[size1] - index[0];

      if(targetSize > size){
        std::cerr << "In List2D.h - "
                  << "ERROR: Requesting too many elements."
                  << "\n Target total size of requested columns: " << targetSize
                  << "\n Preallocated total size: " << size << std::endl;
        SystemModule::my_exit();
      }
#endif
    }//end check if we have columns
  } // end set_ncol

//****************************************************************************80
//! \brief get_ncol : Get the number of columns for a given row
//! \details
//! \nick
//! \version $Rev$
//! \param[in] row The row who's number of columns is desired
//****************************************************************************80
  intT get_ncol(intT row) const
  {

#ifdef DEV_DEBUG
    this->CheckSizeIndex(row);
#endif

    return( index[row + 1] - index[row] );
  }// end get_ncol

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
    return (data[index[i] + j]);
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
    return (data[index[i] + j]);
  } // End

//****************************************************************************80
//!
//! \brief () : Parenthesis operator for acessing array, for const. correctness
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] i First index of the array you want to access with parenthesis
//! \param[in] j Second index of the array you want to access with parenthesis
//****************************************************************************80
  const dataT& operator () (intT i) const
  {
#ifdef DEV_DEBUG
    this->CheckSize(i);
#endif

    //---> The return of this operator is the ith reference of data pointer
    return (data[i]);
  }

//****************************************************************************80
//!
//! \brief () : Parenthesis operator for accessing array with a single value
//! \details  This accessor treats the List2D object as a 1D array.
//! \nick
//! \version
//! \param[in] i Index of the array you want to access with parenthesis
//****************************************************************************80
  dataT& operator () (intT i)
  {
    //calls constant version of operator (i)
    return const_cast<dataT&>(static_cast<const List2D &>(*this).operator()(i));
  } // End

//****************************************************************************80
//!
//! \brief get_lead_size : Get the size of leading index of list
//! \details
//! \nick
//! \version $Rev: 5 $
//****************************************************************************80
  intT get_lead_size( ) const
  {
    return((int)size1);
  }// End get_lead_size

//****************************************************************************80
//!
//! \brief get_total_size : Get the total size of the list
//! \details
//! \nick
//! \version $Rev: 5 $
//****************************************************************************80
  intT get_total_size( ) const
  {
    return((int)size);
  }// End get_lead_size

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
//! \brief operator << : A facility to print the List2D class the screen
//! \details
//! \nick
//! \version $Rev$
//! \param[in] os The ostream operator we stream to
//! \param[in] a The array we wish to stream
//! \return os Returns the modified value of os
//****************************************************************************80
  friend std::ostream& operator << (std::ostream& os, const List2D<dataT>& a)
  {
    //---> Write the data to ostream object
    for (intT i = 0; i < a.size1; i++) {
      os << i << ": ";
      for (intT j = 0; j < a.get_ncol(i); j++){
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
    return(data + index[i] + j);
  }// End get_ptr

//****************************************************************************80
//!
//! \brief get_ptr : Allow acess to pointer at specified address
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] i The first address you want to access
//! \param[in] j The second address you want to access
//****************************************************************************80
  dataT* get_ptr(const intT& i, const intT& j)
  {

#   ifdef DEV_DEBUG
    this->CheckBounds(i,j);
#endif
    return(data + index[i] + j);
    
  }// End get_ptr

//****************************************************************************80
//!
//! \brief set_ptr : Set the pointer to some "external" thing
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] x Pointer you want to assign to data
//****************************************************************************80
  void set_ptr(dataT* x) {
    data = x;
  }

//****************************************************************************80
//! \brief get_index_ptr : Returns the pointer to the indexing array
//! \details
//! \nick
//! \version $Rev$
//! \param[in] i The leading address
//****************************************************************************80
  intT const * get_index_ptr(const intT& i ) const
  {
#ifdef DEV_DEBUG
    this->CheckSizeIndex(i);
#endif
    return(index + i);
  }

//****************************************************************************80
//! \brief get_index_ptr : Returns the pointer to the indexing array
//! \details
//! \nick
//! \version $Rev$
//! \param[in] i The leading address
//****************************************************************************80
  int* get_index_ptr(const intT& i )
  {
    //calls constant version of get_index_ptr
    return const_cast<int*>(static_cast<const List2D&>(*this).get_index_ptr(i));
  }

//****************************************************************************80
//! \brief get_index : Return global index to point i,j in the list
//! \details
//! \nick
//! \version $Rev$
//! \param[in] i Leading dimension index
//! \param[in] j Second dimension index
//****************************************************************************80
  intT get_index(const intT& i, const intT& j) const
  {
#ifdef DEV_DEBUG
    this->CheckBounds(i,j);
#endif
    return( index[i] + j );
  }

//****************************************************************************80
//! \brief begin : Returns pointer to start of the list
//! \nick
//****************************************************************************80
  dataT* begin()
  {
    return(data);
  }// End begin

//****************************************************************************80
//! \brief begin : Returns pointer to start of a row of the list
//! \nick
//****************************************************************************80
  dataT const * begin() const
  {
    return(data);
  }// End begin

//****************************************************************************80
//! \brief end : Returns pointer to end of a row of the list
//! \nick
//****************************************************************************80
  dataT* end()
  {
    if(index !=NULL){
      return(data + index[size1]);
    }else{
      return(data);
    }
  }// End end

//****************************************************************************80
//!
//! \brief end : Returns pointer to end of a row of the list
//! \nick
//****************************************************************************80
  dataT const * end() const
  {
    if(index !=NULL){
      return(data + index[size1]);
    }else{
      return(data);
    }
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
     stream << var_name << "(" << this->get_lead_size() << ")(" 
	    << this->get_total_size() << "):\t "
	    << this->get_mem() << " MB" << std::endl;
     
     //---> However ostringstream is not copyable...so return the string
     return stream.str();

  }
private:
  //+++++++++++++++++++++++++++++++ PRIVATE STUFF ++++++++++++++++++++++++++++++
  intT size1; /*!< Size of the leading dimension of the list */
  intT size; /*!< Total size of the array */
  int* index; /*!< Index pointer that allows for addressing 2D */
  dataT* data; /*!< Data pointer */
  double mem; /*!< Amount of memory in megabytes specified for the array */

//****************************************************************************80
//! \brief List2D : Copy constructor...private so it can't be called, which
//!                  blocks copy construction of this class.
//! \details A copy constructor is a constructor that takes a argument of the
//!      type defined by the class.  In this case List2D is the class.
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
  List2D(const List2D<dataT>&) = delete; //blocked

//****************************************************************************80
//!
//! \brief = : The assignment operator... private so it can't be called, which
//!            blocks assignment of one array to another.
//! \details This is the assignement operator.  We are declaring it private to
//!      block the capability List2D x = List2D y type of behavior.
//!      We are purposely blocking this capability by
//!      intersting this function as a private operator with NO CODE.
//!
//!      NOTE: Just because this function does nothing does not make it
//!      acceptable to omit the documentation and clear detailed commentary.
//!      Note the detail above and follow it as a rule.
//! \nick
//! \version $Rev: 5 $
//****************************************************************************80
  List2D& operator = (const List2D<dataT>&) = delete; // Blocked
  

//****************************************************************************80
//!
//! \brief CheckSize: Debug feature - checks if the input index is allowable
//!                     when list is accessed like a 1D array
//! \details
//! \nick
//! \version
//! \param[in] i - 1D access index
//****************************************************************************80
  void CheckSize(intT i) const{
    if(i >= size){
      std::cerr << "ERROR: In List2D.h - "
                << "Attempting to access index "<< i
                << ".   Number of entries is " << size << std::endl;
      SystemModule::my_exit();
    }
  }// End CheckSize

//****************************************************************************80
//!
//! \brief CheckBounds: Debug feature - checks the bounds of the list
//! \details
//! \nick
//! \version
//! \param[in] i - i index
//****************************************************************************80
  void CheckBounds(intT i, intT j) const{
    if(i >= size1){
      std::cerr << "ERROR: In List2D.h - Over bounds on 1st index. "
		<< "Acessing List2D::("<< i <<","<< j <<"). Size is "
		<< "is " << size1 << ","<< get_ncol(i) << std::endl;
      SystemModule::my_exit();
      
    }
    if( j >= get_ncol(i) ) {
      std::cerr << "ERROR: In List2D.h - Over bounds on 2nd index. "
		<< ".   Number of columns of row " << i
                << " is " << get_ncol(i) << std::endl;
      SystemModule::my_exit();
    }
  }// End CheckBounds


//****************************************************************************80
//!
//! \brief CheckSizeIndex: Debug feature - checks if the input index is 
//!         allowable when list is accessed like a 1D array
//! \details
//! \nick
//! \version
//! \param[in] ix - 1D access index
//****************************************************************************80
  void CheckSizeIndex(intT i) const{
    if(i >= size1 ){
      std::cerr << "ERROR: In List2D.h - "
                << "Attempting to access index "<< i
                << ".   Size of indexing array: " << size1 << std::endl;
      SystemModule::my_exit();
    }
  }// End CheckSizeIndex

}; // End class List2D
#endif
