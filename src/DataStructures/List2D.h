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
#include "SystemUtils/SystemModule.h"
#include "DataStructures/Array1D.h"
#include <sstream>
template <class dataT>
class List2D : public Array1D<dataT>
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
    size1_ = 0;
    size_ = 0;
    index_ = nullptr;
    mem_ = 0.0;
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
  List2D(intT nrow, intT n): List2D(){
    this->initialize(nrow,n);
  } // End List2D

//****************************************************************************80
//!
//! \brief List2D constructor that takes an Array1D with the number of columns
//!        per row;
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] ncol The number of columns per row
//****************************************************************************80
  List2D(const Array1D<intT>& ncol ) : List2D() {
    this->initialize(ncol);
  }// End List2D
//****************************************************************************80
//!
//! \brief List2D : List2D class move constructor
//! \details
//! \nick
//! \param[in] other - rvalue reference that we are "moving" data from
//****************************************************************************80
  List2D(List2D<dataT>&& other) : Array1D<dataT>(std::move(other))
  {
    //--->pilfer other's resource
    size1_ = other.size1_;
    size_  = other.size_;
    index_ = other.index_;
    mem_ = other.mem_;

    //--->reset other
    other.size1_ = 0;
    other.size_  = 0;
    other.index_ = NULL;
    other.mem_ = 0.0;
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
    Array1D<dataT>::operator=(std::move(other));
    //--->pilfer other's resource
    size1_ = other.size1_;
    size_  = other.size_;
    index_ = other.index_;
    mem_ = other.mem_;

    //--->reset other
    other.size1_ = 0;
    other.size_  = 0;
    other.index_ = nullptr;
    other.mem_ = 0.0;
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
    if(index_ !=NULL) delete[] index_;

    //---> Reset data pointer to NULL
    index_ = NULL;

    //---> Reset the size variable
    size1_ = 0;
    size_ = 0;
    mem_ = 0.0;
  } // End ~List2D
//****************************************************************************80
//!
//! \brief CopyPattern : Take an existing List2D and copy it's pattern
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] other The other list 2D
//****************************************************************************80
  template <class otherT>
  inline void initialize_copy_pattern(const List2D<otherT>& other)
  {
    if (this->index_ == NULL && this->size1_ == 0){
      //---> Call the initialize function
      this->initialize(other.get_lead_size(), other.get_total_size());
      //---> Set the number of columns to be the same
      for(intT i = 0; i < size1_; i++) {this->set_ncol(i,other.get_ncol(i));}
    }
    
    else {
      std::cerr << "ERROR: In List2D.h - "
                << "Attemping to initalize and copy pattern for a List2D "
		<< "that has already been intialized.  List2D data follows."
		<< std::endl 
                << "This List2D, nrow: " << size1_ << " nnz: " << size_
		<< std::endl
		<< "Attemping to pattern from List2D, nrow: " << size1_
		<< " nnz: " << size_ << std::endl;
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
  inline void initialize(intT nrow, intT n){
    /*-->only perform index initialization if the index is null
    and input row request is positive*/
    if (index_ == NULL && nrow > 0){
      //---> Set the size
      size1_ = nrow;
      size_ = n;
      //---> Allocate
      mem_ = SystemModule::alloc_mem< int, int, double>(index_, size1_ + 1);
      //---> Set up columns of each row to be zero
      index_[0] = 0;
      for(intT i = 0; i < size1_; i++){
        index_[i + 1] = index_[i] + 0;
      }
    }// end check input sizes

    Array1D<dataT>::initialize(size_);
  } // End initialize

//****************************************************************************80
//!
//! \brief Initialization that takes an Array1D with the number of columns
//!        per row;
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] ncol The number of columns per row
//****************************************************************************80
  inline void initialize(const Array1D<intT>& ncol ) {

    intT nnz = 0;
    for(intT i = 0; i < ncol.get_size(0); i++){nnz += ncol(i);}
    this->initialize(ncol.get_size(0), nnz);
    this->set_ncol(ncol);

  }// End List2D

//****************************************************************************80
//! \brief setup_ncol : Set number of columns and create index array for row
//! \details NOTE: must be called in lexigraphical order!
//!                JKU - i think this function is scary.
//! \nick
//! \version $Rev$
//! \param[in] row The row number having it's number of columns set
//****************************************************************************80
  inline void set_ncol(intT row, intT nc)
  {
#ifdef DEV_DEBUG
    this->CheckSizeIndex(row);
#endif

    index_[row + 1] = index_[row] + nc;
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
  inline void set_ncol(Array1D<int> const & ncol)
  {
#ifdef DEV_DEBUG
    //debug check - checking if size(Array1D) == nrow
    if(ncol.get_size(0) != size1_){
      std::cerr << "In List2D.h - "
                << "ERROR: Mismatch of num rows and num of col specification."
                << "\n Number of rows: "<< size1_
                << "\n Number of ncol specs: " << ncol.get_size(0)
                << std::endl;
      SystemModule::my_exit();
    }
#endif

    if(size1_ > 0) {//check if we have columns
      //setup index array based on ncol input
      index_[0] = 0;
      for(intT irow = 0; irow < size1_; irow++){
        //---> Index;
        index_[irow+1] = index_[irow] + ncol(irow);
      }

#ifdef DEV_DEBUG
      //debug check - checking if total num of entries <= entire size of List
      intT targetSize = index_[size1_] - index_[0];

      if(targetSize > size_){
        std::cerr << "In List2D.h - "
                  << "ERROR: Requesting too many elements."
                  << "\n Target total size of requested columns: " << targetSize
                  << "\n Preallocated total size: " << size_ << std::endl;
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
  inline intT get_ncol(intT row) const
  {

#ifdef DEV_DEBUG
    this->CheckSizeIndex(row);
#endif

    return( index_[row + 1] - index_[row] );
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
  inline const dataT& operator () (intT i, intT j) const
  {
#ifdef DEV_DEBUG
    this->CheckBounds(i,j);
#endif

    //---> The return of this operator is the ith reference of data pointer
    return Array1D<dataT>::operator()(index_[i] + j);
  }

//****************************************************************************80
//!
//! \brief () : Parenthesis operator for accessing array
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] i First index of the array you want to access with parenthesis
//! \param[in] j Second index of the array you want to access with parenthesis
//****************************************************************************80
  inline dataT& operator () (intT i, intT j)
  {
#ifdef DEV_DEBUG
    this->CheckBounds(i,j);
#endif

    //---> The return of this operator is the ith reference of data pointer
    return Array1D<dataT>::operator()(index_[i] + j);
  } // End

//****************************************************************************80
//!
//! \brief () : Parenthesis operator for acessing array, for const. correctness
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] i The index of the array you want to acess with parenthesis
//****************************************************************************80
  inline const dataT& operator () (intT i) const
  {

#ifdef DEV_DEBUG
    this->CheckSize(i);
#endif

    //---> The return of this operator is the ith reference of data pointer
    return Array1D<dataT>::operator()(i);
  } // End ()

//****************************************************************************80
//!
//! \brief () : Parenthesis operator for acessing array
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] i The index of the array you want to acess with parenthesis
//****************************************************************************80
  inline dataT& operator () (intT i)
  {
    //calls constant version of operator (i)
    return const_cast<dataT&>(static_cast<const List2D&>(*this).
                              operator()(i));
  } // End ()

//****************************************************************************80
//!
//! \brief get_lead_size : Get the size of leading index of list
//! \details
//! \nick
//! \version $Rev: 5 $
//****************************************************************************80
  inline intT get_lead_size( ) const
  {
    return((int)size1_);
  }// End get_lead_size

//****************************************************************************80
//!
//! \brief get_total_size : Get the total size of the list
//! \details
//! \nick
//! \version $Rev: 5 $
//****************************************************************************80
  inline intT get_total_size( ) const
  {
    return((int)size_);
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
    return(Array1D<dataT>::mem_ + mem_);
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
    for (intT i = 0; i < a.size1_; i++) {
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
  inline dataT const * get_ptr(const intT& i, const intT& j) const
  {
#ifdef DEV_DEBUG
    this->CheckBounds(i,j);
#endif
    return Array1D<dataT>::get_ptr(index_[i] + j);
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
  inline dataT* get_ptr(const intT& i, const intT& j)
  {

#   ifdef DEV_DEBUG
    this->CheckBounds(i,j);
#endif
    return Array1D<dataT>::get_ptr(index_[i] + j);
    
  }// End get_ptr
//****************************************************************************80
//!
//! \brief get_ptr : Allow access to pointer at specified single address
//! \details Constant correct version
//! \nick
//! \version
//! \param[in] i The address you want to access
//****************************************************************************80
  inline dataT const * get_ptr(const intT& i) const
  {
#ifdef DEV_DEBUG
    this->CheckSize(i);
#endif
    return Array1D<dataT>::get_ptr(i);
  }// End get_ptr

//****************************************************************************80
//!
//! \brief get_ptr : Allow access to pointer at specified single address
//! \details
//! \nick
//! \version $Rev: 5 $
//! \param[in] i The address you want to access
//****************************************************************************80
  inline dataT* get_ptr(const intT& i)
  {

#   ifdef DEV_DEBUG
    this->CheckSize(i);
#endif
    return Array1D<dataT>::get_ptr(i);

  }// End get_ptr
//****************************************************************************80
//! \brief get_index_ptr : Returns the pointer to the indexing array
//! \details
//! \nick
//! \version $Rev$
//! \param[in] i The leading address
//****************************************************************************80
  inline intT const * get_index_ptr(const intT& i ) const
  {
#ifdef DEV_DEBUG
    this->CheckSizeIndex(i);
#endif
    return(index_ + i);
  }

//****************************************************************************80
//! \brief get_index_ptr : Returns the pointer to the indexing array
//! \details
//! \nick
//! \version $Rev$
//! \param[in] i The leading address
//****************************************************************************80
  inline int* get_index_ptr(const intT& i )
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
  inline intT get_index(const intT& i, const intT& j) const
  {
#ifdef DEV_DEBUG
    this->CheckBounds(i,j);
#endif
    return( index_[i] + j );
  }

//****************************************************************************80
//! \brief RowBegin() : The index that begins a row; 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] i The row you want to begin
//! \return index[i] Returns the indexing array value for that row
//****************************************************************************80
  inline intT RowBegin(const intT& i) 
  {
      return index_[i];
  }// End RowBegin 

//****************************************************************************80
//! \brief RowEnd() : The index that ends a row; 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] i The row you want to end
//! \return index[i+1] Returns the indexing array value for that row
//****************************************************************************80
  inline intT RowEnd(const intT& i) 
  {
      return index_[i+1];
  }// End RowEnd
 
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
  std::string MemoryDiagnostic(std::string const & var_name) const {
   
     std::ostringstream stream;
     stream << var_name << "(" << this->get_lead_size() << ")(" 
	    << this->get_total_size() << "):\t "
	    << this->get_mem() << " MB" << std::endl;
     
     //---> However ostringstream is not copyable...so return the string
     return stream.str();

  }
private:
  //+++++++++++++++++++++++++++++++ PRIVATE STUFF ++++++++++++++++++++++++++++++
  intT size1_; /*!< Size of the leading dimension of the list */
  intT size_; /*!< Total size of the array */
  intT* index_; /*!< Index pointer that allows for addressing 2D */
  realT mem_;
//****************************************************************************80
//! \brief List2D : Copy constructor...private so it can't be called, which
//!                  blocks copy construction of this class.
//! \details A copy constructor is a constructor that takes a argument of the
//!      type defined by the class.  In this case List2D is the class.
//!      It then seeks to define operators for constructing a new instance
//!      using the old one.  We are purposely blocking this capability by
//!      defining this function as a private function with NO CODE.
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
    if(i >= size_){
      std::cerr << "ERROR: In List2D.h - "
                << "Attempting to access index "<< i
                << ".   Number of entries is " << size_ << std::endl;
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
    if(i >= size1_){
      std::cerr << "ERROR: In List2D.h - Over bounds on 1st index. "
		<< "Acessing List2D::("<< i <<","<< j <<"). Size is "
		<< "is " << size1_ << ","<< get_ncol(i) << std::endl;
      SystemModule::my_exit();
      
    }
    if( j >= get_ncol(i) ) {
      std::cerr << "ERROR: In List2D.h - Over bounds on 2nd index. "
		<< "Acessing List2D::("<< i <<","<< j <<"). Size is "
		<< "is " << size1_ << ","<< get_ncol(i) << std::endl;
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
    if(i >= size1_ ){
      std::cerr << "ERROR: In List2D.h - "
                << "Attempting to access index "<< i
                << ".   Size of indexing array: " << size1_ << std::endl;
      SystemModule::my_exit();
    }
  }// End CheckSizeIndex

}; // End class List2D
#endif
