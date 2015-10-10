// -*-c++-*-
#ifndef HEAP_H
#define HEAP_H
#include "my_incl.h"
#include "Array1D.h"
#include <cmath>

//****************************************************************************80
//! \class Heap
//! \brief Class contains data for building a heap
//! \nick 
//! \version $Rev: 5 $
//! \tparam intT Template argument meant to mimic integer
//! \tparam realT Template argument meant to mimic real numbers
//****************************************************************************80

template <class dataT >  
class Heap
{
private:
  const intT hsize; //!< Heap size
  double mem; //!< Memory for the heap class
  Array1D<dataT> metric; //!< Metric array 
  Array1D<intT> front; //!< From array defining the thing to be heaped
  Array1D<intT> location;//!< Location tracker of stuff in the heap
//****************************************************************************80
//!
//! \brief Heap : Default constructor
//! \details
//! \nick 
//! \version $Rev$
//****************************************************************************80
  Heap()=delete;
  
public:
  
//****************************************************************************80
//!
//! \brief Heap : Constructor that takes the size of the heap
//! \details
//! \nick 
//! \version $Rev$
//****************************************************************************80
  Heap(const int& n) : hsize(n){
    metric.initialize(hsize);
    front.initialize(hsize);
    location.initialize(hsize);
       
    for(intT i = 0; i < hsize; i++){
      front(i) = i;
      location(i) = i;
    }
  }// End Heap

//****************************************************************************80
//!
//! \brief set_metric : Set value of metric for specified index
//! \details
//! \nick 
//! \version $Rev$
//! \param[in] i Index value
//! \param[in] val Value to set
//****************************************************************************80
  void set_metric(const intT& i, const dataT& val)
  {
    metric(i) = val;
  }// End set_metric

//****************************************************************************80
//!
//! \brief get_metric : get value of metric for specified index
//! \details
//! \nick 
//! \version $Rev$
//! \param[in] i Index value
//****************************************************************************80
  dataT get_metric(const intT& i)
  {
    return(metric(i));
  }// End get_metric

//****************************************************************************80
//!
//! \brief accum_metric : Accum value of metric for specified index
//! \details
//! \nick 
//! \version $Rev$
//! \param[in] i Index value
//! \param[in] val Value to set
//****************************************************************************80
  void accum_metric(const intT& i, const dataT& val)
  {
    metric(i) += val;
  }// End set_metric

//****************************************************************************80
//!
//! \brief sift_down : Takes a node and sifts it down to it's spot in the heap
//! \details
//! \nick 
//! \version $Rev$
//! \param[in] list_size The size of the list to work with != hsize 
//! \param[in] in_node Index to be inserted
//****************************************************************************80
  void sift_down(const intT& list_size, const intT& in_node)
  {
    int node = in_node;
    intT w = 0;
    
    while( w < list_size - 1) { // list traverse
      w = 2*node + 1;
      
      if( w < list_size - 1){
        //---> Get the reference length of the face at position node
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        intT j = front(w);
        dataT len1 = metric(j);
        //---> Face2 being the front face w + 1
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
        dataT comp_len = len1;
        
        j = front(w + 1);
        dataT len2 = metric(j);
        
        //---> If the metric of the w + 1 face if larger then
        //     use this as the comparison length
        if( len2 < len1) {
          w += 1;
          comp_len = len2;
        }
        
        j = front(node);
        dataT ref_len = metric(j);
        
        if( comp_len < ref_len){
          intT temp = front(node);
          front(node) = front(w);
          location(front(w)) = node;
          location(j) = w;
          front(w) = temp;
          
        }
            
        node = w;
        
      }
      else if(w == list_size - 1) {
        //---> Get the reference length of the j at position node
        intT j = front(w);
        dataT len1 = metric(j);
        
        dataT comp_len = len1;
        
        j = front(node);
        dataT ref_len = metric(j);
        
        if( comp_len < ref_len) {
          intT temp = front(node);
          front(node) = front(w);
          location(front(w)) = node;
          location(temp) = w;
          front(w) = temp;
            
        }
        node = w;
        
      } // end Check
                                           
                                           
    }// End list traverse
      

  }// End sift_down

//****************************************************************************80
//!
//! \brief sift_up : Takes a node and sifts it up to it's spot in the heap
//! \details
//! \nick 
//! \version $Rev$
//! \param[in] list_size The size of the list to work with != hsize 
//! \param[in] in_node Index to be inserted
//****************************************************************************80
  void sift_up(const intT& list_size, const intT& in_node)
  {
    intT node = in_node;
    intT parent = 0;
    intT cont = 1;
    if( list_size > 0 ) {
    }
    while( (node - 1) > -1  && cont == 1) {// traverse list
      
      parent = (intT)floor((node - 1 - 1)/2) + 1;
      
      //---> Look at the child's metric
      intT j = front(node);
      dataT ref_len = metric(j);
      
      //---> Now look at the parent's metric
      j = front(parent);
      dataT comp_len = metric(j);
       
      if( ref_len < comp_len) {
        intT temp = front(node);
        front(node) = front(parent);
        location(temp) = parent;
        location(front(parent)) = node;
        front(parent) = temp;
        
        node = parent;
      }
      else {
        //---> if the child in question in the not greater then it's parent 
        //     then it has found it's spot in the queue and we're done
        cont = 0;
      }

    }// End traverse list

  }// End sift_up


//****************************************************************************80
//!
//! \brief heapify : Makes a heap of the class
//! \details
//! \nick 
//! \version $Rev$
//! \param[in] list_size The size of the list to work with != hsize 
//****************************************************************************80
  void heapify(const intT& list_size) {
    //---> Find halfway point
    intT nhalf = (intT)round((realT)(list_size/2));
    
    for(intT i = (list_size - 1) - nhalf; i > -1; i--){
      //---> Call siftdown on all nodes below half way point
      sift_down(list_size, i);
    }
    
  }
//****************************************************************************80
//!
//! \brief pop_off_top : Gets the data element off the top of the heap
//! \details
//! \nick 
//! \version $Rev$
//! \param[in] list_size The size of the list to work with != hsize 
//****************************************************************************80
  dataT pop_off_top(const intT& list_size) {
    //---> Grab the top value
    dataT top = front(0);
    
    //---> Put the last element on top and sift down 
    front(0) = front(list_size - 1);
    location(front(list_size - 1)) = 0;
    sift_down(list_size, 0);
    front(list_size - 1) = 0;
    
    //---> Return to user
    return(top);
  }// End pop_off_top

//****************************************************************************80
//!
//! \brief remove_node : Removes specified node from the heap
//! \details
//! \nick 
//! \version $Rev$
//! \param[in] list_size The size of the list to work with != hsize
//****************************************************************************80
  void remove_node(const intT& list_size, const intT& node)
  {
    //---> Get location and put the bottom element on the heap into that
    // location.  
    intT loc = location(node);
    front(loc) = front(node);
    location(front(list_size-1)) = loc;
    //---> Now sift down to re-establish the heap
    sift_down(list_size, loc);
    //---> Last element in heap in -1 to mark deletion...we don't resize
    front(list_size - 1) = -1;
       
  } //End remove_node

//****************************************************************************80
//!
//! \brief add_node : Adds specified node to the heap
//! \details
//! \nick 
//! \version $Rev$
//! \param[in] list_size The size of the list to work with != hsize
//****************************************************************************80
  void add_node(const intT& list_size, const intT& node)
  {
    front(list_size) = node;
    location(node) = list_size;
    sift_up(list_size + 1, list_size); 

  }// End add_node

}; // End class Heap

#endif
