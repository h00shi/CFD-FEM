//-*-c++-*-
#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H
#include "my_incl.h"
#include "Array1D.h"
#include "Array2D.h"
#include "List2D.h"

//****************************************************************************80
//! \brief This is the base class for sparse matricies.  
//! \details This class will use the curiously recurring template pattern to
//!           generate base classes that make their own bases as necessary. 
//!           Data structure-wise this class and it's children make 
//!           sparse matricies from specified graphs.  These graphs are 
//!           made up of nodes and edges.  The matrix is built of blocks that 
//!           vary in size >= 1 at each node.   
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
template <class dataT, template <class = dataT> class DerivedMatrixType >
class SparseMatrix {

protected:
  intT nnode_ ; //!< Number of graph nodes that sparse matrix is built from
  intT nedge_;//!< Number of graph edges
  intT nrow_; //!< Number of rows in the matrix 
  intT ncol_; //!< Number of columns in the matrix
  intT nnz_;
  
  realT mem_; //!< Number of 
  Array1D<dataT> data_; //!< Data array 
  Array1D<intT> node2diag_; //!< Given node point to diagonal for that node
  Array2D<intT> off_diag_index_; //!< Given a graph edge point to off-diagonal
 
//****************************************************************************80
//! \brief SparseMatrix : Constructor of SpareMatrix From Graph
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] adjacency The adjacency list of graph from which to create
//! \param[in] edge2node The edg2node of graph 
//! \param[in] nrow_per_node The number of rows per node
//****************************************************************************80
  SparseMatrix(const List2D<intT>& adjacency, const Array2D<intT>& edge2node, 
	       const Array1D<intT>& nrow_per_node)
  {
    //---> Get the graph info
    nnode_ = adjacency.get_lead_size();
    nedge_ = edge2node.get_size(0);
    
    //---> Compute the number of rows/columns
    nrow_ = 0;
    for(intT n = 0; n < nnode_; n++){
      nrow_+= nrow_per_node(n);
    }
    ncol_ = nrow_;
    
    //---> Get total number of non-zero entries
    nnz_ = 0;
    for(intT n = 0; n < nnode_; n++){ // Node loop 
      for(intT j = 0; j < adjacency.get_ncol(n); j++) { // Neighbor Loop 
	intT node = adjacency(n,j);
	nnz_ += nrow_per_node(n)*nrow_per_node(node);
      } // End Neighbor Loop 
    } // End Node loop 
    
    //---> Initialize variables
    data_.initialize(nnz_);
    node2diag_.initialize(nnode_);
    off_diag_index_.initialize(nedge_, 2);
    

  } //End SparseMatrix

public:

//****************************************************************************80
//! \brief operator() : Parenthetical operator that accesses based on graph node
//!                     and neighbor index.
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] node The graph node
//! \param[in] j_neighbor The jth neighbor of node i, j = 0,1,2,3...N_neighbors
//! \param[in] block_row The row of the block at node (node, j_neighbor)
//! \param[in] block_col The column of the block at node (node, j_neighbor)
//****************************************************************************80
  inline dataT& operator()(const intT& node, const intT& j_neighbor, 
			    const intT& block_row, const intT& block_col)
  {
    return static_cast< DerivedMatrixType<dataT>* >(this)->
      operator()(node, j_neighbor, block_row, block_col);
  }// End operator ()


//****************************************************************************80
//! \brief Diagonal : Returns reference to the diagonal element  
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  inline dataT& Diagonal(const intT& node, const intT& block_row, 
			 const intT& block_col)
  {
    
    return static_cast< DerivedMatrixType<dataT>* >(this)-> 
      Diagonal(node, block_row, block_col);
    
  }// End Diagonal

//****************************************************************************80
//! \brief OffDiagonal : Given an node and edge in the graph fill off-diagonal 
//!                      entries
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  inline dataT& OffDiagonal(const intT& node, const intT& edge, 
			    const intT& side, const intT& block_row,
			    const intT& block_col)
  {
    return static_cast< DerivedMatrixType<dataT>* >(this)-> 
      OffDiagonal(node, edge, side, block_row, block_col);

  }// End OffDiagonal

//****************************************************************************80
//! \brief get_nrow_ : Gets the number of rows of the matrix 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  inline intT get_nrow() const { return nrow_;}

//****************************************************************************80
//! \brief get_ncol_ : Gets the number of columns of the matrix 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  inline intT get_ncol() const { return ncol_;}

//****************************************************************************80
//! \brief get_data : Gets the reference to the data Array 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  inline Array1D<dataT>& get_data() { return data_;}

//****************************************************************************80
//! \brief 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  virtual ~SparseMatrix() {  } 

private:
//****************************************************************************80
//! \brief SparseMatrix : Default constructor...deleted so you can't call it 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  SparseMatrix() = delete;
};
#endif
