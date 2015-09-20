//-*-c++-*-
#ifndef CSRMATRIX_H
#define CSRMATRIX_H
#include "my_incl.h"
#include "SparseMatrix/SparseMatrix.h"

//****************************************************************************80
//! \brief This is the base class for csr formatted sparse matricies.  
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
template<class dataT> 
class CSRMatrix : public SparseMatrix<dataT>
{
protected:
  Array1D<intT> row_offset_; //!< Row offset pointer as specified by CSR
  List2D<intT>  column_idx_; //!< Column offset point as specified by CSR
  List2D<intT>  adj_data_offset_; //!< Pointer from adj to data indexing 

public:

//****************************************************************************80
//! \brief CSRMatrix : Constructor for CSR Matrix from adjacency and edge defs
//! \details  Deleted
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] adjacency The adjacency list of graph from which to create
//! \param[in] edge2node The edg2node of graph 
//! \param[in] nrow_per_node The number of rows per node
//****************************************************************************80
  CSRMatrix(const List2D<intT>& adjacency, const Array2D<intT>& edge2node, 
            const Array1D<intT>& nrow_per_node) 
    : SparseMatrix<dataT>::SparseMatrix(adjacency, edge2node, nrow_per_node)
  {
    //---> Some temporary storage...makes for shorter code;
    intT nrow  = SparseMatrix<dataT>::nrow_;
    intT nnz   = SparseMatrix<dataT>::nnz_;
    intT nnode = SparseMatrix<dataT>::nnode_;
   
    //---> Initialize indexing arrays for CSR
    row_offset_.initialize(nrow + 1);
    column_idx_.initialize(nrow, nnz);
    adj_data_offset_.intialize_copy_pattern(adjacency);
    
    row_offset_(0) = 0;
    intT row = 0;
    for(intT n = 0; n < nnode; n++){// Node loop 
      //---> Count up number of columns
      intT ncol = 0;
      for(intT j = 0; j < adjacency.get_ncol(n); j++){ // Neighbor loop 
	ncol += nrow_per_node(adjacency(n,j));
      }// End Neighbor loop 
      
      //---> Now form offset for rows irow to irow + nrow_per_node(n);
      for(intT v = 0; v < nrow_per_node(n); v++){ // Row var
        //---> Set row_offset
        row_offset_(row + 1) = row_offset_(row) + ncol;
        //---> Set the number of non-zero columns for this row
        column_idx_.set_ncol(row, ncol);
               
        //---> Form column_idx_
        intT col = 0;
        intT node_start = 0;
        intT icol = 0;

        for(intT j = 0; j < adjacency.get_ncol(n); j++){// Neighbor loop
          //---> Find all zero columns in this row up to first neighbor
          intT node_end = adjacency(n,j) - 1;
	 
          //---> Loop over node indices that are zeros in this row
          for(intT k = node_start; k <= node_end; k++){
            col += nrow_per_node(k);
          }
	  //---> Loop over the number of variables for the adjacent node
          for(intT k = 0; k < nrow_per_node(adjacency(n,j)); k++){
            column_idx_(row,icol) = col;
            col++;
            icol++;
          }
	  // std::cout << row << " " << j << " " << adjacency(n,j) << " " 
	  // 	    << node_start << " " << node_end 
	  // 	    << " " << col << " " << icol << std::endl;
	  // SystemModule::pause();
          node_start = adjacency(n,j) + 1;
          
        } // neighbor_loop 
	//---> increment row counter
        row++;

      } // End Row var
     
    } // End node loop 
    
    //---> Setup the adj_data_offset_ and diag pointers
    
    intT index = 0;
    for(intT n = 0; n < nnode; n++){// Node Loop 
      //---> Cycle through neighbors of this node
      for (intT j = 0; j < adjacency.get_ncol(n); j++) {
	adj_data_offset_(n,j) = index;
	index += nrow_per_node(adjacency(n,j));
      }
      /*---> Account for all rows between first row of node n and 
	first row of node n+1.  */
      for (intT j = 0; j < adjacency.get_ncol(n); j++) {
	index += (nrow_per_node(n)-1)*nrow_per_node(adjacency(n,j));
      }
      
    }// End Node Loop 

    SparseMatrix<dataT>::mem_ +=
      row_offset_.get_mem() + column_idx_.get_mem();
      
  } // End CSRMatrix
//****************************************************************************80
//! \brief operator() : A parenthetical operator to access specified spot in 
//!        matrix based on graph node, neighbor index and block (i,j)
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] node The graph node where you want to acess matrix
//! \param[in] j_neighbor The local neighbor index of neighbor 
//!            you want to access
//! \param[in] block_row Row of block you want
//! \param[in] block_col Column of block you want
//****************************************************************************80
  inline dataT& operator()(const intT& node, const intT& j_neighbor, 
			   const intT& block_row, const intT& block_col)
  {
   
    intT i = adj_data_offset_(node, j_neighbor) + 
      (row_offset_(node + block_row + 1) - 
       row_offset_(node + block_row))*block_row + 
      block_col;
    
    return SparseMatrix<dataT>::data_(i);
  }// End operator()

//****************************************************************************80
//! \brief operator() : A parenthetical operator to access specified spot in 
//!        matrix based on graph node, neighbor index and block (i,j)
//!        CONST VERSION
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] node The graph node where you want to acess matrix
//! \param[in] j_neighbor The local neighbor index of neighbor 
//!            you want to access
//! \param[in] block_row Row of block you want
//! \param[in] block_col Column of block you want
//****************************************************************************80
  inline const dataT& operator() (const intT& node,
				  const intT& j_neighbor, 
				  const intT& block_row, 
				  const intT& block_col) const
  {
   
    intT i = adj_data_offset_(node, j_neighbor) + 
      (row_offset_(node + block_row + 1) - 
       row_offset_(node + block_row))*block_row + 
      block_col;
    
    return SparseMatrix<dataT>::data_(i);
  }// End operator()

//****************************************************************************80
//! \brief Diagonal : Returns reference to the diagonal element  
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] node The graph node where you want to acess matrix
//! \param[in] block_row Row of block you want
//! \param[in] block_col Column of block you want
//****************************************************************************80
  inline dataT& Diagonal(const intT& node, 
			 const intT& block_row, 
			 const intT& block_col)
  {
    intT jneighbor = SparseMatrix<dataT>::self_adj_index_(node);
    return this->operator()(node, jneighbor, block_row, block_col);
       
  }// End Diagonal

//****************************************************************************80
//! \brief Diagonal : Returns reference to the diagonal element. CONST VERSION
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] node The graph node where you want to acess matrix
//! \param[in] block_row Row of block you want
//! \param[in] block_col Column of block you want
//****************************************************************************80
  inline const dataT& Diagonal(const intT& node, 
			       const intT& block_row, 
			       const intT& block_col) const
  {
    intT jneighbor = SparseMatrix<dataT>::self_adj_index_(node);
    return this->operator()(node, jneighbor, block_row, block_col);
  }// End Diagonal
//****************************************************************************80
//! \brief OffDiagonal : Returns reference to the diagonal element. 
//!      
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] edge The graph edge you want the off-diagonal of
//! \param[in] side Side indicates if you want the left or right 0=left, 
//!            1 = right
//! \param[in] block_row Row of block you want
//! \param[in] block_col Column of block you want
//****************************************************************************80
  inline dataT& OffDiagonal(const intT& edge, const intT& side, 
			    const intT& block_row, 
			    const intT& block_col)
  {
    intT node = SparseMatrix<dataT>::edge2node_(edge,side);
    intT jneighbor = SparseMatrix<dataT>::edge_adj_index(edge,side);
    return this->operator()(node, jneighbor, block_row, block_col);
    
  }// End OffDiagonal

//****************************************************************************80
//! \brief OffDiagonal : Returns reference to the diagonal element. 
//!        CONST VERSION
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! \param[in] edge The graph edge you want the off-diagonal of
//! \param[in] side Side indicates if you want the left or right 0=left, 
//!            1 = right
//! \param[in] block_row Row of block you want
//! \param[in] block_col Column of block you want
//****************************************************************************80
 inline const dataT& OffDiagonal(const intT& edge, const intT& side, 
			    const intT& block_row, 
			    const intT& block_col) const
  {
    intT node = SparseMatrix<dataT>::edge2node_(edge,side);
    intT jneighbor = SparseMatrix<dataT>::edge_adj_index(edge,side);
    return this->operator()(node, jneighbor, block_row, block_col);
  }// End OffDiagonal

//****************************************************************************80
//! \brief get_row_offset() : Returns the row_offset array;
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  inline const Array1D<intT>&  get_row_offset() const {return row_offset_;}

//****************************************************************************80
//! \brief get_column_idx() : Returns the column idx array
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  inline const List2D<intT>& get_column_idx() const {return column_idx_;}

//****************************************************************************80
//! \brief get_adj_data_offset : Returns the adj_data_offset array
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  inline const List2D<intT>& get_adj_data_offset() const 
  {
    return adj_data_offset_;
  } // End get_adj_data_offset;

//****************************************************************************80
//! \brief  Diagnostic : Returns diagnostic information to specified stream
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  inline void Diagnostic(std::ostream& out_stream)
  {
    out_stream << std::endl;
    out_stream << "--------------------- CSRMatrix Diagnostics ----------------"
            << std::endl;
        
    SparseMatrix<dataT>::Diagnostic(std::cout);
    out_stream << row_offset_.MemoryDiagnostic("row_offset_"); 
    out_stream << column_idx_.MemoryDiagnostic("column_idx_");
    out_stream << adj_data_offset_.MemoryDiagnostic("adj_data_offset_");
    //---> Total Memory:
    out_stream << "CSR Matrix Memory: " 
	       << SparseMatrix<dataT>::mem_ << " MB" 
	       << std::endl;
    out_stream << std::endl;
    out_stream << "------------------- End CSRMatrix Diagnostics --------------"
            << std::endl;
  }// End Diagnostic

private:
//****************************************************************************80
//! \brief CSRMatrix : Default constructor for CSR Matrix
//! \details  Deleted
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  CSRMatrix() = delete;

};
#endif
