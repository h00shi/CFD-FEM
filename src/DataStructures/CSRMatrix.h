//-*-c++-*-
#ifndef CSRMATRIX_H
#define CSRMATRIX_H
#include "my_incl.h"
#include "SparseMatrix.h"

//****************************************************************************80
//! \brief This is the base class for csr formatted sparse matricies.  
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
template<class dataT> 
class CSRMatrix : public SparseMatrix<dataT, CSRMatrix<dataT> >
{
protected:
  Array1D<intT> row_offset_; //!< Row offset pointer as specified by CSR
  List2D<intT> column_idx_; //!< Column offset point as specified by CSR
  Array1D<intT> ncol_per_row_; //!< Number of columns per Row
  List2D<intT> adj_data_offset_; //!< Adjacency based data offset 

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
    : SparseMatrix< dataT, CSRMatrix<dataT> >::SparseMatrix(adjacency, 
                                                            edge2node, 
                                                            nrow_per_node)
  {
    //---> Some temporary storage...makes for shorter code;
    intT nrow  = SparseMatrix< dataT, CSRMatrix<dataT> >::nrow_;
    intT nnz   = SparseMatrix< dataT, CSRMatrix<dataT> >::nnz_;
    intT nnode = SparseMatrix< dataT, CSRMatrix<dataT> >::nnode_;

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

    //---> Setup the diagonal and offset_ pointers


  } // End CSRMatrix
//****************************************************************************80
//! \brief operator() : A parenthetical operator to access specified spot in 
//!        matrix based on graph node, neighbor index and block (i,j)
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  inline dataT& operator()(const intT& node, const intT& j_neighbor, 
			   const intT& block_row, const intT& block_col)
  {
 
    return SparseMatrix< dataT, CSRMatrix<dataT> >::data_(node);
  }// End operator()

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
//! \brief  Diagnostic : Returns diagnostic information to specified stream
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  inline void Diagnostic(std::ostream& out_stream)
  {
    out_stream << "Row Offset: " << std::endl << row_offset_ << std::endl;
    out_stream << "Column Idx: " << std::endl << column_idx_ << std::endl;
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
