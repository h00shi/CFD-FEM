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
  Array1D<intT> ncol_per_row_;//!< Number of columns per Row
  

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
    //---> Some temp storage...makes for shorter code;
    intT nrow  = SparseMatrix< dataT, CSRMatrix<dataT> >::nrow_;
    intT nnz   = SparseMatrix< dataT, CSRMatrix<dataT> > ::nnz_;
    intT nnode = SparseMatrix< dataT, CSRMatrix<dataT> > :: nnode_;

    //---> Initialize indexing arrays for CSR
    row_offset_.initialize(nrow);
    column_idx_.initialize(nrow, nnz);
    
    row_offset_(0) = 0;
    intT row = 0;
    for(intT n = 0; n < nnode; n++){// Node loop 
      intT ncol = 0;
      for(intT j = 0; j < adjacency.get_ncol(n); j++){ // Neighbor loop 
        intT node = adjacency(n,j);
        ncol += nrow_per_node(node);
      }// End Neighbor loop 
      
      //---> Now form offset for rows irow to irow + nrow_per_node(n);
      for(intT v = 0; v < nrow_per_node(n); v++){ // Row var
        //---> Set row_offset
        row_offset_(row + 1) = row_offset_(row) + ncol;
        //---> Set the number of non-zero columns for this row
        column_idx_.set_ncol(row, ncol);
        //---> increment row counter
        row++;
        
        //---> Form column_idx_
        intT col = 0;
        intT node_start = 0;
        intT icol = 0;
        for(intT j = 0; j < adjacency.get_ncol(n); j++){// Neighbor loop
          //---> Find all zero columns in this row upto first neighbor
          intT node_end = adjacency(n,j);
          //---> Loop over node indicies that are zeros in this row
          for(intT k = node_start; k < node_end; k++){
            col += nrow_per_node(k);
          }
          for(intT k = 0; k < nrow_per_node(node_end); k++){
            column_idx_(row,icol) = col;
            col++;
            icol++;
          }
          node_start = node_end;
          
        } // neighbor_loop 
      } // End Row var
      
    } // End node loop 

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
