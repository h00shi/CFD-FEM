//-*-c++-*-
#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H
#include "my_incl.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/List2D.h"

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
template <class dataT>
class SparseMatrix {

protected:
  intT nnode_ ; //!< Number of graph nodes that sparse matrix is built from
  intT nedge_;//!< Number of graph edges
  intT nrow_; //!< Number of rows in the matrix 
  intT ncol_; //!< Number of columns in the matrix
  intT nnz_;
  
  realT mem_; //!< Number of 
  Array1D<dataT> data_; //!< Data array 
  Array1D<intT> self_adj_index_; /*!< For a node i self_adj_index_(i) = j :
				   adjacency(i,j) = i */
  Array2D<intT> edge_adj_index_; /*!< For an edge e edge_adj_index(e,0) = j :
				   adjacency(edge2node(e,0),j) = edge2node(e,1)
				 */
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
	       const Array1D<intT>& nrow_per_node) : 
    adj_(adjacency), edge2node_(edge2node), nrow_per_node_(nrow_per_node)
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
    self_adj_index_.initialize(nnode_);
    edge_adj_index_.initialize(nedge_,2);
    mem_ = data_.get_mem() + self_adj_index_.get_mem() + 
      edge_adj_index_.get_mem();
    
    //---> Setup Self adjacency index
    for(intT n = 0; n < nnode_; n++){// Node Loop
      for(intT j = 0; j < adjacency.get_ncol(n); j++) { // Neighbor Loop 
      	intT node = adjacency(n,j);
	if(node == n){self_adj_index_(n) = j;break;}
      }// End Neighbor loop
    }
    
    //---> Setup edge adjacency index
    for(intT e = 0; e < nedge_; e++){ // Edge loop
      //---> Left and right nodes
      intT nl = edge2node(e,0);
      intT nr = edge2node(e,1);
      
      for(intT j = 0; j < adjacency.get_ncol(nl); j++){// left node neighbor idx
	intT node = adjacency(nl,j);
	if(node == nr){edge_adj_index_(e,0) = j; break;}
      }// End left node neighbor idx
      
      for(intT j = 0; j < adjacency.get_ncol(nr); j++){//right node neighbor idx
	intT node = adjacency(nr,j);
	if(node == nl){edge_adj_index_(e,1) = j; break;}
      }//End right node neighbor idx
      
    }// End Edge loop 
  
  } //End SparseMatrix

public:

//****************************************************************************80
//! \brief get_nrow_ : Gets the number of rows of the matrix 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  inline const intT& get_nrow() const { return nrow_;}

//****************************************************************************80
//! \brief get_ncol_ : Gets the number of columns of the matrix 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  inline const intT& get_ncol() const { return ncol_;}

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
//! \brief get_data : Gets the reference to the data Array 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  inline const Array1D<dataT>& get_data() const { return data_;}

//****************************************************************************80
//! \brief get_self_adj_index : Gets the reference to the self_adj_index array 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  inline const Array1D<intT>& get_self_adj_index() const 
  {
    return self_adj_index_;
  } // End get_self_adj_index
//****************************************************************************80
//! \brief get_edge_adj_index : Gets the reference to the edge_adj_index
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  inline const Array2D<intT>& get_edge_adj_index() const 
  {
    return edge_adj_index_;
  } // end get_edge_adj_index
//****************************************************************************80
//! \brief 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  virtual ~SparseMatrix() {  } 

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
    out_stream << "Number of Graph Nodes: " << nnode_ << std::endl;
    out_stream << "Number of Graph Edges: " << nedge_ << std::endl;
    out_stream << "Number of Rows: " << nrow_ << std::endl;
    out_stream << "Number of Columns: " << ncol_ << std::endl;
    out_stream << std::endl;
    out_stream << data_.MemoryDiagnostic("data_");
    out_stream << self_adj_index_.MemoryDiagnostic("self_adj_index_");
    out_stream << edge_adj_index_.MemoryDiagnostic("edge_adj_index_");
  }// End Diagnostic 
private:
  
  //---> References to adjacency
  const List2D<intT>& adj_;         //!< Adjacency of the graph
  const Array2D<intT>& edge2node_;  //!< Edge ordering
  const Array1D<intT>& nrow_per_node_;  //!< Number of rows per graph node
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
