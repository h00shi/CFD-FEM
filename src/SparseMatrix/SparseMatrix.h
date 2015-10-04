//-*-c++-*-
#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H
#include "my_incl.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/List2D.h"
#include "DataStructures/Graph.h"

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
  const Array1D<intT>& nrow_per_node_;  //!< Number of rows per graph node
  const Graph& graph_;//!< Graph form which we build the csr matrix
  intT nnode_ ; //!< Number of graph nodes that sparse matrix is built from
  intT nedge_;//!< Number of graph edges
  intT nrow_; //!< Number of rows in the matrix 
  intT ncol_; //!< Number of columns in the matrix
  intT nnz_;
  
  realT mem_; //!< Number of 
  Array1D<dataT> data_; //!< Data array 
   
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
  SparseMatrix(const Graph& graph, const Array1D<intT>& nrow_per_node) : 
    graph_(graph), nrow_per_node_(nrow_per_node)
  {
    //---> Reference variables
    const List2D<intT>& adjacency = graph_.get_GraphAdj();
    const Array2D<intT>& edge2node = graph_.get_GraphEdge2Node();
    
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
	intT neighbor_node = adjacency(n,j);
	nnz_ += nrow_per_node(n)*nrow_per_node(neighbor_node);

      } // End Neighbor Loop 
    } // End Node loop 
    
    //---> Initialize variables
    data_.initialize(nnz_);
    mem_ = data_.get_mem();
  
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
//! \brief get_nnz_ : Gets the number of non-zeros of the matrix 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  inline const intT& get_nnz() const { return nnz_;}
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
//! \brief get_data : Gets the reference to the data array 
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
  inline const Array1D<dataT>& get_data() const { return data_;}

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
  }// End Diagnostic 
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
