// -*-c++-*-
#ifndef UNSTMESHGEOM_H
#define UNSTMESHGEOM_H
#include "my_incl.h"
#include "DataStructures/Array2D.h"
#include "SystemUtils/SystemModule.h"
//****************************************************************************80
//! \brief A class for storing and maining the mesh geometry
//! \details 
//! \nick 
//! \version $Rev$ 
//! \date $Date$ 
//! 
//****************************************************************************80
class UnstMeshGeom 
{
public:
  
  //UnstMeshGeom(const intT& ndim, const intT& nnode, Array2D<realT>& x);
  UnstMeshGeom(std::string const & filename, std::string const & file_type);
  ~UnstMeshGeom();

  inline const intT& get_ndim() const {return ndim_;}
  inline const intT& get_nnode() const {return nnode_;}
  inline const Array2D<realT>& get_x() const {return x_;}
  
private:
  intT ndim_; //!< Number of physical dimensions
  intT nnode_;//!< Number of nodes

  Array2D<realT> x_;//!< The coordinates of each node in the vector
  
  
  
}// End UnstMeshGeom

#endif
