// -*-c++-*-
#ifndef UNSTMESHGEOM_H
#define UNSTMESHGEOM_H
#include "my_incl.h"
#include "DataStructures/Array2D.h"
#include "SystemUtils/SystemModule.h"
#include "IO/UnstMeshReader.h"

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
  UnstMeshGeom(UnstMeshReader& mesh_reader);
  ~UnstMeshGeom();

  inline const intT& get_ndim() const {return ndim_;}
  inline const intT& get_nnode() const {return nnode_;}
  inline const Array2D<realT>& get_x() const {return x_;}
  
private:
  intT ndim_; //!< Number of physical dimensions
  intT nnode_;//!< Number of nodes
  
  Array2D<realT> x_; /*!< Coordinates, stores all the coordinates depending on
                number of dimensions.
                \verbatim
                1-D : (x)
                2-D : (x,y)
                3-D : (x,y,z)
                \endverbatim
              */


  
  
};// End UnstMeshGeom

#endif
