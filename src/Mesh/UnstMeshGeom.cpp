//****************************************************************************80
//UnstMeshGeom::UnstMeshGeom(Array2D<realT>&& x) : x_(x);
UnstMeshGeom::UnstMeshGeom(std::string const & filename, 
                           std::string const & file_type)
{
  //MeshIO::NKBGridIO::ReadNodesFromFile(filename, file_type);

  //---> Query coordinate array for number of nodes 
  nnode_ = x_.get_size(0);
  //---> Query coordinate array for number of dimensions
  ndim_  = x_.get_size(1);
}// End UnstMeshGeom::UnstMeshGeom
