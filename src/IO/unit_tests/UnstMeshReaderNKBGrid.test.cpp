#include "IO/UnstMeshReaderNKBGrid.h"
#include "gtest/gtest.h"
TEST(UnstMeshReaderNKBGrid, ReadIn)
{
  UnstMeshReaderNKBGrid grid_reader;

  grid_reader.Open("Square.grid");
  Array2D<realT> x;
  x = grid_reader.ReadNodes();
 
}
