/*
 * DataStructureWriter.test.cpp
 *
 *  Created on: Oct 31, 2015
 *      Author: rabbit
 */

#include "gtest/gtest.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/List2D.h"
#include "IO/DataStructureWriter.h"
#include "IO/DataStructureReader.h"

TEST(DataStructureReadWrite, Array1DBinary)
{
  std::stringstream test(std::ios::in | std::ios::out | std::ios::binary);

  // StreamParser parser;
   DataStructureWriter writer;
   DataStructureReader reader;
   Array1D<realT> x(3);
   x(0) = 1.1;
   x(1) = 2.2;
   x(2) = 3.3;

   writer.WriteArray1D(test, DataStructureIO::IOMode::BINARY, x);
   Array1D<realT> y = reader.ReadArray1D<realT>(test,
                                          DataStructureIO::IOMode::BINARY);

   for(intT i = 0; i < 3; i++){
     EXPECT_DOUBLE_EQ(x(i), y(i));
   }

}
TEST(DataStructureIO, Array2DBinary) {
  std::stringstream test(std::ios::in | std::ios::out | std::ios::binary);
  DataStructureReader reader;
  DataStructureWriter writer;
  Array2D<realT> x(2,3);
  x(0,0) = 1.1;
  x(0,1) = 2.2;
  x(0,2) = 3.3;
  x(1,0) = 4.4;
  x(1,1) = 5.5;
  x(1,2) = 6.6;

  //write out the array out
  writer.WriteArray2D(test, DataStructureIO::IOMode::BINARY, x);

  //read array back in from stream
  Array2D<realT> y = reader.ReadArray2D<realT>
    (test, DataStructureIO::IOMode::BINARY);

  for(intT i = 0; i < 2; i++){
    for(intT j = 0; j < 3; j++){
     EXPECT_DOUBLE_EQ(x(i,j), y(i,j));
    }
  }
}
TEST(DataStructureIO, List2DBinary) {
  std::stringstream test(std::ios::in | std::ios::out | std::ios::binary);
  DataStructureReader reader;
  DataStructureWriter writer;
  List2D<realT> x(2,5);
  Array1D<intT> col(2);
  col(0) = 2;
  col(1) = 3;
  x.set_ncol(col);

  x(0,0) = 1.1;
  x(0,1) = 2.2;
  x(1,0) = 3.3;
  x(1,1) = 4.4;
  x(1,2) = 5.5;

  //write out the array out
  writer.WriteList2D(test, DataStructureIO::IOMode::BINARY, x);

  //read array back in from stream
  List2D<realT> y = reader.ReadList2D<realT>
    (test, DataStructureIO::IOMode::BINARY);

  for(intT i = 0; i < 2; i++){
    for(intT j = 0; j < x.get_ncol(i); j++){
     EXPECT_DOUBLE_EQ(x(i,j), y(i,j));
    }
  }
}



TEST(DataStructureReadWrite, Array1DASCII)
{
  std::stringstream test(std::ios::in | std::ios::out );

  // StreamParser parser;
   DataStructureWriter writer;
   DataStructureReader reader;
   Array1D<realT> x(3);
   x(0) = 1.1;
   x(1) = 2.2;
   x(2) = 3.3;

   writer.WriteArray1D(test, DataStructureIO::IOMode::ASCII, x);
   Array1D<realT> y = reader.ReadArray1D<realT>(test,
                                          DataStructureIO::IOMode::ASCII);

   for(intT i = 0; i < 3; i++){
     EXPECT_DOUBLE_EQ(x(i), y(i));
   }

}
TEST(DataStructureIO, Array2DASCII) {
  std::stringstream test(std::ios::in | std::ios::out );
  DataStructureReader reader;
  DataStructureWriter writer;
  Array2D<realT> x(2,3);
  x(0,0) = 1.1;
  x(0,1) = 2.2;
  x(0,2) = 3.3;
  x(1,0) = 4.4;
  x(1,1) = 5.5;
  x(1,2) = 6.6;

  //write out the array out
  writer.WriteArray2D(test, DataStructureIO::IOMode::ASCII, x);

  //read array back in from stream
  Array2D<realT> y = reader.ReadArray2D<realT>
    (test, DataStructureIO::IOMode::ASCII);

  for(intT i = 0; i < 2; i++){
    for(intT j = 0; j < 3; j++){
     EXPECT_DOUBLE_EQ(x(i,j), y(i,j));
    }
  }
}
TEST(DataStructureIO, List2DASCII) {
  std::stringstream test(std::ios::in | std::ios::out);
  DataStructureReader reader;
  DataStructureWriter writer;
  List2D<realT> x(2,5);
  Array1D<intT> col(2);
  col(0) = 2;
  col(1) = 3;
  x.set_ncol(col);

  x(0,0) = 1.1;
  x(0,1) = 2.2;
  x(1,0) = 3.3;
  x(1,1) = 4.4;
  x(1,2) = 5.5;

  //write out the array out
  writer.WriteList2D(test, DataStructureIO::IOMode::ASCII, x);

  //read array back in from stream
  List2D<realT> y = reader.ReadList2D<realT>
    (test, DataStructureIO::IOMode::ASCII);

  for(intT i = 0; i < 2; i++){
    for(intT j = 0; j < x.get_ncol(i); j++){
     EXPECT_DOUBLE_EQ(x(i,j), y(i,j));
    }
  }
}

