/*
 * DataStructureReader.test.cpp
 *
 *  Created on: Nov 3, 2015
 *      Author: rabbit
 */
#include "gtest/gtest.h"
#include "IO/DataStructureReader.h"
//--------------------------- BINARY READING TESTS -----------------------------

//----------------> NON-PRECONSTRUCTED DATA STRUCTURES
TEST(DataStructureReader_Binary, ReadArray1D) {
  std::stringstream test(std::ios::in | std::ios::out | std::ios::binary);
  Array1D<double> x(3);
  x(0) = 1.1;
  x(1) = 2.2;
  x(2) = 3.3;
  int size = x.get_size(0);

  test.write(reinterpret_cast<char*>(&size), sizeof(size));
  test.write(reinterpret_cast<char*>(&x(0)), sizeof(x(0))*size);

  Array1D<double> y;
  DataStructureReader reader;
  y = reader.ReadArray1D<double>(test, DataStructureIO::IOMode::BINARY);

  EXPECT_DOUBLE_EQ(x(0), y(0));
  EXPECT_DOUBLE_EQ(x(1), y(1));
  EXPECT_DOUBLE_EQ(x(2), y(2));
}

TEST(DataStructureReader_Binary, ReadArray2D) {
  std::stringstream test(std::ios::in | std::ios::out | std::ios::binary);
  Array2D<double> x(2,3);
  x(0,0) = 1.1;
  x(0,1) = 2.1;
  x(0,2) = 3.1;

  x(1,0) = 1.2;
  x(1,1) = 2.2;
  x(1,2) = 3.2;
  int nrow = x.get_size(0);
  int ncol = x.get_size(1);

  test.write(reinterpret_cast<char*>(&nrow), sizeof(nrow));
  test.write(reinterpret_cast<char*>(&ncol), sizeof(ncol));
  test.write(reinterpret_cast<char*>(&x(0,0)), sizeof(x(0,0))*nrow*ncol);

  Array2D<double> y;
  DataStructureReader reader;
  y = reader.ReadArray2D<double>
    (test ,DataStructureIO::IOMode::BINARY);

  EXPECT_DOUBLE_EQ(x(0,0), y(0,0));
  EXPECT_DOUBLE_EQ(x(0,1), y(0,1));
  EXPECT_DOUBLE_EQ(x(0,2), y(0,2));

  EXPECT_DOUBLE_EQ(x(1,0), y(1,0));
  EXPECT_DOUBLE_EQ(x(1,1), y(1,1));
  EXPECT_DOUBLE_EQ(x(1,2), y(1,2));
}

TEST(DataStructureReader_Binary, ReadList2D) {
  std::stringstream test(std::ios::in | std::ios::out | std::ios::binary);
  List2D<double> x(2,5);
  Array1D<int> col(2);
  col(0) = 2;
  col(1) = 3;
  x.set_ncol(col);

  x(0,0) = 1.1;
  x(0,1) = 2.1;

  x(1,0) = 1.2;
  x(1,1) = 2.2;
  x(1,2) = 3.2;

  int nrow = x.get_lead_size();
  int total_size = x.get_total_size();

  test.write(reinterpret_cast<char*>(&nrow), sizeof(nrow));
  test.write(reinterpret_cast<char*>(&total_size), sizeof(total_size));

  for(int i = 0; i < nrow; i++){
    test.write(reinterpret_cast<char*>(&col(i)), sizeof(col(0)));
  }

  test.write(reinterpret_cast<char*>(&x(0,0)), sizeof(x(0,0))*total_size);

  List2D<double> y;
  DataStructureReader reader;
  y = reader.ReadList2D<double>
    (test, DataStructureIO::IOMode::BINARY);

  EXPECT_DOUBLE_EQ(x(0,0), y(0,0));
  EXPECT_DOUBLE_EQ(x(0,1), y(0,1));

  EXPECT_DOUBLE_EQ(x(1,0), y(1,0));
  EXPECT_DOUBLE_EQ(x(1,1), y(1,1));
  EXPECT_DOUBLE_EQ(x(1,2), y(1,2));
}

//--------------------------- ASCII READING TESTS ------------------------------

//----------------> NON-PRECONSTRUCTED DATA STRUCTURES
TEST(DataStructureReader_ASCII, ReadArray1D) {
  std::string txt("#garbage \n"
             "3 \n"
             "1.1 \n"
             "#garbage \n"
             "2.2 \n"
             "3.3 \n");
  std::istringstream test(txt);

  Array1D<double> x;
  DataStructureReader reader;
  x = reader.ReadArray1D<double>
    (test, DataStructureIO::IOMode::ASCII);

  EXPECT_EQ(1.1, x(0));
  EXPECT_EQ(2.2, x(1));
  EXPECT_EQ(3.3, x(2));
}

TEST(DataStructureReader_ASCII, ReadArray2D) {
  std::string txt("2 3 \n 1.2 2.3 3.4 \n 4.5 5.6 6.7");
  std::istringstream test(txt);

  Array2D<double> x;
  DataStructureReader reader;
  x = reader.ReadArray2D<double>
    (test, DataStructureIO::IOMode::ASCII);

  EXPECT_DOUBLE_EQ(1.2, x(0,0));
  EXPECT_DOUBLE_EQ(2.3, x(0,1));
  EXPECT_DOUBLE_EQ(3.4, x(0,2));
  EXPECT_DOUBLE_EQ(4.5, x(1,0));
  EXPECT_DOUBLE_EQ(5.6, x(1,1));
  EXPECT_DOUBLE_EQ(6.7, x(1,2));
}


TEST(DataStructureReader_ASCII, ReadList2D) {
  std::string txt("2 7 \n 3 4 \n 1.2 2.3 3.4 \n 4.5 5.6 6.7 7.8");
  std::istringstream test(txt);

  List2D<double> x;
  DataStructureReader reader;
  x = reader.ReadList2D<double>
    (test, DataStructureIO::IOMode::ASCII);
  EXPECT_DOUBLE_EQ(1.2, x(0,0));
  EXPECT_DOUBLE_EQ(2.3, x(0,1));
  EXPECT_DOUBLE_EQ(3.4, x(0,2));
  EXPECT_DOUBLE_EQ(4.5, x(1,0));
  EXPECT_DOUBLE_EQ(5.6, x(1,1));
  EXPECT_DOUBLE_EQ(6.7, x(1,2));
  EXPECT_DOUBLE_EQ(7.8, x(1,3));
}



//---------------------- BINARY READING DEATH TESTS ----------------------------

//----------------> NON-PRECONSTRUCTED DATA STRUCTURES
TEST(DataStructureReader_Binary_DeathTest, ReadArray1DBinary_NoHeader){
  //death test to make sure program fails when stream has too few rows
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  typedef Array1D<int> Array1D;
  DataStructureReader reader;

  Array1D(DataStructureReader::*func)
    (std::istream &, DataStructureIO::IOMode) =
    &DataStructureReader::ReadArray1D<int>;
  std::string txt("");
  std::istringstream test(txt);

  EXPECT_EXIT(Array1D x = (reader.*func)
              (test, DataStructureIO::IOMode::BINARY);,
               ::testing::ExitedWithCode(EXIT_FAILURE),
              "ERROR:  In DataStructureReader.h - Error while "
              "reading the header of a Array1D from binary stream");
}
TEST(DataStructureReader_Binary_DeathTest, ReadArray2DBinary_NoHeader){
  //death test to make sure program fails when stream has too few rows
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  typedef Array2D<int> Array2D;
  DataStructureReader reader;

  Array2D(DataStructureReader::*func)
    (std::istream &, DataStructureIO::IOMode) =
    &DataStructureReader::ReadArray2D<int>;
  std::string txt("");
  std::istringstream test(txt);

  EXPECT_EXIT(Array2D x = (reader.*func)
              (test, DataStructureIO::IOMode::BINARY);,
               ::testing::ExitedWithCode(EXIT_FAILURE),
              "ERROR:  In DataStructureReader.h - Error while "
              "reading the header of a Array2D from binary stream");
}

//---------------------- ASCII READING DEATH TESTS -----------------------------

//----------------> NON-PRECONSTRUCTED DATA STRUCTURES
TEST(DataStructureReaderDeathTest, ReadArray1D_NoHeader){
  //death test to make sure program fails when stream has too few rows
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  typedef Array1D<int> Array1D;
  DataStructureReader reader;

  Array1D(DataStructureReader::*func)
    (std::istream &, DataStructureIO::IOMode) =
    &DataStructureReader::ReadArray1D<int>;
  std::string txt("");
  std::istringstream test(txt);

  EXPECT_EXIT(Array1D x = (reader.*func)
              (test, DataStructureIO::IOMode::ASCII);,
               ::testing::ExitedWithCode(EXIT_FAILURE),
              "ERROR:  In DataStructureReader.h - Error while "
              "reading the header of a Array1D from stream");
}

TEST(DataStructureReaderDeathTest, ReadArray2D_NoHeader){
  //death test to make sure program fails when stream has too few rows
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  typedef Array2D<int> Array2D;
  DataStructureReader reader;
  Array2D(DataStructureReader::*func)
    (std::istream &, DataStructureIO::IOMode) =
    &DataStructureReader::ReadArray2D<int>;
  std::string txt("");
  std::istringstream test(txt);

  EXPECT_EXIT(Array2D x = (reader.*func)
              (test, DataStructureIO::IOMode::ASCII);,
               ::testing::ExitedWithCode(EXIT_FAILURE),
              "ERROR:  In DataStructureReader.h - Error while "
              "reading the header of a Array2D from stream");
}

TEST(DataStructureReaderDeathTest, ReadList2D_NoHeader){
  //death test to make sure program fails when stream has too few rows
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  typedef List2D<int> List2D;
  DataStructureReader reader;
  List2D(DataStructureReader::*func)
    (std::istream &, DataStructureIO::IOMode) =
    &DataStructureReader::ReadList2D<int>;
  std::string txt("");
  std::istringstream test(txt);

  EXPECT_EXIT(List2D x = (reader.*func)
              (test, DataStructureIO::IOMode::ASCII);,
               ::testing::ExitedWithCode(EXIT_FAILURE),
              "ERROR:  In DataStructureReader.h - Error while "
              "reading the header of a List2D from stream");
}



