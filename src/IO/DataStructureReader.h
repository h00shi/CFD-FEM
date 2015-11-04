/*
 * DataStructureReader.h
 *
 *  Created on: Nov 1, 2015
 *      Author: rabbit
 */

#ifndef DATASTRUCTUREREADER_H_
#define DATASTRUCTUREREADER_H_

#include "IO/DataStructureIO.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/List2D.h"
#include "SystemUtils/SystemModule.h"
#include <type_traits>
#include <algorithm>
//****************************************************************************80
//! \class DataStructureReader
//! \brief Class for reading the base data structures from input streams
//! \nick
//****************************************************************************80
class DataStructureReader: public DataStructureIO{
public:

//****************************************************************************80
//! \brief ReadArray1D: Reads in the given 1D array from the provided
//!                             input stream.  Depending on the IO mode,
//!                             the output read in ascii or binary
//! \nick
//! \param[in] input_stream -output stream
//! \param[in] input_node   -output mode of stream
//! \return array        -Array1D
//****************************************************************************80
  template<typename T>
  Array1D<T> ReadArray1D(std::istream & input_stream,
                    IOMode input_mode)
  {
    Array1D<T> array;
    switch(input_mode){
      case IOMode::BINARY:
        array = ReadArray1DBinary<T>(input_stream);
        break;
      case IOMode::ASCII:
        array = ReadArray1DASCII<T>(input_stream);
        break;
      default:
        SystemModule::my_exit();
    }
    return array;
  }

//****************************************************************************80
//! \brief ReadArray2D: Reads in the given 2D array from the provided
//!                             input stream.  Depending on the IO mode,
//!                             the output read in ascii or binary
//! \nick
//! \param[in] input_stream -output stream
//! \param[in] input_node   -output mode of stream
//! \return array        -Array2D
//****************************************************************************80
  template<typename T>
  Array2D<T> ReadArray2D(std::istream & input_stream,
                    IOMode input_mode)
  {
    Array2D<T> array;
    switch(input_mode){
      case IOMode::BINARY:
        array = ReadArray2DBinary<T>(input_stream);
        break;
      case IOMode::ASCII:
        array = ReadArray2DASCII<T>(input_stream);
        break;
      default:
        SystemModule::my_exit();
    }
    return array;
  }

//****************************************************************************80
//! \brief ReadList2D: Reads in the given 2D List from the provided
  //!                             input stream.  Depending on the IO mode,
  //!                             the output read in ascii or binary
//! \nick
//! \param[in] input_stream -output stream
//! \param[in] input_node   -output mode of stream
//! \return list         -List2D
//****************************************************************************80
  template<typename T>
  List2D<T> ReadList2D(std::istream & input_stream,
                   IOMode input_mode)
  {
    List2D<T> list;
    switch(input_mode){
      case IOMode::BINARY:
        list = ReadList2DBinary<T>(input_stream);
        break;
      case IOMode::ASCII:
        list = ReadList2DASCII<T>(input_stream);
        break;
      default:
        SystemModule::my_exit();
    }
    return list;
  }

private:
  //intT const print_width_ = 12; //!< Width of between entries
 // intT const input_precision_ = 32; //!< Number of digits to write
//****************************************************************************80
//! \brief ReadArray1DBinary: Read in an array from the
//!                             provided input stream in a binary format.
//! \nick
//! \param[in] input_stream -input stream
//****************************************************************************80
  template<typename T>
  Array1D<T> ReadArray1DBinary(std::istream & input_stream)
  {
    intT nrow;
    //---> Read in header
    if(!input_stream.read(reinterpret_cast<char*>(&nrow), sizeof(intT))){
      std::cerr << "ERROR:  In DataStructureReader.h - "
          << "Error while reading the header of a Array1D "
          << "from binary stream" << std::endl;
      SystemModule::my_exit();
    }
    Array1D<T> array(nrow);
    if(!input_stream.read(reinterpret_cast<char*>(array.begin()),
                          nrow*sizeof(T))){
      std::cerr << "ERROR:  In DataStructureReader.h - "
          << "Error while reading Array1D data from binary stream."
          << std::endl;
      SystemModule::my_exit();
    }
    return array;
  }

//****************************************************************************80
//! \brief ReadArray2DBinary: Read in an array from the
//!                             provided output stream in a binary format.
//! \nick
//! \param[in] input_stream -input stream
//****************************************************************************80
  template<typename T>
  Array2D<T> ReadArray2DBinary(std::istream & input_stream)
  {
    intT nrow, ncol;

    //---> Read in header
    bool nrow_bool =
        input_stream.read(reinterpret_cast<char*>(&nrow), sizeof(intT));
    bool ncol_bool =
        input_stream.read(reinterpret_cast<char*>(&ncol), sizeof(intT));

    if(!(nrow_bool || ncol_bool)){
      std::cerr << "ERROR:  In DataStructureReader.h - "
          << "Error while reading the header of a Array2D "
          << "from binary stream" << std::endl;
      SystemModule::my_exit();
    }
    Array2D<T> array(nrow, ncol);
    if(!input_stream.read(reinterpret_cast<char*>(array.begin()),
                          nrow*ncol*sizeof(T))){
      std::cerr << "ERROR:  In Array2D.h - "
          << "Error while reading Array2D data from binary stream."
          << std::endl;
      SystemModule::my_exit();
    }

    return array;
  }

//****************************************************************************80
//! \brief ReadList2DBinary: Reads in the list from the
//!                           provided input stream in a binary format.
//! \nick
//! \param[in] input_stream -output stream
//****************************************************************************80
  template<typename T>
  List2D<T> ReadList2DBinary(std::istream & input_stream)
  {
    intT nrow, total_size;

    //---> Read in header
    bool nrow_bool =
        input_stream.read(reinterpret_cast<char*>(&nrow), sizeof(intT));
    bool total_size_bool =
        input_stream.read(reinterpret_cast<char*>(&total_size), sizeof(intT));

    if(!(nrow_bool || total_size_bool)){
      std::cerr << "ERROR:  In DataStructureReader.h - "
          << "Error while reading the header of a List2D "
          << "from binary stream" << std::endl;
      SystemModule::my_exit();
    }
    Array1D<intT> col_spec(nrow); //column spec for list
    if(!(input_stream.read(reinterpret_cast<char*>(col_spec.begin()),
                           nrow*sizeof(intT)))){
      std::cerr << "ERROR:  In DataStructureReader.h - "
          << "Error while reading List2D column specification "
          << "from binary stream."
          << std::endl;
      SystemModule::my_exit();
    }

    List2D<T> list(col_spec);

    if(!input_stream.read(reinterpret_cast<char*>(list.begin()),
                          total_size*sizeof(T))){
      std::cerr << "ERROR:  In DataStructureReader.h - "
          << "Error while reading List2D data from binary stream."
          << std::endl;
      SystemModule::my_exit();
    }
    return list;
  }
//****************************************************************************80
//! \brief GetLineWithHashComment:
//!        Reads lines using the specified delimiter until a line is read
//!        that does not start with "#".
//! \details By default, lines are read delimited by newline characters
//! \jun
//! \version $Rev$
//! \param[in] input_stream input stream
//! \param[in] line         first string found that does not start with "#"
//****************************************************************************80
  std::istream & GetLineIgnoreHash(std::istream & input_stream,
                                   std::string & line, char delim = '\n'){
    //---> Grab Line continuously until a # character is found
    while(true){
      if(!std::getline(input_stream, line, delim)){
        //---> prematurely end reading lines if we encounter a problem
        //     reading a line
        return input_stream;
      }
      auto first_nonspace = std::find_if(line.begin(), line.end(),
                             std::not1(std::ptr_fun<int, int>(isspace)));
      int count = std::distance(line.begin(), first_nonspace);

      //---> If we find a line that does not start with
      //     # (ignoring white space)...
      if(line[count] != '#'){
        break;
      }
    }
    return input_stream;
  }

//****************************************************************************80
//! \brief ReadArray1DASCII: Reads in the array from the
//!                             provided input stream in a ASCII format.
//! \nick
//! \param[in] input_stream -output stream
//! \param[in] array        -Array1D
//****************************************************************************80
  template<typename T>
  Array1D<T> ReadArray1DASCII(std::istream & input_stream)
  {
    intT nrow;
    std::istringstream line_stream; //string stream for a given line
    std::string line;

    //--->Read the header of the stream and get number of entries
    if(!GetLineIgnoreHash(input_stream, line)){
      std::cerr << "ERROR:  In DataStructureReader.h - "
          << "Error while reading the header of a Array1D "
          << "from stream" << std::endl;
      SystemModule::my_exit();
    }
    line_stream.str(line); //associate stream with line
    line_stream >> nrow;

    //--->Declare Array1D
    Array1D<T> array(nrow);
    //---> Read in rows
    for (intT row = 0; row < nrow; row++){ //read ith entry
      //get row
      if(!GetLineIgnoreHash(input_stream, line)){
        std::cerr << "ERROR:  In DataStructureReader.h - "
            << "Error while reading Array1D stream. "
            << "Can't read line "<<row<<" from stream."
            << std::endl;
        SystemModule::my_exit();
      }

      line_stream.clear(); //clear line stream
      line_stream.str(line); //associate stream with new line
      line_stream >> array(row);
    } //end read loop
    return array;

  }

//****************************************************************************80
//! \brief ReadArray2DASCII: Reads in the array from the
//!                             provided input stream in a ASCII format.
//! \nick
//! \param[in] input_stream -output stream
//****************************************************************************80
  template<typename T>
  Array2D<T> ReadArray2DASCII(std::istream & input_stream)
  {
    intT nrow, ncol;
    std::istringstream line_stream; //string stream for a given line
    std::string line;

    //--->Read the header of the stream and get number of entries
    if(!GetLineIgnoreHash(input_stream, line)){

      std::cerr << "ERROR:  In DataStructureReader.h - "
          << "Error while reading the header of a Array2D "
          << "from stream" << std::endl;
      std::cout << " input_stream: " << input_stream << std::endl;
      SystemModule::my_exit();
    }

    line_stream.str(line); //associate stream with line
    line_stream >> nrow;
    line_stream >> ncol;

    //--->Declare Array2D
    Array2D<T> array(nrow, ncol);
    //---> Read in rows
    for (intT row = 0; row < nrow; row++){ //read rows
      //get row
      if(!GetLineIgnoreHash(input_stream, line)){
        std::cerr << "ERROR:  In DataStructureReader.h - "
            << "Error while reading Array2D stream. "
            << "Can't read line "<<row<<" from stream."
            << std::endl;
        SystemModule::my_exit();
      }

      line_stream.clear(); //clear line stream
      line_stream.str(line); //associate stream with new line

      for (intT col = 0; col < ncol; col++) { //read columns
        bool readSuccess = (line_stream >> array(row,col));

        if(!readSuccess){ //check if read was successful
          std::cerr << "ERROR:  In DataStructureReader.h - "
              << "Error while reading Array2D stream. "
              << "Can't read column "<< col
              << " of row " << row << " from stream"
              << std::endl;
          SystemModule::my_exit();
        } //end read check
      } //end column loop
    } //end row loop
    return array;
  }

//****************************************************************************80
//! \brief ReadList2DASCII: Read in the list from the
//!                           provided in stream in a ASCII format.
//! \nick
//! \param[in] input_stream -output stream
//****************************************************************************80
  template<typename T>
  List2D<T> ReadList2DASCII(std::istream & input_stream)
  {
    intT nrow, size;
    std::istringstream line_stream; //string stream for a given line
    std::string line;

    //--->Read the header of the stream and get number of entries
    if(!GetLineIgnoreHash(input_stream, line)){
      std::cerr << "ERROR:  In DataStructureReader.h - "
          << "Error while reading the header of a List2D "
          << "from stream" << std::endl;
      SystemModule::my_exit();
    }

    line_stream.str(line); //associate stream with line
    line_stream >> nrow;
    line_stream >> size;

    //--->Declare List2D
    List2D<T> list(nrow, size);
    Array1D<intT> col_spec(nrow); //column spec for list
    col_spec.set_value(0);
    if(nrow > 0){
      //---> Read in column specification of list 2d
      if(!GetLineIgnoreHash(input_stream, line)){
        std::cerr << "ERROR:  In DataStructureReader.h - "
            << "Error while reading List2D column specification "
            << "header from stream."
            << std::endl;
        SystemModule::my_exit();
      }
      line_stream.clear(); //clear line stream
      line_stream.str(line); //associate stream with this line

      for(intT i = 0; i < nrow; i++){
        if(!((line_stream >> col_spec(i)))){
          std::cerr << "ERROR:  In DataStructureReader.h - "
              << "Error while reading List2D column specification "
              << "from stream."
              << std::endl;
          SystemModule::my_exit();
        }
      }
    }
    list.set_ncol(col_spec);

    //---> Read in actual data
    for(intT i = 0; i < nrow; i++){
      if(!GetLineIgnoreHash(input_stream, line)){
        std::cerr << "ERROR:  In DataStructureReader.h - "
            << "Error while reading List2D. "
            << "Cannot read line " << i << " "
            << "from stream."
            << std::endl;
        SystemModule::my_exit();
      }
      line_stream.clear(); //clear line stream
      line_stream.str(line); //associate stream with this line
      intT column_size = col_spec(i);
      for(int j = 0; j < column_size; j++){
        T data; //dat

        if(!(line_stream >> data)){
          std::cerr << "ERROR:  In DataStructureReader.h - "
              << "Error while reading List2D. "
              << "Cannot read column " << j << " of "
              << "line " << i << " "
              << "from stream."
              << std::endl;
          SystemModule::my_exit();
        }
        list(i,j) = data;
      }
    }
    return list;
  }
};


#endif /* DATASTRUCTUREREADER_H_ */
