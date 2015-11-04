/*
 * DataStructureWriter.h
 *
 *  Created on: Oct 31, 2015
 *      Author: rabbit
 */

#ifndef DATASTRUCTUREWRITER_H_
#define DATASTRUCTUREWRITER_H_
#include "IO/DataStructureIO.h"
#include "DataStructures/Array1D.h"
#include "DataStructures/Array2D.h"
#include "DataStructures/List2D.h"
#include "SystemUtils/SystemModule.h"
#include <type_traits>
//****************************************************************************80
//! \class DataStructureWriter
//! \brief Class for writing the base data structures to output streams
//! \nick
//****************************************************************************80
class DataStructureWriter : public DataStructureIO{
public:

//****************************************************************************80
//! \brief WriteArray1D: Writes out the given 1d array to the provided
//!                             output stream.  Depending on the output mode,
//!                             the output written in ASCII or binary
//! \nick
//! \param[in] output_stream -output stream
//! \param[in] input_node   -output mode of stream
//! \param[in] array        -Array1D
//****************************************************************************80
  template<typename T>
  void WriteArray1D(std::ostream & output_stream,
                    IOMode output_mode,
                    Array1D<T> const & array) const
  {
    switch(output_mode){
      case IOMode::BINARY:
        WriteArray1DBinary(output_stream, array);
        break;
      case IOMode::ASCII:

        WriteArray1DASCII(output_stream, array);
        break;
      default:
        SystemModule::my_exit();
    }

  }

//****************************************************************************80
//! \brief WriteArray2D: Writes out the given 2d array to the provided
//!                             output stream.  Depending on the output mode,
//!                             the output written in ASCII or binary
//! \nick
//! \param[in] output_stream -output stream
//! \param[in] input_node   -output mode of stream
//! \param[in] array        -Array2D
//****************************************************************************80
  template<typename T>
  void WriteArray2D(std::ostream & output_stream,
                    IOMode output_mode,
                    Array2D<T> const & array) const
  {
    switch(output_mode){
      case IOMode::BINARY:
        WriteArray2DBinary(output_stream, array);
        break;
      case IOMode::ASCII:
        WriteArray2DASCII(output_stream, array);
        break;
      default:
        SystemModule::my_exit();
    }

  }

//****************************************************************************80
//! \brief WriteList2D: Writes out the given 2d list to the provided
//!                             output stream.  Depending on the output mode,
//!                             the output written in ASCII or binary
//! \nick
//! \param[in] output_stream -output stream
//! \param[in] input_node   -output mode of stream
//! \param[in] list         -List2D
//****************************************************************************80
  template<typename T>
  void WriteList2D(std::ostream & output_stream,
                   IOMode output_mode,
                   List2D<T> const & list) const
  {
    switch(output_mode){
      case IOMode::BINARY:
        WriteList2DBinary(output_stream, list);
        break;
      case IOMode::ASCII:
        WriteList2DASCII(output_stream, list);
        break;
      default:
        SystemModule::my_exit();
    }

  }

private:
  intT const print_width_ = 12; //!< Width of between entries
  intT const output_precision_ = 32; //!< Number of digits to write
//****************************************************************************80
//! \brief WriteArray1DBinary: Writes out the array to the
//!                             provided output stream in a binary format.
//! \nick
//! \param[in] output_stream -output stream
//! \param[in] array        -Array1D
//****************************************************************************80
  template<typename T>
  void WriteArray1DBinary(std::ostream & output_stream,
                          Array1D<T> const & array) const
  {
    //---> Write number of rows to make reading easier
    intT nrow = array.get_size(0);

    output_stream.write(reinterpret_cast<char*>(&nrow), sizeof(nrow));
    output_stream.

    write(reinterpret_cast<char const *>(array.begin()), sizeof(T)*nrow);

  }

//****************************************************************************80
//! \brief WriteArray2DBinary: Writes out the array to the
//!                             provided output stream in a binary format.
//! \nick
//! \param[in] output_stream -output stream
//! \param[in] array        -Array2D
//****************************************************************************80
  template<typename T>
  void WriteArray2DBinary(std::ostream & output_stream,
                          Array2D<T> const & array) const
  {
    //---> Write the number of rows and columns to make reading easier
    intT nrow = array.get_size(0);
    intT ncol = array.get_size(1);

    output_stream.write(reinterpret_cast<char*>(&nrow), sizeof(nrow));
    output_stream.write(reinterpret_cast<char*>(&ncol), sizeof(ncol));

    output_stream.
    write(reinterpret_cast<char const *>(array.begin()), sizeof(T)*nrow*ncol);
  }

//****************************************************************************80
//! \brief WriteList2DBinary: Writes out the list to the
//!                           provided output stream in a binary format.
//! \nick
//! \param[in] output_stream -output stream
//! \param[in] list         -List2D
//****************************************************************************80
  template<typename T>
  void WriteList2DBinary(std::ostream & output_stream,
                         List2D<T> const & list) const
  {
    //---> Write leading and total sizes to make allocating easy
    intT nrow = list.get_lead_size();
    intT total_size = list.get_total_size();

    output_stream.write(reinterpret_cast<char*>(&nrow), sizeof(nrow));
    output_stream.write(reinterpret_cast<char*>(&total_size),
                        sizeof(total_size));

    //---> Write the number of columns per row
    for (intT row = 0; row < nrow; row++){
      intT column_size = list.get_ncol(row);
      output_stream.write(reinterpret_cast<char*>(&column_size),
                          sizeof(column_size));
    }

    output_stream.
    write(reinterpret_cast<char const *>(list.begin()), sizeof(T)*total_size);

  }

//****************************************************************************80
//! \brief WriteArray1DASCII: Writes out the array to the
//!                             provided output stream in a ASCII format.
//! \nick
//! \param[in] output_stream -output stream
//! \param[in] array        -Array1D
//****************************************************************************80
  template<typename T>
  void WriteArray1DASCII(std::ostream & output_stream,
                         Array1D<T> const & array) const
  {
    output_stream.precision(output_precision_);
    output_stream.setf(std::ios_base::scientific, std::ios_base::floatfield);

    //---> Write number of rows to make reading easier
    intT nrow = array.get_size(0);
    output_stream << std::setw(print_width_) << nrow << "\n";

    for (intT row = 0; row < nrow; row++){
      output_stream << std::setw(print_width_) << array(row) << "\n";
    }

  }

//****************************************************************************80
//! \brief WriteArray2DASCII: Writes out the array to the
//!                             provided output stream in a ASCII format.
//! \nick
//! \param[in] output_stream -output stream
//! \param[in] array        -Array2D
//****************************************************************************80
  template<typename T>
  void WriteArray2DASCII(std::ostream & output_stream,
                         Array2D<T> const & array) const
  {

    output_stream.precision(output_precision_);
    output_stream.setf(std::ios_base::scientific, std::ios_base::floatfield);

    intT nrow = array.get_size(0);
    intT ncol = array.get_size(1);
    //---> Write the number of rows and columns to make reading easier
    output_stream << std::setw(print_width_) << nrow
        << std::setw(print_width_) << ncol
        <<  "\n";

    for (intT row = 0; row < nrow; row++){
      for (intT col = 0; col < ncol;col++){
        output_stream << std::setw(print_width_) << array(row,col) << " ";
      }
      output_stream << "\n";
    }
  }

//****************************************************************************80
//! \brief WriteList2DASCII: Writes out the list to the
//!                           provided output stream in a ASCII format.
//! \nick
//! \param[in] output_stream -output stream
//! \param[in] list         -List2D
//****************************************************************************80
  template<typename T>
  void WriteList2DASCII(std::ostream & output_stream,
                        List2D<T> const & list) const
  {
    output_stream.precision(output_precision_);
    output_stream.setf(std::ios_base::scientific, std::ios_base::floatfield);

    //---> Write leading and total sizes to make allocating easy
    intT nrow = list.get_lead_size();
    intT total_size = list.get_total_size();

    output_stream << std::setw(print_width_) << nrow
        << std::setw(print_width_) << total_size
        << "\n";

    //---> Write the number of columns per row
    for (intT row = 0; row < nrow; row++){
      output_stream << std::setw(print_width_) << list.get_ncol(row) << " ";
    }
    if(nrow >0){
      output_stream << "\n";
    }

    for (intT row = 0; row < nrow; row++){
      for (intT col = 0; col < list.get_ncol(row); col++){
        output_stream << std::setw(print_width_) << list(row,col) << " ";
      }
      output_stream << "\n";
    }

  }
};


#endif /* DATASTRUCTUREWRITER_H_ */
