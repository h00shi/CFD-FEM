/*
 * DataStructureIO.h
 *
 *  Created on: Oct 31, 2015
 *      Author: rabbit
 */

#ifndef DATASTRUCTUREIO_H_
#define DATASTRUCTUREIO_H_

//****************************************************************************80
//! \class DataStructureIO
//! \brief Base class for DataStructure Input and Output
//! \nick
//****************************************************************************80
class DataStructureIO{
public:
  //enumeration for the stream mode of the input/output - binary or ASCII
  enum class IOMode{BINARY,ASCII};
  virtual ~DataStructureIO(){}
};

#endif /* DATASTRUCTUREIO_H_ */
