#include "IO/UnstMeshReader.h"
//****************************************************************************80
UnstMeshReader::UnstMeshReader(const std::string& filename) : 
  filename_(filename)
{
	nnode_ = 0;
	nelement_= 0;
	nbc_face_ = 0;
	nbc_id_ = 0;
}// End UnstMeshReader::UnstMeshReader
//****************************************************************************80
UnstMeshReader::~UnstMeshReader(){}
