#include  "FadalightMesh/coarseninfo.hpp"
#include  <fstream>
#include  <iostream>

using namespace FadalightMesh;

/*---------------------------------------------------------*/
CoarsenInfo::~CoarsenInfo(){}
CoarsenInfo::CoarsenInfo(){}
CoarsenInfo::CoarsenInfo(const CoarsenInfo& coarseninfo){}
CoarsenInfo& CoarsenInfo::operator=( const CoarsenInfo& coarseninfo)
{
  assert(0);
}

std::string CoarsenInfo::getClassName() const
{
  return "CoarsenInfo";
}

/*---------------------------------------------------------*/

alat::SparsityPattern& CoarsenInfo::getOldNodes()
{
  return _oldnodeids;
}

alat::SparsityPattern& CoarsenInfo::getOldCells()
{
  return _oldcellids;
}

alat::SparsityPattern& CoarsenInfo::getOldSides()
{
  return _oldsideids;
}

alat::SparsityPattern& CoarsenInfo::getOldEdges()
{
  return _oldedgeids;
}

int CoarsenInfo::getNNodes() const
{
  return _oldnodeids.n();
}

int CoarsenInfo::getNOldNodes(int i) const
{
  return _oldnodeids.rowsize(i);
}

int CoarsenInfo::getOldNode(int i, int ii) const
{
  return _oldnodeids.col(_oldnodeids.rowstart(i)+ii);
}

int CoarsenInfo::getNCells() const
{
  return _oldcellids.n();
}

int CoarsenInfo::getNOldCells(int i) const
{
  return _oldcellids.rowsize(i);
}

int CoarsenInfo::getOldCell(int i, int ii) const
{
  return _oldcellids.col(_oldcellids.rowstart(i)+ii);
}

int CoarsenInfo::getNSides() const
{
  return _oldsideids.n();
}

int CoarsenInfo::getNOldSides(int i) const
{
  return _oldsideids.rowsize(i);
}

int CoarsenInfo::getOldSide(int i, int ii) const
{
  return _oldsideids.col(_oldsideids.rowstart(i)+ii);
}

int CoarsenInfo::getNEdges() const
{
  return _oldedgeids.n();
}

int CoarsenInfo::getNOldEdges(int i) const
{
  return _oldedgeids.rowsize(i);
}

int CoarsenInfo::getOldEdge(int i, int ii) const
{
  return _oldedgeids.col(_oldedgeids.rowstart(i)+ii);
}

alat::armaivec& CoarsenInfo::getNodeNewToOld()
{
  return _nodenewtoold;
}

const alat::armaivec& CoarsenInfo::getNodeNewToOld() const
{
  return _nodenewtoold;
}

/*---------------------------------------------------------*/
void CoarsenInfo::load(std::string filename)
{
  std::ifstream file( filename.c_str() );
  assert( file.is_open() );
  _oldnodeids.load(file);
  _oldcellids.load(file);
  _oldsideids.load(file);
  // _oldedgeids.read(file);
  _nodenewtoold.load(file);
  file.close();
}

/*---------------------------------------------------------*/
void CoarsenInfo::save(std::string filename, arma::file_type datatype) const
{
  std::ofstream file( filename.c_str() );
  assert( file.is_open() );
  _oldnodeids.save(file, datatype);
  _oldcellids.save(file, datatype);
  _oldsideids.save(file, datatype);
  // _oldedgeids.write(file, datatype);
  _nodenewtoold.save(file, datatype);
  file.close();
}
