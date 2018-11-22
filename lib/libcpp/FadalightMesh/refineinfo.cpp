#include  "FadalightMesh/refineinfo.hpp"
#include  <fstream>
#include  <iostream>

using namespace FadalightMesh;

/*---------------------------------------------------------*/

RefineInfo::~RefineInfo(){}
RefineInfo::RefineInfo(){}
RefineInfo::RefineInfo(const RefineInfo& refineinfo){}
RefineInfo& RefineInfo::operator=( const RefineInfo& refineinfo)
{
  assert(0);
}

std::string RefineInfo::getClassName() const
{
  return "RefineInfo";
}

/*---------------------------------------------------------*/

alat::SparsityPattern& RefineInfo::getCoarseNodes()
{
  return _coarsenodeids;
}

alat::SparsityPattern& RefineInfo::getCoarseCells()
{
  return _coarsecellids;
}

alat::SparsityPattern& RefineInfo::getCoarseSides()
{
  return _coarsesideids;
}

alat::SparsityPattern& RefineInfo::getCoarseEdges()
{
  return _coarseedgeids;
}

const alat::SparsityPattern& RefineInfo::getCoarseNodes() const
{
  return _coarsenodeids;
}

const alat::SparsityPattern& RefineInfo::getCoarseCells() const
{
  return _coarsecellids;
}

const alat::SparsityPattern& RefineInfo::getCoarseSides() const
{
  return _coarsesideids;
}

const alat::SparsityPattern& RefineInfo::getCoarseEdges() const
{
  return _coarseedgeids;
}

int RefineInfo::getNCoarseCells(int i) const
{
  return _coarsecellids.rowsize(i);
}

int RefineInfo::getNCoarseSides(int i) const
{
  return _coarsesideids.rowsize(i);
}

int RefineInfo::getNCoarseNodes(int i) const
{
  return _coarsenodeids.rowsize(i);
}

int RefineInfo::getNCoarseEdges(int i) const
{
  return _coarseedgeids.rowsize(i);
}

int RefineInfo::getCoarseCellNumber(int i, int ii) const
{
  return _coarsecellids.get(i, ii);
}

int RefineInfo::getCoarseSideNumber(int i, int ii) const
{
  return _coarsesideids.get(i, ii);
}

int RefineInfo::getCoarseNodesNumber(int i, int ii) const
{
  return _coarsenodeids.get(i, ii);
}

int RefineInfo::getCoarseEdgesNumber(int i, int ii) const
{
  return _coarseedgeids.get(i, ii);
}

alat::armaivec& RefineInfo::getNodeIds()
{
  return _nodeids;
}

const alat::armaivec& RefineInfo::getNodeIds() const
{
  return _nodeids;
}

alat::armaivec& RefineInfo::getSideIds()
{
  return _sideids;
}

const alat::armaivec& RefineInfo::getSideIds() const
{
  return _sideids;
}

alat::armaivec& RefineInfo::getCellIds()
{
  return _cellids;
}

const alat::armaivec& RefineInfo::getCellIds() const
{
  return _cellids;
}

/*---------------------------------------------------------*/
void RefineInfo::load(std::string filename)
{
  std::ifstream file( filename.c_str() );
  assert( file.is_open() );
  _coarsenodeids.load(file);
  _coarsecellids.load(file);
  _coarsesideids.load(file);
  // _coarseedgeids.read(file);
  _nodeids.load(file);
  _sideids.load(file);
  _cellids.load(file);
  refinfoinfonode.load(file);
  refinfoinfoside.load(file);
  refinfoinfocell.load(file);
  // std::cerr << "RefineInfo::read() " << filename << " n = " << _coarsenodeids.n() << "\n";
  file.close();
}

/*---------------------------------------------------------*/
void RefineInfo::save(std::string filename, arma::file_type datatype) const
{
  // std::cerr << "RefineInfo::write() " << filename << " n = " << _coarsenodeids.n() << "\n";
  std::ofstream file( filename.c_str() );
  assert( file.is_open() );
  _coarsenodeids.save(file, datatype);
  _coarsecellids.save(file, datatype);
  _coarsesideids.save(file, datatype);
  // _coarseedgeids.write(file, datatype);
  _nodeids.save(file, datatype);
  _sideids.save(file, datatype);
  _cellids.save(file, datatype);
  refinfoinfonode.save(file, datatype);
  refinfoinfoside.save(file, datatype);
  refinfoinfocell.save(file, datatype);
  file.close();
}
