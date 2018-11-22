#include  "FadalightMesh/hangingnodeinfo.hpp"
#include  "FadalightMesh/meshinterface.hpp"
#include  "Alat/sparsitypattern.hpp"
#include  <fstream>
#include  <list>
#include  <map>
#include  <set>

using namespace FadalightMesh;

/*---------------------------------------------------------*/
HangingNodeInfo::~HangingNodeInfo() {}
HangingNodeInfo::HangingNodeInfo() : FadalightMesh::GeometryObject(),  alat::SparsityPattern() {}
HangingNodeInfo::HangingNodeInfo(const HangingNodeInfo& hangingnodeinfo) : FadalightMesh::GeometryObject(hangingnodeinfo), alat::SparsityPattern(hangingnodeinfo) {}
HangingNodeInfo& HangingNodeInfo::operator=(const HangingNodeInfo& hangingnodeinfo)
{
  assert(0);
}

std::string HangingNodeInfo::getClassName() const
{
  return "HangingNodeInfo";
}

int HangingNodeInfo::getCellNumber(int i) const
{
  return col( rowstart(i) );
}

int HangingNodeInfo::getLocalSide(int i) const
{
  return col(rowstart(i)+1);
}

int HangingNodeInfo::getNumberOfHangingNodes(int i) const
{
  return rowstart(i+1)-rowstart(i)-2;
}

int HangingNodeInfo::getHangingNodes(int i, int in) const
{
  return col(rowstart(i)+in+2);
}

int& HangingNodeInfo::getCellNumber(int i)
{
  return alat::SparsityPattern::_col[rowstart(i)];
  // return col( rowstart(i) );
}

int& HangingNodeInfo::getLocalSide(int i)
{
  return alat::SparsityPattern::_col[rowstart(i)+1];
  // return col(rowstart(i)+1);
}

int& HangingNodeInfo::getHangingNodes(int i, int in)
{
  return alat::SparsityPattern::_col[rowstart(i)+in+2];
  // return col(rowstart(i)+in+2);
}

void HangingNodeInfo::set_size(int n_hanging, int n_local_data)
{
  int n_data = n_local_data*n_hanging;
  col().set_size(n_data);
  rowstart().set_size(n_hanging+1);
  for(int i = 0; i <= n_hanging; i++)
  {
    // rowstart(i) = i*n_local_data;
    _rowstart[i] = i*n_local_data;
  }
}

void HangingNodeInfo::load(std::string filename )
{
  std::ifstream file( filename.c_str() );
  assert( file.is_open() );
  col().load(file);
  rowstart().load(file);
  file.close();
}

void HangingNodeInfo::save(std::string filename, arma::file_type datatype) const
{
  std::ofstream file( filename.c_str() );
  assert( file.is_open() );
  col().save(file, datatype);
  rowstart().save(file, datatype);
  file.close();
}
