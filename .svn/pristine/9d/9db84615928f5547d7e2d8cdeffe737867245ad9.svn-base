#include  "FadalightMesh/hangingsideinfo.hpp"
#include  "FadalightMesh/meshinterface.hpp"
#include  <fstream>
#include  <list>
#include  <map>
#include  <set>

using namespace FadalightMesh;

/*---------------------------------------------------------*/
HangingSideInfo::~HangingSideInfo() {}
HangingSideInfo::HangingSideInfo() : FadalightMesh::GeometryObject(), alat::SparsityPattern() {}
HangingSideInfo::HangingSideInfo(const HangingSideInfo& hangingsideinfo) : FadalightMesh::GeometryObject(hangingsideinfo), alat::SparsityPattern(hangingsideinfo) {}
HangingSideInfo& HangingSideInfo::operator=(const HangingSideInfo& hangingsideinfo)
{
  assert(0);
}
std::string HangingSideInfo::getClassName() const
{
  return "HangingSideInfo";
}

/*---------------------------------------------------------*/
int HangingSideInfo::getCellNumber(int i) const
{
  return col( rowstart(i) );
}

int HangingSideInfo::getLocalSide(int i) const
{
  return col(rowstart(i)+1);
}

int HangingSideInfo::getNumberOfHangingSides(int i) const
{
  return rowstart(i+1)-rowstart(i)-2;
}

int HangingSideInfo::getHangingSides(int i, int in) const
{
  return col(rowstart(i)+in+2);
}

int& HangingSideInfo::getCellNumber(int i)
{
  return alat::SparsityPattern::_col[rowstart(i)];
  // return col( rowstart(i) );
}

int& HangingSideInfo::getLocalSide(int i)
{
  return alat::SparsityPattern::_col[rowstart(i)+1];
  // return col(rowstart(i)+1);
}

int& HangingSideInfo::getHangingSides(int i, int in)
{
  return alat::SparsityPattern::_col[rowstart(i)+in+2];
  // return col(rowstart(i)+in+2);
}

void HangingSideInfo::set_size(int n_hanging, int n_local_data)
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

std::vector<int> HangingSideInfo::getNeighbourCellsOfSide(int iK, int iis, const FadalightMesh::MeshInterface* mesh)
{
  std::vector<int> _neighbourcells(0);
  int iS = mesh->getSideIdOfCell(iK, iis);
  for(int ih = 0; ih < n(); ih++)
  {
    int iKh = getCellNumber(ih);
    int ils = getLocalSide(ih);
    int nhs = getNumberOfHangingSides(ih);
    if( ( iKh == iK ) &&( ils == iis ) ) // iK est la maille grossière
    {
      for(int ihs = 0; ihs < nhs; ihs++)
      {
        int isfine = getHangingSides(ih, ihs);
        int iKneighbour = mesh->getCellIdOfSide(isfine, 0);
        assert(iKneighbour >= 0);
        _neighbourcells.push_back(iKneighbour);
      }
      return _neighbourcells;
    }
    else
    {
      for(int ihs = 0; ihs < nhs; ihs++)
      {
        int isfine = getHangingSides(ih, ihs);
        int iKneighbour = mesh->getCellIdOfSide(isfine, 0);
        assert(iKneighbour >= 0);
        if( ( iKneighbour == iK ) &&( isfine == iS ) ) //iK est une maille raffinée
        {
          _neighbourcells.push_back(iKh);
          return _neighbourcells;
        }
      }
    }
  }
  //assert(0);
  return _neighbourcells;
}

void HangingSideInfo::load(std::string filename )
{
  std::ifstream file( filename.c_str() );
  assert( file.is_open() );
  col().load(file);
  rowstart().load(file);
  file.close();
}

void HangingSideInfo::save(std::string filename, arma::file_type datatype) const
{
  std::ofstream file( filename.c_str() );
  assert( file.is_open() );
  col().save(file, datatype);
  rowstart().save(file, datatype);
  file.close();
}
