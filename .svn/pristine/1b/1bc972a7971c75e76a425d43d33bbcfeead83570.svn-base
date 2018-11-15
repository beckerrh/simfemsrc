#include  "FadalightMesh/curvedinteriorsideinfo.hpp"
#include  "Alat/getlinesplit.hpp"
#include  <cassert>
#include  <iostream>

using namespace FadalightMesh;
using namespace std;

/*---------------------------------------------------------*/
CurvedInteriorSideInfo::~CurvedInteriorSideInfo(){}
CurvedInteriorSideInfo::CurvedInteriorSideInfo(){}
CurvedInteriorSideInfo::CurvedInteriorSideInfo(const CurvedInteriorSideInfo& boundaryinfo)
{
  _index_of_col = boundaryinfo._index_of_col;
  _colors = boundaryinfo._colors;
  _curvedsides = boundaryinfo._curvedsides;
}

CurvedInteriorSideInfo::CurvedInteriorSideInfo(const std::string& filename)
{
  read(filename);
}

CurvedInteriorSideInfo&  CurvedInteriorSideInfo::operator=(const CurvedInteriorSideInfo& boundaryinfo)
{
  assert(0);
  return *this;
}

/*---------------------------------------------------------*/
int CurvedInteriorSideInfo::_IndexOfColor(int color) const
{
  alat::IntMap::const_iterator p = _index_of_col.find(color);
  if( p == _index_of_col.end() )
  {
    std::cerr<<"*** not found color "<<color<<std::endl;
    std::cerr<<"of "<<_index_of_col<<std::endl;
    assert(0);
  }
  return p->second;
}

/*---------------------------------------------------------*/

void CurvedInteriorSideInfo::_ConstructIntOfColor()
{
  _index_of_col.clear();
  for(int ic = 0; ic < _colors.size(); ic++)
  {
    _index_of_col[_colors[ic]] = ic;
  }
}

/*---------------------------------------------------------*/

void CurvedInteriorSideInfo::set_size(const std::map<int, int>& size_of_color)
{
  int ncolors = size_of_color.size();
  _colors.set_size(ncolors);
  // _cells.set_size(ncolors);
  _curvedsides.set_size(ncolors);
  // _sideids_of_cells.set_size(ncolors);
  int count = 0;
  for(map<int, int>::const_iterator p = size_of_color.begin(); p != size_of_color.end(); p++)
  {
    int color = p->first;
    int size = p->second;

    _colors[count] = color;
    // _cells[count].set_size(size);
    _curvedsides[count].set_size(size);
    // _sideids_of_cells[count].set_size(size);
    count++;
  }
  _ConstructIntOfColor();
}

/*---------------------------------------------------------*/

void CurvedInteriorSideInfo::read(const std::string& filename)
{
  ifstream file;
  string name = filename;
  file.open( name.c_str() );
  assert( file.is_open() );
  _colors.load(file);
  int ncolors = _colors.size();
  _curvedsides.set_size(ncolors);
  for(int ic = 0; ic < ncolors; ic++)
  {
    _curvedsides[ic].load(file);
  }
  _ConstructIntOfColor();
}

/*---------------------------------------------------------*/

void CurvedInteriorSideInfo::write(const std::string& filename, arma::file_type datatype) const
{
  ofstream file;
  string name = filename;
  file.open( name.c_str() );

  _colors.save(file, datatype);
  for(int ic = 0; ic < _colors.size(); ic++)
  {
    _curvedsides[ic].save(file, datatype);
  }

  file.close();
}

/*---------------------------------------------------------*/
int CurvedInteriorSideInfo::getNSides() const
{
  unsigned int n = 0;
  for(unsigned int i = 0; i < _curvedsides.size(); i++)
  {
    n += _curvedsides[i].size();
  }
  return n;
}

/*---------------------------------------------------------*/
alat::armaivec&  CurvedInteriorSideInfo::getColors()
{
  return _colors;
}

/*---------------------------------------------------------*/

const alat::armaivec&  CurvedInteriorSideInfo::getColors() const
{
  return _colors;
}
//
// /*---------------------------------------------------------*/
//
// alat::armaivec&  CurvedInteriorSideInfo::getCellsOfColor(int color)
// {
//   return _cells[_IndexOfColor(color)];
// }
//
// /*---------------------------------------------------------*/
//
// const alat::armaivec&  CurvedInteriorSideInfo::getCellsOfColor(int color) const
// {
//   return _cells[_IndexOfColor(color)];
// }

/*---------------------------------------------------------*/

alat::armaivec&  CurvedInteriorSideInfo::getSidesOfColor(int color)
{
  return _curvedsides[_IndexOfColor(color)];
}
const alat::armaivec&  CurvedInteriorSideInfo::getSidesOfColor(int color) const
{
  return _curvedsides[_IndexOfColor(color)];
}
//
// /*---------------------------------------------------------*/
//
// alat::armaivec&  CurvedInteriorSideInfo::getSidesIdOfCellsOfColor(int color)
// {
//   return _sideids_of_cells[_IndexOfColor(color)];
// }
//
// /*---------------------------------------------------------*/
//
// const alat::armaivec&  CurvedInteriorSideInfo::getSidesIdOfCellsOfColor(int color) const
// {
//   return _sideids_of_cells[_IndexOfColor(color)];
// }
