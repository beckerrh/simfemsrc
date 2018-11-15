#include  "FadalightMesh/boundaryinfo.hpp"
#include  "Alat/getlinesplit.hpp"
#include  <cassert>
#include  <iostream>

using namespace FadalightMesh;
using namespace std;

/*---------------------------------------------------------*/
BoundaryInfo::~BoundaryInfo(){}
BoundaryInfo::BoundaryInfo(){}
BoundaryInfo::BoundaryInfo(const BoundaryInfo& boundaryinfo)
{
  _index_of_col = boundaryinfo._index_of_col;
  _colors = boundaryinfo._colors;
  _cells = boundaryinfo._cells;
  _sides = boundaryinfo._sides;
  _sideids_of_cells = boundaryinfo._sideids_of_cells;
}
BoundaryInfo::BoundaryInfo(const std::string& filename)
{
  read(filename);
}
BoundaryInfo& BoundaryInfo::operator=(const BoundaryInfo& boundaryinfo)
{
  assert(0);
  return *this;
}

/*---------------------------------------------------------*/

int BoundaryInfo::_IndexOfColor(int color) const
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

void BoundaryInfo::_ConstructIntOfColor()
{
  _index_of_col.clear();
  for(int ic = 0; ic < _colors.size(); ic++)
  {
    _index_of_col[_colors[ic]] = ic;
  }
}

/*---------------------------------------------------------*/

void BoundaryInfo::set_size(const std::map<int, int>& size_of_color)
{
  int ncolors = size_of_color.size();
  _colors.set_size(ncolors);
  _cells.set_size(ncolors);
  _sides.set_size(ncolors);
  _sideids_of_cells.set_size(ncolors);
  int count = 0;
  for(map<int, int>::const_iterator p = size_of_color.begin(); p != size_of_color.end(); p++)
  {
    int color = p->first;
    int size = p->second;

    _colors[count] = color;
    _cells[count].set_size(size);
    _sides[count].set_size(size);
    _sideids_of_cells[count].set_size(size);
    count++;
  }
  _ConstructIntOfColor();
}

/*---------------------------------------------------------*/

void BoundaryInfo::read(const std::string& filename)
{
  ifstream file;
  string name = filename;
//    name+=".boundary";
  file.open( name.c_str() );
  assert( file.is_open() );

  // vector<string> words = getLineSplit(file);
  // assert(words.size()==2);
  // string type(words[0]);
  // int ncolors=atoi(words[1].c_str());
  // std::cerr << "BoundaryInfo::read() " << filename << std::endl;
  // std::cerr << "BoundaryInfo::read() " << type << std::endl;
  // std::cerr << "BoundaryInfo::read() " << ncolors << std::endl;

  _colors.load(file);
  int ncolors = _colors.size();
  // std::cerr << "BoundaryInfo::read() colors: " << _colors << "\n";

  _cells.set_size(ncolors);
  _sides.set_size(ncolors);
  _sideids_of_cells.set_size(ncolors);
  for(int ic = 0; ic < ncolors; ic++)
  {
    _cells[ic].load(file);
    _sides[ic].load(file);
    _sideids_of_cells[ic].load(file);
  }
  _ConstructIntOfColor();
}

/*---------------------------------------------------------*/

void BoundaryInfo::write(const std::string& filename, arma::file_type datatype) const
{
  assert( _cells.size() == _colors.size() );
  assert( _cells.size() == _sides.size() );
  assert( _cells.size() == _sideids_of_cells.size() );
  ofstream file;
  string name = filename;
  file.open( name.c_str() );

  _colors.save(file, datatype);
  for(int ic = 0; ic < _colors.size(); ic++)
  {
    _cells[ic].save(file, datatype);
    _sides[ic].save(file, datatype);
    _sideids_of_cells[ic].save(file, datatype);
  }

  file.close();
}
//
// /*---------------------------------------------------------*/
//
// void BoundaryInfo::removeColor(int color)
// {
//   map<int, int>::iterator it = _index_of_col.find(color);
//
//   _index_of_col.erase(it);
//   int i = it->second;
//
//   _colors.erase(_colors.begin()+i);
//   _cells.erase(_cells.begin()+i);
//   _sides.erase(_sides.begin()+i);
//   _sideids_of_cells.erase(_sideids_of_cells.begin()+i);
//
//   _ConstructIntOfColor();
// }
//
// /*---------------------------------------------------------*/
//
// int BoundaryInfo::getNBoundaries()
// {
//   return _cells.size();
// }

/*---------------------------------------------------------*/

int BoundaryInfo::getNSides() const
{
  unsigned int n = 0;
  for(unsigned int i = 0; i < _sides.size(); i++)
  {
    n += _sides[i].size();
  }
  return n;
}

/*---------------------------------------------------------*/

alat::armaivec& BoundaryInfo::getColors()
{
  return _colors;
}

/*---------------------------------------------------------*/

const alat::armaivec& BoundaryInfo::getColors() const
{
  return _colors;
}

/*---------------------------------------------------------*/

alat::armaivec& BoundaryInfo::getCellsOfColor(int color)
{
  return _cells[_IndexOfColor(color)];
}

/*---------------------------------------------------------*/

const alat::armaivec& BoundaryInfo::getCellsOfColor(int color) const
{
  return _cells[_IndexOfColor(color)];
}

/*---------------------------------------------------------*/

alat::armaivec& BoundaryInfo::getSidesOfColor(int color)
{
  return _sides[_IndexOfColor(color)];
}

/*---------------------------------------------------------*/

const alat::armaivec& BoundaryInfo::getSidesOfColor(int color) const
{
  return _sides[_IndexOfColor(color)];
}

/*---------------------------------------------------------*/

alat::armaivec& BoundaryInfo::getSidesIdOfCellsOfColor(int color)
{
  return _sideids_of_cells[_IndexOfColor(color)];
}

/*---------------------------------------------------------*/

const alat::armaivec& BoundaryInfo::getSidesIdOfCellsOfColor(int color) const
{
  return _sideids_of_cells[_IndexOfColor(color)];
}
