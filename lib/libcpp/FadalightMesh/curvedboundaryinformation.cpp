#include  "FadalightMesh/boundaryinfo.hpp"
#include  "FadalightMesh/curvedinteriorsideinfo.hpp"
#include  "FadalightMesh/meshinterface.hpp"
#include  "FadalightMesh/curvedboundarydescriptionconstructor.hpp"
#include  "FadalightMesh/curvedboundaryinformation.hpp"
#include  "Alat/stringvector.hpp"
#include  "Alat/tokenize.hpp"
#include  <stdlib.h>

using namespace FadalightMesh;
using namespace std;

/*--------------------------------------------------------------------------*/
CurvedBoundaryInformation::~CurvedBoundaryInformation() {}
CurvedBoundaryInformation::CurvedBoundaryInformation() {}
CurvedBoundaryInformation::CurvedBoundaryInformation( const CurvedBoundaryInformation& curvedboundaryinformation)
{
  assert(0);
}

CurvedBoundaryInformation& CurvedBoundaryInformation::operator=( const CurvedBoundaryInformation& curvedboundaryinformation)
{
  const alat::Map<int, FadalightMesh::CurvedBoundaryDescriptionInterface*>& curvedboundaries = curvedboundaryinformation._curvedboundaries;
  for(alat::Map<int, FadalightMesh::CurvedBoundaryDescriptionInterface*>::const_iterator p = curvedboundaries.begin(); p != curvedboundaries.end(); p++)
  {
    _curvedboundaries[p->first] = p->second->clone();
  }
  const alat::armaivec& celliscurved = curvedboundaryinformation._celliscurved;
  _celliscurved.set_size( celliscurved.size() );
  _celliscurved = celliscurved;
  const alat::SparsityPatternFixArray<2>& spfa = curvedboundaryinformation._spfa;
  std::cerr << "my _spfa " << _spfa.rowstart() << "\n";
  std::cerr << "my _spfa " << _spfa.col() << "\n";
  std::cerr << "your spfa " << spfa.rowstart()<< "\n";
  std::cerr << "your spfa " << spfa.col()<< "\n";
  _spfa = spfa;
}

std::string CurvedBoundaryInformation::getClassName() const
{
  return "CurvedBoundaryInformation";
}

/*--------------------------------------------------------------------------*/
alat::SparsityPatternFixArray<2>& CurvedBoundaryInformation::curvedCellIndices()
{
  return _spfa;
}

const alat::Map<int, FadalightMesh::CurvedBoundaryDescriptionInterface*>& CurvedBoundaryInformation::get() const
{
  return _curvedboundaries;
}

const FadalightMesh::CurvedBoundaryDescriptionInterface* CurvedBoundaryInformation::get(int color) const
{
  return _curvedboundaries[color];
}

FadalightMesh::CurvedBoundaryDescriptionInterface* CurvedBoundaryInformation::get(int color)
{
  return _curvedboundaries[color];
}

FadalightMesh::CurvedBoundaryDescriptionInterface*& CurvedBoundaryInformation::getPointer(int color)
{
  return _curvedboundaries[color];
}

bool CurvedBoundaryInformation::boundaryColorIsCurved(int color) const
{
  return _curvedboundaries.find(color) != _curvedboundaries.end();
}

const alat::SparsityPatternFixArray<2>& CurvedBoundaryInformation::curvedCellIndices() const
{
  return _spfa;
}

// bool CurvedBoundaryInformation::cellIsCurved(int iK) const
// {
//   if( iK >= _celliscurved.size() )
//   {
//     _error_string("cellIsCurved", "iK,_celliscurved.size() =", iK, _celliscurved.size());
//   }
//   return _celliscurved[iK] != -1;
// }

int CurvedBoundaryInformation::rowstart(int iK) const
{
  if(_celliscurved[iK] == -1)
  {
    return 0;
  }
  return _spfa.rowstart(_celliscurved[iK]);
}

int CurvedBoundaryInformation::rowstop(int iK) const
{
  if(_celliscurved[iK] == -1)
  {
    return 0;
  }
  return _spfa.rowstop(_celliscurved[iK]);
}

const alat::FixArray<2, int>& CurvedBoundaryInformation::getCurvedInfo(int pos) const
{
  return _spfa.col(pos);
}

/*--------------------------------------------------------------------------*/
void CurvedBoundaryInformation::constructBoundaryInformation(const FadalightMesh::MeshInterface* mesh)
{
  _celliscurved.clear();
  _celliscurved.set_size(mesh->getNCells());
  _celliscurved.fill(-1);

  alat::IntMap  iK_to_count;
  const FadalightMesh::BoundaryInfo* boundaryinfo = mesh->getBoundaryInfo();
  assert(boundaryinfo);
  int count = 0;
  {
    const alat::armaivec& colors = boundaryinfo->getColors();
    for(int icol = 0; icol < colors.size(); icol++)
    {
      int color = colors[icol];
      if( _curvedboundaries.find(color) == _curvedboundaries.end() )
      {
        continue;
      }
      const alat::armaivec& cells_of_color = boundaryinfo->getCellsOfColor(color);
      for(int i = 0; i < cells_of_color.size(); i++)
      {
        int iK = cells_of_color[i];
        if( iK_to_count.find(iK) == iK_to_count.end() )
        {
          iK_to_count[iK] = count++;
        }
      }
    }
  }
  const FadalightMesh::CurvedInteriorSideInfo* curvedinteriorsideinfo = mesh->getCurvedInteriorSideInfo();
  {
    const alat::armaivec& colors = curvedinteriorsideinfo->getColors();
    for(int icol = 0; icol < colors.size(); icol++)
    {
      int color = colors[icol];
      if( _curvedboundaries.find(color) == _curvedboundaries.end() )
      {
        continue;
      }
      const alat::armaivec& sides_of_color = curvedinteriorsideinfo->getSidesOfColor(color);
      for(int i = 0; i < sides_of_color.size(); i++)
      {
        int iside = sides_of_color[i];
        int iKL = mesh->getCellIdOfSide(iside, 0);
        int iKR = mesh->getCellIdOfSide(iside, 1);
        if( iK_to_count.find(iKL) == iK_to_count.end() )
        {
          iK_to_count[iKL] = count++;
        }
        if( iK_to_count.find(iKR) == iK_to_count.end() )
        {
          iK_to_count[iKR] = count++;
        }
      }
    }
  }


  alat::SparsityPatternFixArraySoft<2> sparsitypatternfixarraysoft(count);
  {
    const alat::armaivec& colors = boundaryinfo->getColors();
    for(int icol = 0; icol < colors.size(); icol++)
    {
      int color = colors[icol];
      if( _curvedboundaries.find(color) == _curvedboundaries.end() )
      {
        continue;
      }
      const alat::armaivec& cells_of_color = boundaryinfo->getCellsOfColor(color);
      const alat::armaivec& localside_of_color = boundaryinfo->getSidesIdOfCellsOfColor(color);
      for(int i = 0; i < cells_of_color.size(); i++)
      {
        int iK = cells_of_color[i];
        alat::FixArray<2, int> f;
        f[0] = color;
        f[1] = localside_of_color[i];
        int count = iK_to_count[iK];
        sparsitypatternfixarraysoft[count].insert(f);
        _celliscurved[iK] = count;
      }
    }
  }
  {
    const alat::armaivec& colors = curvedinteriorsideinfo->getColors();
    for(int icol = 0; icol < colors.size(); icol++)
    {
      int color = colors[icol];
      bool found = _curvedboundaries.find(color) != _curvedboundaries.end();
      if( _curvedboundaries.find(color) == _curvedboundaries.end() )
      {
        continue;
      }
      const alat::armaivec& sides_of_color = curvedinteriorsideinfo->getSidesOfColor(color);
      // std::cerr << "sides_of_color="<<sides_of_color.size() << "\n";
      for(int i = 0; i < sides_of_color.size(); i++)
      {
        int iside = sides_of_color[i];
        int iKL = mesh->getCellIdOfSide(iside, 0);
        int iKR = mesh->getCellIdOfSide(iside, 1);
        // std::cerr << "iside="<< iside<< "iKL="<<iKL << " iKR="<<iKR<<"\n";
        {
          alat::FixArray<2, int> f;
          f[0] = color;
          f[1] = mesh->getLocalIndexOfSideInCell(iKL, iside);
          int count = iK_to_count[iKL];
          sparsitypatternfixarraysoft[count].insert(f);
          _celliscurved[iKL] = count;
        }
        {
          alat::FixArray<2, int> f;
          f[0] = color;
          f[1] = mesh->getLocalIndexOfSideInCell(iKR, iside);
          int count = iK_to_count[iKR];
          sparsitypatternfixarraysoft[count].insert(f);
          _celliscurved[iKR] = count;
        }
      }
    }
  }
  // std::cerr << "iK_to_count="<<iK_to_count<<"\n";
  // std::cerr << "_celliscurved="<<_celliscurved<<"\n";
  // assert(0);


  _spfa.set_size(sparsitypatternfixarraysoft);
}

/*---------------------------------------------------*/
void CurvedBoundaryInformation::readCurvedBoundaryDescription(std::istream& in)
{
  int n = 0;
  streampos nodestart;
  string toto;
  while( !in.eof() )
  {
    nodestart = in.tellg();
    getline(in, toto);
    // std::cerr << "toto="<<toto<<"\n";
    alat::StringVector totos = alat::Tokenize(toto, " ");
    if(totos.size() != 2)
    {
      continue;
    }
    if(totos[0] == "CurvedSides")
    {
      n = atoi( totos[1].c_str() );
      break;
    }
  }
  if( !in.eof() )
  {
    in.seekg(nodestart);
    // std::cerr << "file.eof()="<<file.eof()<<"\n";
    getline(in, toto);
    // std::cerr << "first="<<toto<<"\n";
  }

  // std::cerr << "CurvedBoundaryInformation n="<<n<<"\n";
  // assert(0);
  // std::cerr << "CurvedBoundaryInformation::read() n = " << n << "\n";
  CurvedBoundaryDescriptionConstructor curvedboundarydescriptionconstructor;
  alat::StringVector names(n);
  alat::armaivec colors(n);
  std::string name;
  int color;
  for(int i = 0; i < n; i++)
  {
    in >> color >> name;
    // std::cerr << "color="<<color<<" name"<<name<<"\n";
    _curvedboundaries[color] = curvedboundarydescriptionconstructor.newDescription(name);
    _curvedboundaries[color]->read(in);
    // file >> colors[i] >> names[i];
    // std::cerr << "CurvedBoundaryInformation::read() i = " << i << " " << colors[i] << " " << names[i] << "\n";
  }
  // CurvedBoundaryDescriptionConstructor BDC(*this, names, colors);
  // for(iterator p = _curvedboundaries.begin(); p != _curvedboundaries.end(); p++)
  // {
  //   FadalightMesh::CurvedBoundaryDescriptionInterface* BDI = p->second;
  //   BDI->read(file);
  //   // std::cerr << "CurvedBoundaryInformation::read() p = " << p->first << " " << p->second->getParameters() << "\n";
  // }
  // write("toto","ascii");
  // assert(0);
}

/*---------------------------------------------------*/

void CurvedBoundaryInformation::writeCurvedBoundaryDescription(std::ostream& out) const
{
  int n = _curvedboundaries.size();
  if(n > 0)
  {
    out << "CurvedSides " << n << "\n";
    for(const_iterator p = _curvedboundaries.begin(); p != _curvedboundaries.end(); p++)
    {
      out << p->first << " " << p->second->getClassName() << "\n";
      p->second->write(out);
    }
    // for(const_iterator p = _curvedboundaries.begin(); p != _curvedboundaries.end(); p++)
    // {
    //   p->second->write(file);
    // }
  }
}

/*---------------------------------------------------*/
void CurvedBoundaryInformation::load(std::string filename)
{
  std::ifstream file( filename.c_str() );
  readCurvedBoundaryDescription(file);
  _celliscurved.load(file);
  curvedCellIndices().load(file);
  file.close();
}

/*---------------------------------------------------*/
void CurvedBoundaryInformation::save(std::string filename, arma::file_type datatype) const
{
  // std::cerr << "CurvedBoundaryInformation::write() _curvedboundaries.size() = " << _curvedboundaries.size() << "\n";
  std::ofstream file( filename.c_str() );
  file << "CurvedSides "<< _curvedboundaries.size() << "\n";
  for(const_iterator p = _curvedboundaries.begin(); p != _curvedboundaries.end(); p++)
  {
    file << p->first << " " << p->second->getClassName() << "\n";
    p->second->write(file, datatype);
  }
  // for(const_iterator p = _curvedboundaries.begin(); p != _curvedboundaries.end(); p++)
  // {
  //   p->second->write(file, datatype);
  // }
  file << "\n";
  _celliscurved.save(file, datatype);
  curvedCellIndices().save(file, datatype);
  file.close();
}

/*---------------------------------------------------*/
void CurvedBoundaryInformation::set_size(const CurvedBoundaryInformation& curvedboundaryinformation)
{
  // const CurvedBoundaryInformation* curvedboundaryinformation = dynamic_cast<const CurvedBoundaryInformation*>(&curvedboundaryinformationinterface);
  // assert(curvedboundaryinformation);
  *this = curvedboundaryinformation;
}
