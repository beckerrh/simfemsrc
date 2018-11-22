#include  "FadalightMesh/quadrilateralmesh.hpp"
#include  "Alat/stringvector.hpp"
#include  "Alat/tokenize.hpp"
#include  <algorithm>
#include  <fstream>

using namespace FadalightMesh;
using namespace std;

/*--------------------------------------------------------------------------*/
QuadrilateralMesh::~QuadrilateralMesh()
{}

QuadrilateralMesh::QuadrilateralMesh() : FadalightMesh::FadalightMeshBase2d<4>()
{}

QuadrilateralMesh::QuadrilateralMesh(const QuadrilateralMesh& quadrilateralmesh) : FadalightMesh::FadalightMeshBase2d<4>(quadrilateralmesh)
{}

std::string QuadrilateralMesh::getClassName() const
{
  return "FadalightMesh::QuadrilateralMesh";
}

std::string QuadrilateralMesh::getCellType() const
{
  return "Quadrilateral";
}

FadalightMeshEnums::meshtype QuadrilateralMesh::getType() const
{
  return FadalightMeshEnums::QuadrilateralMesh;
}

/*--------------------------------------------------------------------------*/

int QuadrilateralMesh::getCouplingOffset(int iS) const
{
  return -1;
}

/*--------------------------------------------------------------------------*/

int QuadrilateralMesh::getNodeIdOfSideOfCell(int iK, int iis, int ii) const
{
  int S;
  if(iis == 0)
  {
    if(ii == 0)
    {
      S = getNodeIdOfCell(iK, 0);
    }
    else
    {
      S = getNodeIdOfCell(iK, 1);
    }
  }
  else if(iis == 1)
  {
    if(ii == 0)
    {
      S = getNodeIdOfCell(iK, 1);
    }
    else
    {
      S = getNodeIdOfCell(iK, 2);
    }
  }
  else if(iis == 2)
  {
    if(ii == 0)
    {
      S = getNodeIdOfCell(iK, 2);
    }
    else
    {
      S = getNodeIdOfCell(iK, 3);
    }
  }
  else if(iis == 3)
  {
    if(ii == 0)
    {
      S = getNodeIdOfCell(iK, 3);
    }
    else
    {
      S = getNodeIdOfCell(iK, 0);
    }
  }
  return S;
}

/*--------------------------------------------------------------------------*/

void QuadrilateralMesh::writeQuad(std::string filename) const
{
  string name = filename;
  // name += ".quad";

  ofstream file( name.c_str() );
  assert( file.is_open() );

  file << getNNodes() << " ascii\n";
  for(int i = 0; i < getNNodes(); i++)
  {
    file << getNode(i) << "\n";
  }
  file << "\n";
  file << getNCells() << " ascii\n";
  for(int i = 0; i < getNCells(); i++)
  {
    for(int ii = 0; ii < getNNodesPerCell(i); ii++)
    {
      file << getNodeIdOfCell(i, ii) << " ";
    }
    file << "\n";
  }
  file << "\n";

  {
    const FadalightMesh::BoundaryInfo* boundaryinfo = getBoundaryInfo();
    const alat::armaivec& colors = boundaryinfo->getColors();
    int nbdrysides2 = 0;
    for(int i = 0; i < colors.size(); i++)
    {
      int color = colors[i];
      nbdrysides2 += boundaryinfo->getSidesOfColor(color).size();
    }
    int nbdrysides = boundaryinfo->getNSides();
    assert(nbdrysides == nbdrysides2);
    file << nbdrysides << " ascii\n";
    for(int i = 0; i < colors.size(); i++)
    {
      int color = colors[i];
      const alat::armaivec& sides = boundaryinfo->getSidesOfColor(color);
      for(int j = 0; j < sides.size(); j++)
      {
        for(int ii = 0; ii < getNNodesPerSide(0); ii++)
        {
          file << getNodeIdOfSide(sides[j], ii) << " ";
        }
        file << color << std::endl;
      }
    }
  }

  {
    const FadalightMesh::CurvedInteriorSideInfo* curvedinteriorsideinfo = getCurvedInteriorSideInfo();
    const alat::armaivec& colors = curvedinteriorsideinfo->getColors();
    int nintcurvedsides = curvedinteriorsideinfo->getNSides();
    file << "CurvedInteriorSides " << nintcurvedsides << "\n";
    for(int i = 0; i < colors.size(); i++)
    {
      int color = colors[i];
      const alat::armaivec& sides = curvedinteriorsideinfo->getSidesOfColor(color);
      for(int j = 0; j < sides.size(); j++)
      {
        for(int ii = 0; ii < getNNodesPerSide(0); ii++)
        {
          file << getNodeIdOfSide(sides[j], ii) << " ";
        }
        file << color << std::endl;
      }
    }
  }

  const CurvedBoundaryInformation* BD = getCurvedBoundaryInformation();
  assert(BD);
  BD->writeCurvedBoundaryDescription(file);
}

/*--------------------------------------------------------------------------*/
void QuadrilateralMesh::readQuad(std::string filename)
{
  string name = filename;
  // name += ".quad";

  ifstream file( name.c_str() );
  if( !file.is_open() )
  {
    std::cerr<<"*** QuadrilateralMesh::readQuad() : cannot read file \""<<filename<<"\"\n";
    assert(0);
  }

  streampos nodestart;
  string toto;
  while(1)
  {
    nodestart = file.tellg();
    getline(file, toto);
    if(toto[0] != '#')
    {
      break;
    }
  }
  file.seekg(nodestart);

  getAllNodes().load(file);

//  double scaling=10.0;
  // double scaling = 1.0;
  // for(int i = 0; i < getNodes().size(); i++)
  // {
  //   alat::Node& n = _nodes[i];
  //   n.scale(scaling);
  // }
//  std::cerr   << "nodes " << getNodes() << "\n";

  // on n'utilse pas ma fonction read pour les quads, car on veut lire que les nodeids !
  alat::Vector<Quadrilateral>& cells = getCells();
  int nc;
  std::string datatype;
  file>>nc>>datatype;
  cells.set_size(nc);
  for(int i = 0; i < nc; i++)
  {
    for(int ii = 0; ii < 4; ii++)
    {
      file>>cells[i][ii];
    }
  }
  // std::cerr << "cells " << cells << "\n";

  // lecture infos bord
  QuadrilateralMesh::BoundarySideToColor bstc;
  QuadrilateralMesh::Side S;
  int nb;
  file>>nb>>datatype;
  for(int i = 0; i < nb; i++)
  {
    for(int ii = 0; ii < 2; ii++)
    {
      file>>S[ii];
    }
    int color;
    file>>color;
    sort( S.begin(), S.end() );
    bstc[S] = color;
  }

  streampos afterboundaries = file.tellg();

  // lecture infos interior curved sides
  QuadrilateralMesh::BoundarySideToColor istc;
  int nistc = 0;
  while( !file.eof() )
  {
    nodestart = file.tellg();
    getline(file, toto);
    alat::StringVector totos = alat::Tokenize(toto, " ");
    if(totos.size() != 2)
    {
      continue;
    }
    if(totos[0] == "CurvedInteriorSides")
    {
      nistc = std::atoi( totos[1].c_str() );
      break;
    }
  }
  // std::cerr << "QuadrilateralMesh::readQuad() filename="<<filename << " nistc="<<nistc<<"\n";
  
  if( !file.eof() )
  {
    file.seekg(nodestart);
    getline(file, toto);
    // std::cerr << "first="<<toto<<"\n";
  }
  for(int i = 0; i < nistc; i++)
  {
    for(int ii = 0; ii < 2; ii++)
    {
      file>>S[ii];
    }
    int color;
    file>>color;
    sort( S.begin(), S.end() );
    istc[S] = color;
  }


  CurvedBoundaryInformation* curvedboundaryinformation = getCurvedBoundaryInformation();
  assert(curvedboundaryinformation);

  file.clear();
  file.seekg(afterboundaries);

  curvedboundaryinformation->readCurvedBoundaryDescription(file);
  file.close();

  constructSidesFromCells(bstc, istc);

  curvedboundaryinformation->constructBoundaryInformation(this);
}

/*--------------------------------------------------------------------------*/

void QuadrilateralMesh::check() const
{
  for(int ic = 0; ic < getNCells(); ic++)
  {
    for(int ii = 0; ii < getNSidesPerCell(ic); ii++)
    {
      int is = getSideIdOfCell(ic, ii);
      if( ( getCellIdOfSide(is, 0) != ic )&&( getCellIdOfSide(is, 1) != ic ) )
      {
        std::cerr << "*** ERROR in QuadrilateralMesh::check(): cell not found in side "<< ic << " " <<is << "\n";
        std::cerr << "side: " << getNodeIdOfSide(is, 0) << " " << getNodeIdOfSide(is, 1) << "\n";
        std::cerr << "cells of side: " << getCellIdOfSide(is, 0) << " " << getCellIdOfSide(is, 1) << "\n";
        assert(0);
        exit(1);
      }
    }
  }
}

/*--------------------------------------------------------------------------*/

const QuadrilateralMesh::Side& QuadrilateralMesh::_getSideOfCell(int iK, int iis) const
{
  if(iis == 0)
  {
    _S[0] = getNodeIdOfCell(iK, 0);
    _S[1] = getNodeIdOfCell(iK, 1);
  }
  else if(iis == 1)
  {
    _S[0] = getNodeIdOfCell(iK, 1);
    _S[1] = getNodeIdOfCell(iK, 2);
  }
  else if(iis == 2)
  {
    _S[0] = getNodeIdOfCell(iK, 2);
    _S[1] = getNodeIdOfCell(iK, 3);
  }
  else if(iis == 3)
  {
    _S[0] = getNodeIdOfCell(iK, 3);
    _S[1] = getNodeIdOfCell(iK, 0);
  }
  return _S;
}

// /*--------------------------------------------------------------------------*/
//
// void QuadrilateralMesh::getMeshSizeForStabilization(double& hs, int iS, int iK, int iil) const
// {
//   const QuadrilateralMesh::Side& side = _getSideOfCell(iK, iil);
//   Node N = getNode(side[0]);
//   N -= getNode(side[1]);
//   hs = N.norm();
// }

/*--------------------------------------------------------------------------*/

// alat::Vector<QuadrilateralMesh::Quadrilateral>& QuadrilateralMesh::getQuadrilaterals()
// {
//   return getCells();
// }
// 
// const alat::Vector<QuadrilateralMesh::Quadrilateral>& QuadrilateralMesh::getQuadrilaterals() const
// {
//   return getCells();
// }

// const QuadrilateralMesh::Quadrilateral& QuadrilateralMesh::getQuadrilateral(int i) const
// {
//   return getCell(i);
// }

// const QuadrilateralMesh::Edge& QuadrilateralMesh::getEdge(int i) const
// {
//   return getSide(i);
// }

void QuadrilateralMesh::read(std::string filename)
{
  // std::cerr << "QuadrilateralMesh::read() filename " << filename << "\n";
  // assert(0);
  // readFadalightMesh(filename);
  checkGeoFile(filename);
  readQuad(filename);
}

void QuadrilateralMesh::writeSimpleMesh(std::string filename, std::string type) const
{
  writeQuad(filename);
}

void QuadrilateralMesh::write(std::string filename) const
{
  // std::cerr<< "QuadrilateralMesh::write() filename="<<filename<<"\n";
  // assert(0);
  writeQuad(filename);
}

void QuadrilateralMesh::getLocalIndicesOfSidesInCell(alat::armaivec& sideindex_a, alat::armaivec& sideindex_e) const
{
  sideindex_a.set_size(4);
  sideindex_e.set_size(4);
  sideindex_a[0] = 0;
  sideindex_e[0] = 1;
  sideindex_a[1] = 1;
  sideindex_e[1] = 2;
  sideindex_a[2] = 2;
  sideindex_e[2] = 3;
  sideindex_a[3] = 3;
  sideindex_e[3] = 0;
}

void QuadrilateralMesh::getLocalIndicesOfSidesAndDiagonalsInCell(alat::armaivec& sideindex_a, alat::armaivec& sideindex_e) const
{
  sideindex_a.set_size(6);
  sideindex_e.set_size(6);
  sideindex_a[0] = 0;
  sideindex_e[0] = 1;
  sideindex_a[1] = 1;
  sideindex_e[1] = 2;
  sideindex_a[2] = 2;
  sideindex_e[2] = 3;
  sideindex_a[3] = 3;
  sideindex_e[3] = 0;
  sideindex_a[4] = 0;
  sideindex_e[4] = 2;
  sideindex_a[5] = 1;
  sideindex_e[5] = 3;
}
