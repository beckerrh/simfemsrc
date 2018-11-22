#include  "FadalightMesh/fadalightmeshbase3d.hpp"
#include  <algorithm>
#include  <fstream>

using namespace FadalightMesh;
using namespace std;

template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
CellBase<2> FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::_dummyedge;


/*---------------------------------------------------------*/

template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::FadalightMeshBase3d() : FadalightMeshBase<3, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>( )
{}
template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
const typename FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::Edge& FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::_getEdgeOfCell(int i, int ii) const
{
  FadalightMesh::MeshInterface::_notWritten("_getEdgeOfCell");
  return _dummyedge;
}

template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
const typename FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::Edge& FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::_getEdgeOfSide(int i, int ii) const
{
  FadalightMesh::MeshInterface::_notWritten("_getEdgeOfSide");
  return _dummyedge;
}

template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
alat::Vector<typename FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::Edge>& FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::getEdges()
{
  return _edges;
}

template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
alat::Vector<typename FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::EdgeCell>& FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::getEdgesOfCells()
{
  return _edges_of_cells;
}

template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
alat::Vector<typename FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::EdgeSide>&  FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::getEdgesOfSides()
{
  return _edges_of_sides;
}

template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
const alat::Vector<typename FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::Edge>& FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::getEdges() const
{
  return _edges;
}

template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
const alat::Vector<typename FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::EdgeCell>& FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::getEdgesOfCells() const
{
  return _edges_of_cells;
}

template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
const alat::Vector<typename FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::EdgeSide>& FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::getEdgesOfSides() const
{
  return _edges_of_sides;
}

template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
const typename FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::Edge& FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::getEdge(int i) const
{
  return _edges[i];
}

template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
const typename FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::EdgeCell& FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::getEdgesOfCell(int i) const
{
  return _edges_of_cells[i];
}

template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
const typename FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::EdgeSide& FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::getEdgesOfSides(int i) const
{
  return _edges_of_sides[i];
}

template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
std::string FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::getClassName() const
{
  return "FadalightMesh::FadalightMeshBase3d";
}

template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
int FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::getNEdgesPerCell(int i) const
{
  return EDGESPERCELL;
}

template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
int FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::getNEdgesPerSide(int i) const
{
  return NODESPERSIDE;
}

template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
int FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::getNEdges() const
{
  return _edges.size();
}

template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
int FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::getNodeIdOfEdge(int i, int ii) const
{
  return _edges[i][ii];
}

template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
int FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::getEdgeIdOfCell(int i, int ii) const
{
  return _edges_of_cells[i][ii];
}

template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
int FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::getEdgeIdOfSide(int i, int ii) const
{
  return _edges_of_sides[i][ii];
}
//
// template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
// std::string FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::getEnsightType() const
// {
//   return "hexa8";
// }

/*-------------------------------------------------------*/

template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
alat::Node FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::getNodeOfEdge(int ie) const
{
  alat::Node v;
  double d = 0.5;
  for(int ii = 0; ii < 2; ii++)
  {
    int in = getNodeIdOfEdge(ie, ii);
    const alat::Node& vii = FadalightMeshBase<3, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNode(in);
    v.x() += d*vii.x();
    v.y() += d*vii.y();
    v.z() += d*vii.z();
  }
  return v;
}

/*---------------------------------------------------------*/
template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
int FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::getVtkType() const
{
  if(NODESPERCELL == 8)
  {
    return 12;
  }
  else
  {
    assert(0);
    return 9;
  }
}

/*---------------------------------------------------------*/

template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
void FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::readFadalightMesh(const std::string& basefilename)
{
  FadalightMeshBase<3, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::readFadalightMesh(basefilename);

  string dirname = basefilename + ".fadalightmesh";

  string filename = dirname;
  std::ifstream file;
  // //! read decription file
  //  string filename = dirname + "/DESC";
  //  std::ifstream file(filename.c_str());
  //  assert(file.is_open());
  //  // file << getClassName() << " " << type << endl;
  //  // file << getInfo() << endl;
  //  // file << "nodes cells sides sides_of_cells boundaryinfo";
  //  file.close();

  //! read edges
  filename = dirname + "/edges";
  file.open( filename.c_str() );
  assert( file.is_open() );
  getEdges().load(file);
  file.close();

  //! read EdgesCells
  filename = dirname + "/edges_of_cells";
  file.open( filename.c_str() );
  assert( file.is_open() );
  getEdgesOfCells().load(file);
  file.close();

  //! read EdgesSides
  filename = dirname + "/edges_of_sides";
  file.open( filename.c_str() );
  assert( file.is_open() );
  getEdgesOfSides().load(file);
  file.close();
}

/*---------------------------------------------------------*/

template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
void FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::writeFadalightMesh(const std::string& basefilename, arma::file_type datatype) const
{
  FadalightMeshBase<3, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::writeFadalightMesh(basefilename, datatype);

  string filename;
  std::ofstream file;
  string dirname = basefilename + ".fadalightmesh";
  // //! write decription file
  // string filename = dirname + "/DESC";
  // std::ofstream file(filename.c_str());
  // assert(file.is_open());
  // file << getClassName() << " " << datatype << endl;
  // file << getInfo() << endl;
  // file << "nodes cells sides sides_of_cells boundaryinfo";
  // file.close();
  //
  //! write edges
  filename = dirname + "/edges";
  file.open( filename.c_str() );
  assert( file.is_open() );
  getEdges().save(file, datatype);
  file.close();

  //! write edgescells
  filename = dirname + "/edges_of_cells";
  file.open( filename.c_str() );
  assert( file.is_open() );
  getEdgesOfCells().save(file, datatype);
  file.close();

  //! write edgessides
  filename = dirname + "/edges_of_sides";
  file.open( filename.c_str() );
  assert( file.is_open() );
  getEdgesOfSides().save(file, datatype);
  file.close();
}

/*---------------------------------------------------------*/

template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
void FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>::constructSidesFromCells(BoundarySideToColor& bsides, int color_default)
{
  FadalightMeshBase<3, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::constructSidesFromCells(bsides, color_default);

  // Construction des 'edges'
  std::map<Edge, int> edgetoid;
  int count = 0;
  for(int i = 0; i < FadalightMeshBase<3, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNCells(); i++)
  {
    for(int ii = 0; ii < EDGESPERCELL; ii++)
    {
      Edge e = _getEdgeOfCell(i, ii);
      sort( e.begin(), e.end() );
      typename std::map<Edge, int>::iterator p = edgetoid.find(e);
      if( p == edgetoid.end() )
      {
        edgetoid.insert( make_pair(e, count++) );
      }
    }
  }
  _edges.set_size(count);
  for(typename std::map<Edge, int>::const_iterator p = edgetoid.begin(); p != edgetoid.end(); p++)
  {
    _edges[p->second] = p->first;
  }

  _edges_of_cells.set_size( FadalightMeshBase<3, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNCells() );
  for(int i = 0; i < FadalightMeshBase<3, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNCells(); i++)
  {
    for(int ii = 0; ii < EDGESPERCELL; ii++)
    {
      Edge e = _getEdgeOfCell(i, ii);
      Edge esort = e;
      sort( esort.begin(), esort.end() );
      typename std::map<Edge, int>::iterator p = edgetoid.find(esort);
      if( p == edgetoid.end() )
      {
        std::cerr << "edge not found " << e << "\n";
        std::cerr << "edge not found " << esort << "\n";
        assert(0);
      }
      _edges_of_cells[i][ii] = p->second;
    }
  }
  _edges_of_sides.set_size( FadalightMeshBase<3, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNSides() );
  for(int i = 0; i < FadalightMeshBase<3, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNSides(); i++)
  {
    for(int ii = 0; ii < NODESPERSIDE; ii++)
    {
      Edge e = _getEdgeOfSide(i, ii);
      Edge esort = e;
      sort( esort.begin(), esort.end() );
      typename std::map<Edge, int>::iterator p = edgetoid.find(esort);
      if( p == edgetoid.end() )
      {
        std::cerr << "edge not found " << e << "\n";
        std::cerr << "edge not found " << esort << "\n";
        std::cerr << "side " << FadalightMeshBase<3, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getSide(i) << "\n";
        assert(0);
      }
      _edges_of_sides[i][ii] = p->second;
    }
  }
}

/*---------------------------------------------------------*/
/*---------------------------------------------------------*/

// tetrahedral mesh
#define NODESPERCELL 4
#define SIDESPERCELL 4
#define EDGESPERCELL 6
#define NODESPERSIDE 3
template class FadalightMesh::FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>;
#undef NODESPERCELL
#undef SIDESPERCELL
#undef EDGESPERCELL
#undef NODESPERSIDE

// hexahedral mesh
#define NODESPERCELL  8
#define SIDESPERCELL  6
#define EDGESPERCELL 12
#define NODESPERSIDE  4
template class FadalightMesh::FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>;
#undef NODESPERCELL
#undef SIDESPERCELL
#undef EDGESPERCELL
#undef NODESPERSIDE
