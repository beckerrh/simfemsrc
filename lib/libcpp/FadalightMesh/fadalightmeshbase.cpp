#include  "Alat/directoryandfiles.hpp"
#include  "Alat/sparsitypatternsoft.hpp"
#include  "FadalightMesh/coarseninfo.hpp"
#include  "FadalightMesh/fadalightmeshbase.hpp"
#include  "FadalightMesh/geometryobjectsconstructor.hpp"
#include  "FadalightMesh/hangingnodeinfo.hpp"
#include  "FadalightMesh/hangingsideinfo.hpp"
#include  "FadalightMesh/refineinfo.hpp"
#include  "Alat/stringvector.hpp"
#include  "Alat/tokenize.hpp"
// #include  "H5Cpp.hpp"
#include  <algorithm>
#include  <fstream>
#include  <stdio.h>
#include  <string.h>

using namespace FadalightMesh;
using namespace std;

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
CellBase<NODESPERSIDE> FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::_dummyside;

/*---------------------------------------------------------*/

template<int N>
int CellBase<N>::node(int i) const
{
  assert(i < N);
  return alat::FixArray<N, int>::operator[](i);
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::~FadalightMeshBase()
{
  for(std::map<std::string, FadalightMesh::GeometryObject*>::iterator p = _geometryobjects.begin(); p != _geometryobjects.end(); p++)
  {
    if(p->second)
    {
      delete p->second;
      p->second = NULL;
    }
  }
  delete _geometryobjects_constructor;
  _geometryobjects_constructor = NULL;
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::FadalightMeshBase() : FadalightMesh::MeshInterface(), _infilename("none"), _visutype("cg")
{
  _geometryobjects_constructor = new GeometryObjectsConstructor;
  createGeometryObject("HangingNodeInfo");
  createGeometryObject("HangingSideInfo");
  // createGeometryObject("CurvedBoundaryInformation");
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::FadalightMeshBase(const FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>& fadalightmeshbase) : FadalightMesh::MeshInterface(fadalightmeshbase), _infilename("none")
{
  assert(0);
}


/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
void FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::setResolution(int level) const
{}


template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
int FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getResolution() const
{
  return 0;
}

/*--------------------------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
bool FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::isMultilevel() const
{
  return false;
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
void FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNodesOfCell(int iK, alat::Vector<alat::Node>& F) const
{
  F.set_size(NODESPERCELL);
  for(int ii = 0; ii < NODESPERCELL; ii++)
  {
    F[ii] = getNode( getNodeIdOfCell(iK, ii) );
  }
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
bool FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::geometryObjectExists(std::string name) const
{
  bool exists = _geometryobjects.find(name) != _geometryobjects.end();
  return exists;
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
const FadalightMesh::GeometryObject* FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getGeometryObject(std::string name) const
{
  for(std::map<std::string, FadalightMesh::GeometryObject*>::const_iterator p = _geometryobjects.begin(); p != _geometryobjects.end(); p++)
  {
    if(p->first == name)
    {
      return p->second;
    }
  }
  std::cerr << "*** FadalightMesh::FadalightMeshBase::getGeometryObject() : object \""<<name<<"\" not found\n";
  std::cerr << _geometryobjects.size() << "\n";
  for(std::map<std::string, FadalightMesh::GeometryObject*>::const_iterator p = _geometryobjects.begin(); p != _geometryobjects.end(); p++)
  {
    std::cerr << p->first << "\n";
  }
  return NULL;
  //      assert(0);
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
FadalightMesh::GeometryObject* FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getGeometryObject(std::string name)
{
  if( _geometryobjects.find(name) == _geometryobjects.end() )
  {
    std::cerr << "*** FadalightMesh::FadalightMeshBase::getGeometryObject() : object \""<<name<<"\" not found\n";
    assert(0);
    exit(1);
  }
  return _geometryobjects[name];
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
alat::Vector<alat::Node>& FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getAllNodes()
{
  return _nodes;
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
alat::Vector<typename FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::Cell>& FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getCells()
{
  return _cells;
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
alat::Vector<typename FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::Side>& FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getSides()
{
  return _sides;
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
alat::Vector<typename FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::SideCell>& FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getSidesOfCells()
{
  return _sides_of_cells;
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
alat::Vector<typename FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::CellSide>&  FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getCellsOfSides()
{
  return _cells_of_sides;
}

/*---------------------------------------------------------*/
//
// template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
// CurvedBoundaryInformation* FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getCurvedBoundaryInformation()
// {
//   FadalightMesh::GeometryObject* GO = getGeometryObject("CurvedBoundaryInformation");
//   CurvedBoundaryInformation* MO = dynamic_cast<CurvedBoundaryInformation*>( GO );
//   assert(MO);
//   return MO;
// }
//
// /*---------------------------------------------------------*/
//
// template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
// const CurvedBoundaryInformation* FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getCurvedBoundaryInformation() const
// {
//   const FadalightMesh::GeometryObject* GO = getGeometryObject("CurvedBoundaryInformation");
//   const CurvedBoundaryInformation* MO = dynamic_cast<const CurvedBoundaryInformation*>( GO );
//   assert(MO);
//   return MO;
// }

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
const alat::Vector<alat::Node>&   FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNodes() const
{
  return _nodes;
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
const alat::Vector<typename FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::Cell>&   FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getCells() const
{
  return _cells;
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
const alat::Vector<typename FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::Side>&   FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getSides() const
{
  return _sides;
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
const alat::Vector<typename FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::SideCell>&  FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getSidesOfCells() const
{
  std::cerr<<"getSidesOfCells "<<DIM<<" "<<NODESPERCELL<<" "<<SIDESPERCELL<<" "<<NODESPERSIDE<<'\n';
  return _sides_of_cells;
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
const alat::Vector<typename FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::CellSide>&  FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getCellsOfSides() const
{
  return _cells_of_sides;
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
const typename FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::Cell& FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getCell(int i) const
{
  return _cells[i];
}


/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
const typename FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::Side& FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getSide(int i) const
{
  return _sides[i];
}


/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
const typename FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::SideCell& FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getSidesOfCell(int i) const
{
  return _sides_of_cells[i];
}

/*---------------------------------------------------------*/
template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
const std::string& FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getInfileName() const
{
  return _infilename;
}

/*---------------------------------------------------------*/
template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
std::string FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getClassName() const
{
  return "FadalightMesh::FadalightMeshBase";
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
std::string FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getInfo() const
{
  std::string info;
  std::stringstream ss1, ss2, ss3;
  ss1 << getNNodes();
  ss2 << getNCells();
  ss3 << getNSides();
  return _infilename + "_" + ss1.str()+"_"+ss2.str()+"_"+ss3.str();
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
int FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getDimension() const
{
  return DIM;
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
int FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNNodesPerCell(int i) const
{
  return NODESPERCELL;
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
int FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNNodesPerSide(int i) const
{
  return NODESPERSIDE;
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
int FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNSidesPerCell(int i) const
{
  return SIDESPERCELL;
}

// /*---------------------------------------------------------*/
//
// template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
// const alat::SparsityPatternFixArray<2>& FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getCellNeighbours() const
// {
//   return _neighbors_of_cells;
// }

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
int FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNNodes() const
{
  return _nodes.size();
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
int FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNCells() const
{
  return _cells.size();
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
int FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNSides() const
{
  return _sides.size();
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
int FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNEdges() const
{
  return 0;
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
int FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNodeIdOfCell(int iK, int ii) const
{
  return _cells[iK][ii];
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
int FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNodeIdOfSide(int iS, int ii) const
{
  return _sides[iS][ii];
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
int FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getSideIdOfCell(int iK, int ii) const
{
  return _sides_of_cells[iK][ii];
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
int FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getCellIdOfSide(int iS, int ii) const
{
  return _cells_of_sides[iS][ii];
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
const alat::Node& FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNode(int i) const
{
  return _nodes[i];
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
const alat::Node& FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNodeOfSide(int is, int ii) const
{
  int i =  getNodeIdOfSide(is, ii);
  return getNode(i);
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
const alat::Node& FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNodeOfCell(int iK, int ii) const
{
  int i =  getNodeIdOfCell(iK, ii);
  assert(i >= 0);
  assert( i < getNNodes() );
  return getNode(i);
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
int FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getLocalIndexOfSideInCell(int iK, int iS) const
{
  int ii = -1;
  for(int jj = 0; jj < getNSidesPerCell(iK); jj++)
  {
    if(getSideIdOfCell(iK, jj) == iS)
    {
      ii = jj;
    }
  }
  assert(ii != -1);
  return ii;
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
const FadalightMesh::BoundaryInfo* FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getBoundaryInfo() const
{
  return &_boundaryinfo;
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
FadalightMesh::BoundaryInfo* FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getBoundaryInfo()
{
  return &_boundaryinfo;
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
const FadalightMesh::CurvedInteriorSideInfo* FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getCurvedInteriorSideInfo() const
{
  return &_curvedinteriorsideinfo;
}
//
// /*---------------------------------------------------------*/
//
// template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
// FadalightMesh::CurvedInteriorSideInfo* FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getCurvedInteriorSideInfo()
// {
//   return &_boundaryinfo;
// }

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
const alat::armaivec& FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getBoundaryColors() const
{
  return getBoundaryInfo()->getColors();
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
const alat::armaivec& FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getBoundarySides(int color) const
{
  return getBoundaryInfo()->getSidesOfColor(color);
}

// /*---------------------------------------------------------*/
//
// template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
// bool FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::cellIsCurved(int iK) const
// {
//   return getCurvedBoundaryInformation()->cellIsCurved(iK);
// }

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
const typename FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::Side& FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::_getSideOfCell(int i, int ii) const
{
  _notWritten("_getSideOfCell");
  return _dummyside;
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
int FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNLevels() const
{
  return 1;
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
int FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::findNeightborHangingCells(int iK, int iS, alat::Node pt)
{
  assert(getDimension() == 2);
  const FadalightMesh::GeometryObject* geo = getGeometryObject("HangingSideInfo");
  const FadalightMesh::HangingSideInfo* hsinfo = dynamic_cast<const FadalightMesh::HangingSideInfo*>( geo );
  assert(hsinfo);
  int iKr = hsinfo->getCellNumber(0);
  int ils = hsinfo->getLocalSide(0);
  int iSr = getSideIdOfCell(iKr, ils);
  int ih = 0;
  while( iKr != iK && iS != iSr && ih < hsinfo->n() )
  {
    iKr = hsinfo->getCellNumber(ih);
    ils = hsinfo->getLocalSide(ih);
    iSr = getSideIdOfCell(iKr, ils);
    ih++;
  }
  assert(iK == iKr && iS == iSr);
  int nbhs = hsinfo->getNumberOfHangingSides(ih-1);
  for(int js = 0; js < nbhs; js++)
  {
    int iSl = hsinfo->getHangingSides(ih-1, js);
    int iKl = getCellIdOfSide(iSl, 0);
    std::cerr<<"FadalightMeshBase::findNeightborHangingCells iKl "<<iKl<<'\n';
    alat::Node extrem1, extrem2;
    extrem1 = getNodeOfSide(iSl, 0);
    extrem2 = getNodeOfSide(iSl, 1);
    extrem1.add(-1.0, pt);
    extrem2.add(-1.0, pt);
    if(extrem1*extrem2 < 0)
    {
      return iKl;
    }
  }
  assert(0);
  return -1;
}

/*---------------------------------------------------------*/
template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
void FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::createGeometryObject(std::string name)
{
  if( geometryObjectExists(name) )
  {
    /// Changement par Robert et Roland 17/1/2011 (on appelle la fonction plusieurs fois, certains GO sont contruits dans le constructeur)
    // std::cerr<<"*** FadalightMeshBase::createGeometryObject(): ::GeometryObject exists already "<<name<<"\n";
    return;

    assert(0);
  }

  // *****
  _geometryobjects_constructor->constructGeometryObject(_geometryobjects, name);
}

/*-------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
alat::Node FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNodeOfCell(int iK) const
{
  alat::Node v;
  double d = 1.0/double(NODESPERCELL);
  for(int ii = 0; ii < NODESPERCELL; ii++)
  {
    const alat::Node& vii = getNodeOfCell(iK, ii);
    v.x() += d*vii.x();
    v.y() += d*vii.y();
    v.z() += d*vii.z();
  }
  return v;
}

/*-------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
alat::Node FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNodeOfSide(int is) const
{
  alat::Node v;
  double d = 1.0/double(NODESPERSIDE);
  for(int ii = 0; ii < NODESPERSIDE; ii++)
  {
    const alat::Node& vii = getNodeOfSide(is, ii);
    v.x() += d*vii.x();
    v.y() += d*vii.y();
    v.z() += d*vii.z();
  }
  return v;
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
void FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::readFadalightMesh(const std::string& basefilename)
{
  _infilename = basefilename;
  string dirname = basefilename+".fadalightmesh";
  // std::cerr << "@@@@@@@@@@@ _infilename " << _infilename << "\n";
  string filename;
  std::ifstream file;

  //! check name !
  filename = dirname+"/name";
  file.open( filename.c_str() );
  if( !file.is_open() )
  {
    std::cerr<<"*** ERROR in FadalightMeshBase::readFadalightMesh() : cannot open file "<<filename<<"\n";
    assert(0);
  }
  std::string name, datatype;
  file>>name>>datatype;
  if( name != getClassName() )
  {
    std::cerr<<"*** ERROR in FadalightMeshBase::readFadalightMesh() : meshname(file) is \""<<name<<"\" but I am \""<<getClassName()<<"\"\n";
    assert(0);
  }
  if( ( datatype != "ascii" )&&( datatype != "binary" ) )
  {
    std::cerr<<"*** FadalightMeshBase::readFadalightMesh() unknown datatype "<<datatype<<"\n";
    assert(0);
  }
  file.close();
  if(datatype == "ascii")
  {
    std::cout<<"FadalightMeshBase::readFadalightMesh() \'" << getClassName()<< "\' reading "<<basefilename<<"\n";//" : "<<datatype<<"\n"
  }

  //! read geometry-objects decription file and geometry objects
  filename = dirname+"/geometry_objects.desc";
  file.open( filename.c_str() );
  if( !file.is_open() )
  {
    std::cerr<<"*** ERROR in FadalightMeshBase::readFadalightMesh() : cannot open file "<<filename<<"\n";
    assert(0);
  }
  int n_geometric_objects;
  file>>n_geometric_objects;
  // std::cout << "FadalightMeshBase::readFadalightMesh() : n_geometric_objects= " << n_geometric_objects << "\n";
  // name;
  for(int i = 0; i < n_geometric_objects; i++)
  {
    file>>name;
    std::cout << "FadalightMeshBase::readFadalightMesh() : creating geometric object: " << name << "\n";
    createGeometryObject(name);
    getGeometryObject(name)->load(dirname+"/"+name);
  }
  file.close();

  //! read nodes
  filename = dirname+"/nodes";
  file.open( filename.c_str() );
  if( !file.is_open() )
  {
    std::cerr<<"*** FadalightMeshBase::readFadalightMesh() : cannot open file "<<filename<<std::endl;
    assert(0);
  }
  _nodes.load(file);
  file.close();
  // std::cerr << "geom\n";
  // std::cerr << "***nodes ok\n";

  //! read cells
  // std::cerr << "===== reading cells\n";
  filename = dirname+"/cells";
  file.open( filename.c_str() );
  assert( file.is_open() );
  _cells.load(file);
  // _measureofcell.read(file);
  file.close();
  // std::cerr << "***cells ok\n";

  //! read sides
  // std::cerr << "===== reading sides\n";
  filename = dirname+"/sides";
  file.open( filename.c_str() );
  assert( file.is_open() );
  _sides.load(file);
  file.close();
  // std::cerr << "***sides ok\n";

  //! read sides_of_cells
  // std::cerr << "===== reading sides_of_cells\n";
  filename = dirname+"/sides_of_cells";
  file.open( filename.c_str() );
  assert( file.is_open() );
  _sides_of_cells.load(file);
  file.close();
  // std::cerr << "***sides_of_cells ok\n";

  //! read cells_of_sides
  // std::cerr << "===== reading cells_of_sides\n";
  filename = dirname+"/cells_of_sides";
  file.open( filename.c_str() );
  assert( file.is_open() );
  _cells_of_sides.load(file);
  file.close();
  // std::cerr << "***cells_of_sides ok\n";

  //! read boundaryinfo
  // std::cerr << "===== reading boundaryinfo\n";
  filename = dirname+"/boundaryinfo";
  _boundaryinfo.read(filename);

  //! read boundaryinfo
  // std::cerr << "===== reading curvedinteriorsideinfo\n";
  // filename = dirname+"/curvedinteriorsideinfo";
  // _curvedinteriorsideinfo.read(filename);
  // std::cerr << "***boundaryinfo ok\n";
  std::cerr << "===== reading finished\n";
}

/*---------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
void FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::writeFadalightMesh(const std::string& basefilename, arma::file_type datatype) const
{
  if(datatype == arma::arma_ascii)
  {
    std::cout<<"FadalightMeshBase::writeFadalightMesh() \'" << getClassName()<< "\' writing "<<basefilename<<"\n";//" : "<<datatype<<"\n"
  }
  //! create directory
  string dirname = basefilename+".fadalightmesh";
  string cmd = "rm -rf "+dirname;
  int error = system( cmd.c_str() );
  assert(!error);
  cmd = "mkdir "+dirname;
  error = system( cmd.c_str() );
  if(error)
  {
    std::cerr<<"*** FadalightMeshBase::writeFadalightMesh(): command failed "<<cmd<<"\n";
    exit(1);
  }
  assert(!error);

  string filename;
  std::ofstream file;

  //! write decription file
  filename = dirname+"/name";
  file.open( filename.c_str() );
  assert( file.is_open() );
  file<<getClassName()<<" ";
  if(datatype == arma::arma_ascii)
  {
    file<<"ascii"<<endl;
  }
  else if(datatype == arma::arma_binary)
  {
    file<<"binary"<<endl;
  }
  else
  {
    // std::cerr << "unknown datatype " << datatype << "\n";
    std::cerr << "unknown datatype "  << "\n";
    exit(1);
  }
  file<<getInfo()<<endl;
  file.close();

  //! write geometry-objects decription file and geometry objects
  filename = dirname+"/geometry_objects.desc";
  file.open( filename.c_str() );
  assert( file.is_open() );
  file<<_geometryobjects.size()<<"\n";
  for(std::map<std::string, FadalightMesh::GeometryObject*>::const_iterator p = _geometryobjects.begin(); p != _geometryobjects.end(); p++)
  {
    std::cerr << "FadalightMeshBase writing " << p->first << " =? " << p->second->getClassName() << "\n";
    file<<p->first<<"\n";
    p->second->save(dirname+"/"+p->second->getClassName(), datatype);
  }
  file.close();

  /*** Temporary process for block structured meshes */
  // filename = dirname+"/CurvedBoundaryInformation";
  // if( !alat::_FileExists(filename) )
  // {
  //   assert(0);
  //   alat::armaivec curvedinfo;
  //   std::cerr<<"getNCells :::"<<getNCells()<<std::endl;
  //   curvedinfo.set_size(getNCells());
  //   curvedinfo.fill(-1);
  //   file.open( filename.c_str() );
  //   file << "0\n\n";
  //   curvedinfo.save(file, datatype);
  //   file << "0 ascii\n";
  //   file << "1 ascii\n";
  //   file << 0;
  //   file.close();
  //
  //   filename = dirname+"/geometry_objects.desc";
  //   file.open( filename.c_str() );
  //   file<<( _geometryobjects.size()+1 )<<"\n";
  //   for(std::map<std::string, FadalightMesh::GeometryObject*>::const_iterator p = _geometryobjects.begin(); p != _geometryobjects.end(); p++)
  //   {
  //     file<<p->first<<"\n";
  //     p->second->save(dirname+"/"+p->second->getClassName(), datatype);
  //   }
  //   file<<"CurvedBoundaryInformation\n";
  //   file.close();
  // }
  /*************************************************/

  //! write nodes
  filename = dirname+"/nodes";
  file.open( filename.c_str() );
  assert( file.is_open() );
  file.precision(12);
  file.setf(ios::scientific);
  _nodes.save(file, datatype);
  file.close();
  //! write cells
  filename = dirname+"/cells";
  file.open( filename.c_str() );
  assert( file.is_open() );
  file.precision(12);
  file.setf(ios::scientific);
  _cells.save(file, datatype);
  // _measureofcell.write(file, datatype);
  file.close();

  //! write sides
  filename = dirname+"/sides";
  file.open( filename.c_str() );
  assert( file.is_open() );
  _sides.save(file, datatype);
  file.close();

  //! write sides_of_cells
  filename = dirname+"/sides_of_cells";
  file.open( filename.c_str() );
  assert( file.is_open() );
  _sides_of_cells.save(file, datatype);
  file.close();

  //! write cells_of_sides
  filename = dirname+"/cells_of_sides";
  file.open( filename.c_str() );
  assert( file.is_open() );
  _cells_of_sides.save(file, datatype);
  file.close();

  //! write boundaryinfo
  filename = dirname+"/boundaryinfo";
  _boundaryinfo.write(filename, datatype);

  //! write curvedinteriorsideinfo
  // filename = dirname+"/curvedinteriorsideinfo";
  // _curvedinteriorsideinfo.write(filename, datatype);
}

/*--------------------------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
void FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::writeSimpleMesh(std::string filename, arma::file_type datatype) const
{
  assert(0);
  string name = filename;
  // name += "."+type;

  ofstream file( name.c_str() );
  if( !file.is_open() )
  {
    std::cerr<<"*** FadalightMeshBase::writeSimpleMesh() : cannot open file \""<<filename<<"\"";
    assert(0);
  }
  getNodes().save(file, datatype);
  getCells().save(file, datatype);
  const FadalightMesh::BoundaryInfo* BI = getBoundaryInfo();
  int nsides = BI->getNSides();
  const alat::armaivec& colors = BI->getColors();
  file<<nsides<<" ascii\n";
  for(int i = 0; i < colors.size(); i++)
  {
    int color = colors[i];
    const alat::armaivec& sides = BI->getSidesOfColor(color);
    for(int j = 0; j < sides.size(); j++)
    {
      for(int ii = 0; ii < NODESPERSIDE; ii++)
      {
        file<<getNodeIdOfSide(sides[j], ii)<<" "<<'\n';
      }
      file<<color<<'\n';
    }
  }
  // const CurvedBoundaryInformation* BD = getCurvedBoundaryInformation();
  // assert(BD);
  // BD->writeCurvedBoundaryDescription(file);
  // file.close();
}

//     assert(disc);
//     int ic=disc->getLocalIdOfInterceptedCell(iK);
//     if (ic>=0) return 0;
//     else return ic+2;
// }

/*---------------------------------------------------------*/
template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
void FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::checkGeoFile(std::string filename) const
{
  ifstream file( filename.c_str() );
  if( !file.is_open() )
  {
    std::cerr<<"*** FadalightMeshBase::checkGeoFile() : cannot read file \""<<filename<<"\"\n";
    assert(0);
  }

}

/*---------------------------------------------------------*/
template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
void FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::constructSidesFromCells(BoundarySideToColor& bstc, const BoundarySideToColor& icsides, int color_default)
{
  constructSidesFromCells(bstc, color_default);

  alat::IntMap size_of_color;
  for(typename BoundarySideToColor::const_iterator p = icsides.begin(); p != icsides.end(); p++)
  {
    int col = p->second;
    if( size_of_color.find(col) == size_of_color.end() )
    {
      size_of_color[col] = 0;
    }
    size_of_color[col]++;
  }
  _curvedinteriorsideinfo.set_size(size_of_color);

  int count=0;
  for(int iside = 0; iside < _sides.size(); iside++)
  {
    Side s = _sides[iside];
    Side ssorted=s;
    sort(ssorted.begin(), ssorted.end());
    for(typename BoundarySideToColor::const_iterator p = icsides.begin(); p != icsides.end(); p++)
    {
      Side s2 = p->first;
      Side s2sorted=s2;
      sort(s2sorted.begin(), s2sorted.end());
      // std::cerr << "ssorted="<< ssorted << " s2sorted="<< s2sorted<< "\n";
      if(ssorted==s2sorted)
      {
        // std::cerr << "FOUND  s="<<s<<" iside="<<iside<<"\n";
        _curvedinteriorsideinfo.getSidesOfColor(p->second)[count++]=iside;
      }
    }
  }
  // assert(0);
}

/*---------------------------------------------------------*/
template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
void FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::constructSidesFromCells(BoundarySideToColor& bstc, int color_default)
{
  //! Assumption: bstc_key is sorted !!
  //!
  //! Assumption: no hanging nodes !!
  //!
  //! BoundarInfo::reInit is done here !


  std::cerr << "constructSidesFromCells() bstc = " << bstc << "\n";

  // helpers
  typedef map<Side, alat::FixArray<2, int> >  SideToInfo;
  SideToInfo _found;

  // first round : how many interor sides ?
  // interior == in two cells
  // boundary == only in one cell
  int ninteriorsides = 0;
  for(int i = 0; i < _cells.size(); i++)
  {
    for(int ii = 0; ii < SIDESPERCELL; ii++)
    {
      Side s = _getSideOfCell(i, ii);
      sort( s.begin(), s.end() );
      typename SideToInfo::iterator p = _found.find(s);
      if( p == _found.end() )
      {
        alat::FixArray<2, int> index;
        _found.insert( make_pair(s, index) );
      }
      else
      {
        ninteriorsides++;
        _found.erase(p);
      }
    }
  }
  // boundary sides
  int nbsides = _found.size();
  // for(typename SideToInfo::iterator p = _found.begin(); p != _found.end(); p++)  nbsides++;
  // cout<<"constructSidesFromCells(): # interior sides, # boundary sides : "<<ninteriorsides<<" "<<nbsides<<endl;

  int nsides = ninteriorsides+nbsides;
  // cout<<"nsides = "<<nsides<<"\n";
  _sides.set_size(nsides);

  // second round : insert data
  _found.clear();
  _sides_of_cells.set_size( getNCells() );
  int count = 0;
  for(int i = 0; i < _cells.size(); i++)
  {
    for(int ii = 0; ii < SIDESPERCELL; ii++)
    {
      // The side plus a sorted copy
      // The sorted copy is necessary in order to retrieve the information
      Side s = _getSideOfCell(i, ii);
      Side ssort = s;
      sort( ssort.begin(), ssort.end() );
      typename SideToInfo::iterator p = _found.find(ssort);
      if( p == _found.end() )
      {
        alat::FixArray<2, int> index;
        index[0] = i;
        index[1] = ii;
        _found.insert( make_pair(ssort, index) );
      }
      else
      {
        // we put in the NON-SORTED !
        // cerr << nsides << " ?  " << count << "\n";
        // std::cerr << "side found " << s << "\n";
        _sides[count] = s;
        int k = p->second[0];
        int kk = p->second[1];
        // cerr << getNCells() << " ?>  " << k << " " << i << " -- " << kk << " " << ii << "\n";
        _sides_of_cells[k][kk] = count;
        _sides_of_cells[i][ii] = count;
        _found.erase(p);
        count++;
      }
    }
  }
  // for(typename SideToInfo::const_iterator p = _found.begin(); p != _found.end(); p++)  cerr << p->first << " -- > " << p->second << "\n";
  //
  // for(typename BoundarySideToColor::iterator pp = bstc.begin(); pp != bstc.end(); pp++)  cerr << " , " << pp->first << " -- > " << pp->second << "\n";

  // now on "_found" are only sides which have been found ounce
  // We suppose the read mesh does not contain "hanging nodes", so the only-ounce found sides
  // are the boundary sides !

  // we add all boundary-sides, which are not yet defined, by giving it the color 0 !!

  // std::cout<<"bstc.size()="<<bstc.size()<<"\n";
  // std::cout<<"_found.size()="<<_found.size()<<"\n";
  int nadditionalboundarysides = 0;
  for(int i = 0; i < _cells.size(); i++)
  {
    for(int ii = 0; ii < SIDESPERCELL; ii++)
    {
      Side s = _getSideOfCell(i, ii);
      Side ssort = s;
      sort( ssort.begin(), ssort.end() );
      typename SideToInfo::iterator p = _found.find(ssort);
      if( p != _found.end() )
      {
        typename map<Side, int>::const_iterator pb = bstc.find(ssort);
        if( pb == bstc.end() )
        {
          cerr<<"***constructSidesFromCells*** not found "<<ssort<<"\n";
          bstc[ssort] = color_default;
          nadditionalboundarysides++;
        }
      }
    }
  }
  // cout<<"constructSidesFromCells(): # boundary sides without color (given color " <<  color_default << " ): "<<nadditionalboundarysides<<endl;

  std::map<int, int> size_of_color;
  for(typename BoundarySideToColor::const_iterator p = bstc.begin(); p != bstc.end(); p++)
  {
    int col = p->second;
    if( size_of_color.find(col) == size_of_color.end() )
    {
      size_of_color[col] = 0;
    }
    size_of_color[col]++;
  }
  FadalightMesh::BoundaryInfo* BI = getBoundaryInfo();
  BI->set_size(size_of_color);
  // boundary sides
  map<int, int> col2size_bsides, col2size_bcells, col2size_bsideids;
  const alat::armaivec& colors = _boundaryinfo.getColors();
  for(int ii = 0; ii < colors.size(); ii++)
  {
    int color = colors[ii];
    col2size_bsides[color] = 0;
    col2size_bcells[color] = 0;
    col2size_bsideids[color] = 0;
  }

  for(int i = 0; i < _cells.size(); i++)
  {
    for(int ii = 0; ii < SIDESPERCELL; ii++)
    {
      // The side plus a sorted copy
      // The sorted copy is necessary in order to retrieve the information
      Side s = _getSideOfCell(i, ii);
      Side ssort = s;
      sort( ssort.begin(), ssort.end() );
      typename SideToInfo::iterator p = _found.find(ssort);
      if( p != _found.end() )
      {
        typename map<Side, int>::const_iterator pb = bstc.find(p->first);
        assert( pb != bstc.end() );
        int color = pb->second;
        alat::armaivec& bsides = _boundaryinfo.getSidesOfColor(color);
        alat::armaivec& bcells = _boundaryinfo.getCellsOfColor(color);
        alat::armaivec& bsidesid = _boundaryinfo.getSidesIdOfCellsOfColor(color);

        _sides[count] = s;
        // std::cerr << "boundary side found " << s << "\n";
        int i = p->second[0];
        int ii = p->second[1];
        _sides_of_cells[i][ii] = count;

        bsides[col2size_bsides[color]++] = count;
        bcells[col2size_bcells[color]++] = i;
        bsidesid[col2size_bsideids[color]++] = ii;

        count++;
      }
    }
  }
  _cells_of_sides.set_size(nsides);
  for(int is = 0; is < nsides; is++)
  {
    _cells_of_sides[is][0] = -1;
    _cells_of_sides[is][1] = -1;
  }
  for(int i = 0; i < _cells.size(); i++)
  {
    for(int ii = 0; ii < SIDESPERCELL; ii++)
    {
      int is = _sides_of_cells[i][ii];
      // std::cerr << i << " " << ii << " " << is << "@\n";
      if(_cells_of_sides[is][0] == -1)
      {
        _cells_of_sides[is][0] = i;
      }
      else
      {
        int help = _cells_of_sides[is][0];
        _cells_of_sides[is][0] = i;
        _cells_of_sides[is][1] = help;
      }
    }
  }
}

/*---------------------------------------------------------*/
template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
void FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::setVisuType(const std::string& visutype) const
{
  _visutype = visutype;
}

/*---------------------------------------------------------*/
template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
void FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::writeH5(H5::H5File& file) const
{
  assert(0);
  // // setResolution(0);
  // H5::IntType datatype( H5::PredType::NATIVE_DOUBLE );
  // datatype.setOrder( H5T_ORDER_LE );
  // int nnodes = getNNodes();
  // int ncells = getNCells();
  // int nnodespercell = NODESPERCELL;
  // int dimension = getDimension();
  //
  //
  // alat::armavec node_coord;
  // alat::armaivec connectivity;
  // hsize_t dimscoord[2];
  //
  // if(_visutype == "dg")
  // {
  //   dimscoord[0] = nnodespercell*ncells;
  //   dimscoord[1] = dimension;
  //   node_coord.set_size(dimension*ncells*nnodespercell);
  //   connectivity.set_size(nnodespercell*ncells);
  //   // coordinates
  //   // connectivities
  //   alat::Node N;
  //   int count = 0;
  //   for(int iK = 0; iK < ncells; iK++)
  //   {
  //     for(int ii = 0; ii < nnodespercell; ii++)
  //     {
  //       N = getNode( getNodeIdOfCell(iK, ii) );
  //       for(int idim = 0; idim < dimension; idim++)
  //       {
  //         node_coord[dimension*count+idim] = N[idim];
  //       }
  //       connectivity[nnodespercell*iK+ii] = count;
  //       count++;
  //     }
  //   }
  // }
  // else if(_visutype == "cg")
  // {
  //   dimscoord[0] = nnodes;
  //   dimscoord[1] = dimension;
  //   node_coord.set_size(dimension*nnodes);
  //   connectivity.set_size(nnodespercell*ncells);
  //   // coordinates
  //   alat::Node N;
  //   for(int iN = 0; iN < nnodes; iN++)
  //   {
  //     N = getNode(iN);
  //     for(int idim = 0; idim < dimension; idim++)
  //     {
  //       node_coord[dimension*iN+idim] = N[idim];
  //     }
  //   }
  //   for(int icell = 0; icell < ncells; icell++)
  //   {
  //     for(int ii = 0; ii < nnodespercell; ii++)
  //     {
  //       connectivity[nnodespercell*icell+ii] = getNodeIdOfCell(icell, ii);
  //     }
  //   }
  // }
  // else
  // {
  //   assert(0);
  // }
  //
  //
  //
  // H5::DataSpace fcoordspace(2, dimscoord );   //file dataspace
  // // write nodes coordinates
  // std::string nodes_coord_name;
  // if(dimension == 2)
  // {
  //   nodes_coord_name = "XY";
  // }
  // else
  // {
  //   nodes_coord_name = "XYZ";
  // }
  // const H5std_string nodes_coord_datasetname(nodes_coord_name);
  // H5::DataSet nodes_coord_dataset = file.createDataSet( nodes_coord_datasetname, datatype, fcoordspace );
  // nodes_coord_dataset.write( reinterpret_cast<const double*>( &( *node_coord.begin() ) ), H5::PredType::NATIVE_DOUBLE, fcoordspace, fcoordspace);
  // nodes_coord_dataset.close();
  // hsize_t dimsconnect[2];
  // dimsconnect[0] = ncells;
  // dimsconnect[1] = nnodespercell;
  // H5::DataSpace fconnectspace(2, dimsconnect );
  // // write connectivities
  // const H5std_string connectivity_datasetname("connectivities");
  // H5::DataSet connectivity_dataset = file.createDataSet( connectivity_datasetname, datatype, fconnectspace );
  // connectivity_dataset.write( reinterpret_cast<const double*>( &( *connectivity.begin() ) ), H5::PredType::NATIVE_INT, fconnectspace, fconnectspace);
  // connectivity_dataset.close();
}

/*---------------------------------------------------------*/
template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
void FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::writeMeshInfo(std::string filename) const
{
  // setResolution(0);
  std::ofstream file( filename.c_str() );
  // mesh type
  FadalightMeshEnums::meshtype meshtype = getType();
  if(meshtype == FadalightMeshEnums::TriangleMesh)
  {
    file<<"Triangle"<< "\n";
  }
  else if(meshtype == FadalightMeshEnums::QuadrilateralMesh)
  {
    file<<"Quadrilateral"<< "\n";
  }
  else if(meshtype == FadalightMeshEnums::HexahedralMesh)
  {
    file<<"Hexahedron"<< "\n";
  }
  else
  {
    _error_string( "writeMeshInfo", "unknown mesh", FadalightMeshEnums::meshTypeToString(meshtype) );
  }
  if(_visutype == "cg")
  {
    // topology dimensions
    file<<getNCells()<<'\n';
    // // meshdata dimensions
    // file<<getNNodes()<<" "<<getDimension()<<'\n';
    // nodes dimensions
    file<<getNNodes()<<'\n';
  }
  else if(_visutype == "dg")
  {
    // topology dimensions
    file<<getNCells()<<'\n';
    // // meshdata dimensions
    // file<<getNCells()*NODESPERCELL<<" "<<getDimension()<<'\n';
    // nodes dimensions
    file<<getNCells()*NODESPERCELL<<'\n';
  }
  else
  {
    assert(0);
  }
  // cells dimensions
  file<<getNCells()<<'\n';
  // nnodes per cell
  file<<NODESPERCELL<<'\n';
  file.close();
}
//
// /*---------------------------------------------------------*/
// template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
// void FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::writeEnsightGeometry(std::string filename) const
// {
//   return;
//
//   string name = filename+".case";
//   vector<string> v = alat::Tokenize(filename, "/");
//   string geoname = v[v.size()-1]+".geo";
//
//   ofstream casefile( name.c_str() );
//   casefile<<"FORMAT"<<'\n';
//   casefile<<"type: ensight gold"<<'\n';
//   casefile<<"GEOMETRY"<<'\n';
//   casefile<<"model: "<<geoname<<'\n';
//   casefile.close();
//
//   name = filename+".geo";
//   FILE* file = fopen(name.c_str(), "wb");
//
//   char buffer[80];
//   strcpy(buffer, "C Binary");
//   fwrite(buffer, sizeof( char ), 80, file);
//   strcpy(buffer, "Ensight geometry file");
//   fwrite(buffer, sizeof( char ), 80, file);
//   strcpy(buffer, "description");
//   fwrite(buffer, sizeof( char ), 80, file);
//   strcpy(buffer, "node id off");
//   fwrite(buffer, sizeof( char ), 80, file);
//   strcpy(buffer, "element id off");
//   fwrite(buffer, sizeof( char ), 80, file);
//   strcpy(buffer, "part");
//   fwrite(buffer, sizeof( char ), 80, file);
//   int part_number = 1;
//   fwrite(&part_number, sizeof( int ), 1, file);
//   strcpy(buffer, "description");
//   fwrite(buffer, sizeof( char ), 80, file);
//   strcpy(buffer, "coordinates");
//   fwrite(buffer, sizeof( char ), 80, file);
//   int nn = getNNodes();
//   fwrite(&nn, sizeof( int ), 1, file);
//   typename alat::Vector<alat::Node>::const_iterator firstn = getNodes().begin();
//   typename alat::Vector<alat::Node>::const_iterator lastn = getNodes().end();
//   while(firstn != lastn)
//   {
//     float vx = float( ( *firstn ).x() );
//     fwrite(&vx, sizeof( float ), 1, file);
//     firstn++;
//   }
//   firstn = getNodes().begin();
//   while(firstn != lastn)
//   {
//     float vy = float( ( *firstn ).y() );
//     fwrite(&vy, sizeof( float ), 1, file);
//     firstn++;
//   }
//   firstn = getNodes().begin();
//   while(firstn != lastn)
//   {
//     float vz = float( ( *firstn ).z() );
//     fwrite(&vz, sizeof( float ), 1, file);
//     firstn++;
//   }
//   strcpy( buffer, getEnsightType().c_str() );
//   fwrite(buffer, sizeof( char ), 80, file);
//   int nel = getNCells();
//   fwrite(&nel, sizeof( int ), 1, file);
//   for(int i = 0; i < nel; i++)
//   {
//     int npc = getNNodesPerCell(i);
//     for(int ii = 0; ii < npc; ii++)
//     {
//       int in = getNodeIdOfCell(i, ii)+1;
//       fwrite(&in, sizeof( int ), 1, file);
//     }
//   }
// }

/*--------------------------------------------------------------------------*/
template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
void FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::writeVtk(std::string filename) const
{
  string name = filename;
  name += ".vtk";

  ofstream file( name.c_str() );
  assert( file.is_open() );


  file<<"# vtk DataFile Version 4.0 "<<std::endl;
  file<<"output from QuadrilateralMesh"<<std::endl;
  file<<"ASCII"<<std::endl;
  //     file << "binary" << std::endl;
  file<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
  file<<std::endl;

  int nn = getNNodes();

  file<<"POINTS "<<nn;
  file<<" FLOAT"<<std::endl;
  for(int i = 0; i < nn; i++)
  {
    const alat::Node& v = getNode(i);
    file<<v.x()<<" "<<v.y()<<" "<<v.z()<<" "<<std::endl;
  }
  file<<std::endl;

  int ne = getNCells();
  int nle = NODESPERCELL;
  int lenght = ne*( nle+1 );

  file<<std::endl<<"CELLS "<<ne<<" "<<lenght<<std::endl;

  for(int ie = 0; ie < ne; ie++)
  {
    file<<nle<<" ";
    for(int ii = 0; ii < nle; ii++)
    {
      file<<getNodeIdOfCell(ie, ii)<<" ";
    }
    file<<std::endl;
  }
  file<<std::endl<<"CELL_TYPES "<<ne<<std::endl;
  for(int c = 0; c < ne; c++)
  {
    file<<getVtkType()<<" ";
  }
  file<<std::endl;

  file.close();
}

/*---------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
void FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::writeBoundaryVtk(std::string filename) const
{
  string name = filename;
  name += "-boundary.vtk";

  ofstream file( name.c_str() );
  assert( file.is_open() );

  file<<"# vtk DataFile Version 4.0 "<<std::endl;
  file<<"output from QuadrilateralMesh"<<std::endl;
  file<<"ASCII"<<std::endl;
  //     file << "binary" << std::endl;
  file<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
  file<<std::endl;

  int nn = getNNodes();

  file<<"POINTS "<<nn;
  file<<" FLOAT"<<std::endl;
  for(int i = 0; i < nn; i++)
  {
    const alat::Node& v = getNode(i);
    file<<v.x()<<" "<<v.y()<<" "<<v.z()<<" "<<std::endl;
  }
  file<<std::endl;

  const FadalightMesh::BoundaryInfo* BI = getBoundaryInfo();
  int nsides = BI->getNSides();
  const alat::armaivec& colors = BI->getColors();

  int nle = NODESPERSIDE;
  int lenght = nsides*( nle+1 );
  file<<std::endl<<"CELLS "<<nsides<<" "<<lenght<<std::endl;

  for(int i = 0; i < colors.size(); i++)
  {
    int color = colors[i];
    const alat::armaivec& sides = BI->getSidesOfColor(color);
    for(int j = 0; j < sides.size(); j++)
    {
      file<<nle<<" ";
      for(int ii = 0; ii < nle; ii++)
      {
        file<<getNodeIdOfSide(sides[j], ii)<<" ";
      }
      file<<std::endl;
    }
  }
  file<<std::endl<<"CELL_TYPES "<<nsides<<std::endl;
  for(int c = 0; c < nsides; c++)
  {
    file<<getBoundaryVtkType()<<" ";
  }
  file<<std::endl;
  file<<std::endl<<"CELL_DATA "<<nsides<<std::endl;
  file<<std::endl<<"SCALARS "<<" bdry_colors "<<" int "<<1<<std::endl;
  file<<std::endl<<"LOOKUP_TABLE default"<<std::endl;
  for(int i = 0; i < colors.size(); i++)
  {
    int color = colors[i];
    const alat::armaivec& sides = BI->getSidesOfColor(color);
    for(int j = 0; j < sides.size(); j++)
    {
      file<<color<<" ";
      file<<std::endl;
    }
    file<<std::endl;
  }

  file.close();
}

/*--------------------------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
int FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::localIndexOfNodeInCell(int iK, int in) const
{
  int ii = -1;
  for(int jj = 0; jj < getNNodesPerCell(iK); jj++)
  {
    if(getNodeIdOfCell(iK, jj) == in)
    {
      ii = jj;
      break;
    }
  }
  assert(ii != -1);
  return ii;
}

/*--------------------------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
int FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::localIndexOfSideInCell(int iK, int is) const
{
  int ii = -1;
  for(int jj = 0; jj < getNSidesPerCell(iK); jj++)
  {
    if(getSideIdOfCell(iK, jj) == is)
    {
      ii = jj;
      break;
    }
  }
  assert(ii != -1);
  return ii;
}

/*--------------------------------------------------------------------------*/

template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
int FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::localIndexOfEdgeInCell(int iK, int ie) const
{
  int ii = -1;
  for(int jj = 0; jj < getNEdgesPerCell(iK); jj++)
  {
    if(getEdgeIdOfCell(iK, jj) == ie)
    {
      ii = jj;
      break;
    }
  }
  assert(ii != -1);
  return ii;
}

/*--------------------------------------------------------------------------*/
// line mesh
#define DIM          1
#define NODESPERCELL 2
#define SIDESPERCELL 2
#define NODESPERSIDE 1
template class FadalightMesh::FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>;
#undef DIM
#undef NODESPERCELL
#undef SIDESPERCELL
#undef NODESPERSIDE

// triangle mesh
#define DIM          2
#define NODESPERCELL 3
#define SIDESPERCELL 3
#define NODESPERSIDE 2
template class FadalightMesh::FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>;
#undef DIM
#undef NODESPERCELL
#undef SIDESPERCELL
#undef NODESPERSIDE

// quadrilateral mesh
#define DIM          2
#define NODESPERCELL 4
#define SIDESPERCELL 4
#define NODESPERSIDE 2
template class FadalightMesh::FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>;
#undef DIM
#undef NODESPERCELL
#undef SIDESPERCELL
#undef NODESPERSIDE

// tetrahedral mesh
#define DIM          3
#define NODESPERCELL 4
#define SIDESPERCELL 4
#define NODESPERSIDE 3
template class FadalightMesh::FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>;
#undef DIM
#undef NODESPERCELL
#undef SIDESPERCELL
#undef NODESPERSIDE

// hexahedral mesh
#define DIM          3
#define NODESPERCELL 8
#define SIDESPERCELL 6
#define NODESPERSIDE 4
template class FadalightMesh::FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>;
#undef DIM
#undef NODESPERCELL
#undef SIDESPERCELL
#undef NODESPERSIDE
