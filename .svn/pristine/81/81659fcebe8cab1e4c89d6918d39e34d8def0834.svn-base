#include  "FadalightMesh/meshinterface.hpp"
#include  <cassert>

using namespace FadalightMesh;

alat::Node FadalightMesh::MeshInterface::_dummynode;

/*--------------------------------------------------------------------------*/
MeshInterface::~MeshInterface()
{}
MeshInterface::MeshInterface() : alat::InterfaceBase()
{}
MeshInterface::MeshInterface( const MeshInterface& meshinterface) : alat::InterfaceBase(meshinterface)
{
  assert(0);
}
MeshInterface& MeshInterface::operator=( const MeshInterface& meshinterface)
{
  InterfaceBase::operator=(meshinterface);
  assert(0);
  return *this;
}
std::string MeshInterface::getInterfaceName() const
{
  return "MeshInterface";
}


/*--------------------------------------------------------------------------*/
const std::string& MeshInterface::getInfileName() const
{
  _notWritten("getInfileName");
}

/*--------------------------------------------------------------------------*/
std::string MeshInterface::getCellType() const
{
  _notWritten("getCellType");
}

/*--------------------------------------------------------------------------*/

void MeshInterface::getNodesOfCell(int iK, alat::Vector<alat::Node>& F) const
{
  _notWritten("getNodesOfCell");
}

/*--------------------------------------------------------------------------*/

void MeshInterface::setResolution(int level) const
{
  _notWritten("setResolution");
}
int MeshInterface::getResolution() const
{
  _notWritten("getResolution");
}

/*--------------------------------------------------------------------------*/

int MeshInterface::getNLevels() const
{
  _notWritten("getNLevels");
}

/*--------------------------------------------------------------------------*/

alat::Node MeshInterface::getNodeOfSide(int iS) const
{
  _notWritten("getNodeOfSide");
  return _dummynode;
}

/*--------------------------------------------------------------------------*/

int MeshInterface::getNodeIdOfSideOfCell(int iK, int iis, int ii) const
{
  _notWritten("getNodeIdOfSideOfCell");
  return 0;
}

/*--------------------------------------------------------------------------*/

alat::Node MeshInterface::getNodeOfCell(int iK) const
{
  _notWritten("getNodeOfCell");
  return _dummynode;
}

/*--------------------------------------------------------------------------*/

alat::Node MeshInterface::getNodeOfEdge(int iK) const
{
  _notWritten("getNodeOfEdge");
  return _dummynode;
}

/*--------------------------------------------------------------------------*/

int MeshInterface::getNNodes() const
{
  _notWritten("getNNodes");
  return 0;
}

/*--------------------------------------------------------------------------*/

int MeshInterface::getNEdges() const
{
  _notWritten("getNEdges");
  return 0;
}

/*--------------------------------------------------------------------------*/

int MeshInterface::getNSides() const
{
  _notWritten("getNSides");
  return 0;
}

/*--------------------------------------------------------------------------*/

int MeshInterface::getNCells() const
{
  _notWritten("getNCells");
  return 0;
}

/*--------------------------------------------------------------------------*/

const alat::Node& MeshInterface::getNode(int i) const
{
  _notWritten("getNode");
  return _dummynode;
}

/*--------------------------------------------------------------------------*/

const alat::Node& MeshInterface::getNodeOfCell(int iK, int ii) const
{
  return getNode( getNodeIdOfCell(iK, ii) );
}

/*--------------------------------------------------------------------------*/

const alat::Node& MeshInterface::getNodeOfEdge(int iE, int ii) const
{
  return getNode( getNodeIdOfEdge(iE, ii) );
}

/*--------------------------------------------------------------------------*/

const alat::Node& MeshInterface::getNodeOfSide(int iS, int ii) const
{
  return getNode( getNodeIdOfSide(iS, ii) );
}

/*--------------------------------------------------------------------------*/

int MeshInterface::getNNodesPerCell(int i) const
{
  _notWritten("getNNodesPerCell");
  return 0;
}

/*--------------------------------------------------------------------------*/

int MeshInterface::getNNodesPerSide(int i) const
{
  _notWritten("getNNodesPerSide");
  return 0;
}

/*--------------------------------------------------------------------------*/

int MeshInterface::getNSidesPerCell(int i) const
{
  _notWritten("getNSidesPerCell");
  return 0;
}

/*--------------------------------------------------------------------------*/

int MeshInterface::getNEdgesPerCell(int i) const
{
  _notWritten("getNEdgesPerCell");
  return 0;
}

/*--------------------------------------------------------------------------*/

int MeshInterface::getNEdgesPerSide(int i) const
{
  _notWritten("getNEdgesPerSide");
  return 0;
}

/*--------------------------------------------------------------------------*/

const FadalightMesh::BoundaryInfo* MeshInterface::getBoundaryInfo() const
{
  _notWritten("getBoundaryInfo");
  return NULL;
}
const FadalightMesh::CurvedInteriorSideInfo* MeshInterface::getCurvedInteriorSideInfo() const
{
  _notWritten("getCurvedInteriorSideInfo");
}

/*--------------------------------------------------------------------------*/

int MeshInterface::getNodeIdOfCell(int i, int ii) const
{
  _notWritten("getNodeIdOfCell");
  return 0;
}

/*--------------------------------------------------------------------------*/

int MeshInterface::getNodeIdOfSide(int i, int ii) const
{
  _notWritten("getNodeIdOfSide");
  return 0;
}

/*--------------------------------------------------------------------------*/

int MeshInterface::getNodeIdOfEdge(int i, int ii) const
{
  _notWritten("getNodeIdOfEdge");
  return 0;
}

/*--------------------------------------------------------------------------*/

int MeshInterface::getSideIdOfCell(int i, int ii) const
{
  _notWritten("getSideIdOfCell");
  return 0;
}

/*--------------------------------------------------------------------------*/

int MeshInterface::getEdgeIdOfCell(int i, int ii) const
{
  _notWritten("getEdgeIdOfCell");
  return 0;
}

/*--------------------------------------------------------------------------*/

int MeshInterface::getEdgeIdOfSide(int i, int ii) const
{
  _notWritten("getEdgeIdOfSide");
  return 0;
}

/*--------------------------------------------------------------------------*/

int MeshInterface::getCellIdOfSide(int i, int ii) const
{
  _notWritten("getCellIdOfSide");
  return 0;
}

/*--------------------------------------------------------------------------*/

int MeshInterface::getLocalIndexOfSideInCell(int iK, int iS) const
{
  _notWritten("getLocalIndexOfSideInCell");
  return 0;
}

/*--------------------------------------------------------------------------*/

bool MeshInterface::geometryObjectExists(std::string name) const
{
  _notWritten("geometryObjectExists");
  return 0;
}

/*--------------------------------------------------------------------------*/

const FadalightMesh::GeometryObject* MeshInterface::getGeometryObject(std::string name) const
{
  _notWritten("getGeometryObject");
  return NULL;
}

/*--------------------------------------------------------------------------*/

FadalightMesh::GeometryObject* MeshInterface::getGeometryObject(std::string name)
{
  _notWritten("getGeometryObject");
  return NULL;
}

// /*--------------------------------------------------------------------------*/
//
// const FadalightMesh::GeometryObjectsConstructorInterface* MeshInterface::getGeometryObjectsConstructor() const
// {
//   _notWritten("getGeometryObjectsConstructor");
//   return NULL;
// }
//
// /*--------------------------------------------------------------------------*/
//
// FadalightMesh::GeometryObjectsConstructorInterface*& MeshInterface::getGeometryObjectsConstructor()
// {
//   _notWritten("getGeometryObjectsConstructor");
//   return dummy_geometryobjects_constructor;
// }

/*--------------------------------------------------------------------------*/

void MeshInterface::createGeometryObject(std::string name)
{
  _notWritten("CreateGeometryObject");
}

/*--------------------------------------------------------------------------*/

const FadalightMesh::CurvedBoundaryInformation* MeshInterface::getCurvedBoundaryInformation() const
{
  _notWritten("getCurvedBoundaryInformation");
  return NULL;
}

/*--------------------------------------------------------------------------*/

void MeshInterface::read(std::string filename)
{
  _notWritten("read");
}

/*--------------------------------------------------------------------------*/

void MeshInterface::write(std::string filename) const
{
  _notWritten("write");
}

/*--------------------------------------------------------------------------*/

void MeshInterface::readFadalightMesh(const std::string& basefilename)
{
  _notWritten("readFadalightMesh");
}

/*--------------------------------------------------------------------------*/
void MeshInterface::writeFadalightMesh(const std::string& basefilename, arma::file_type datatype) const
{
  _notWritten("writeFadalightMesh");
}

/*--------------------------------------------------------------------------*/

void MeshInterface::writeSimpleMesh(std::string basefilename, arma::file_type datatype) const
{
  _notWritten("writeSimpleMesh");
}

/*--------------------------------------------------------------------------*/

void MeshInterface::writeVtk(std::string filename) const
{
  _notWritten("writeVtk");
}
//
// /*--------------------------------------------------------------------------*/
//
// void MeshInterface::writeEnsightGeometry(std::string filename) const
// {
//   _notWritten("writeEnsightGeometry");
// }

/*--------------------------------------------------------------------------*/
void MeshInterface::writeMeshInfo(std::string filename) const
{
  _notWritten("writeMeshInfo");
}

/*--------------------------------------------------------------------------*/

void MeshInterface::writeBoundaryVtk(std::string filename) const
{
  _notWritten("writeBoundaryVtk");
}

/*--------------------------------------------------------------------------*/

void MeshInterface::writeH5(H5::H5File& file) const
{
  _notWritten("writeH5");
}

/*--------------------------------------------------------------------------*/

int MeshInterface::getVtkType() const
{
  _notWritten("getVtkType");
  return 0;
}

// /*--------------------------------------------------------------------------*/
//
// std::string MeshInterface::getEnsightType() const
// {
//   _notWritten("getEnsightType");
//   return "none";
// }

/*--------------------------------------------------------------------------*/

int MeshInterface::getBoundaryVtkType() const
{
  _notWritten("getBoundaryVtkType");
  return 0;
}

// /*--------------------------------------------------------------------------*/
//
// void MeshInterface::computeCellConnectivity(alat::SparsityPattern& SPC) const
// {
//   _notWritten("computeCellConnectivity");
// }
//
// /*--------------------------------------------------------------------------*/
//
// void MeshInterface::computeCellNeighbours()
// {
//   _notWritten("computeCellNeighbours(");
// }

/*--------------------------------------------------------------------------*/

bool MeshInterface::cellIsCurved(int iK) const
{
  _notWritten("cellIsCurved");
  return 0;
}
//
// /*--------------------------------------------------------------------------*/
//
// void MeshInterface::addGeometryObjects(const std::string& filename, const alat::StringVector& names_of_geometry_objects, std::string datatype)
// {
//   _notWritten("addGeometryObjects");
// }
//
// /*--------------------------------------------------------------------------*/
//
// void MeshInterface::addGeometryObject(const std::string& name, FadalightMesh::GeometryObject* geo)
// {
//   _notWritten("addGeometryObject");
// }

// /*--------------------------------------------------------------------------*/
//
// void MeshInterface::writeEnsightGeometryObjectDescFile(const std::string& basefilename)
// {
//   _notWritten("writeEnsightGeometryObjectDescFile");
// }

/*--------------------------------------------------------------------------*/

int MeshInterface::findNeightborHangingCells(int iK, int iS, alat::Node pt)
{
  _notWritten("findNeightborHangingCells");
  return 0;
}

/*--------------------------------------------------------------------------*/

int MeshInterface::getLocalNodeIndiceOfSide(int ii, int isl) const
{
  _notWritten("getLocalNodeIndiceOfSide");
  return 0;
}

//
// /*--------------------------------------------------------------------------*/
//
// void MeshInterface::getMeshSizeForStabilization(double& hs, int iS, int iK, int iil) const
// {
//   _notWritten("getMeshSizeForStabilization");
// }

/*--------------------------------------------------------------------------*/

int MeshInterface::getCouplingOffset(int iS) const
{
  _notWritten("getCouplingOffset");
}

/*--------------------------------------------------------------------------*/

alat::Vector<alat::Node>& MeshInterface::getAllNodes()
{
  _notWritten("getAllNodes");
}

/*--------------------------------------------------------------------------*/

int MeshInterface::localIndexOfNodeInCell(int iK, int in) const
{
  _notWritten("localIndexOfNodeInCell");
}

int MeshInterface::localIndexOfSideInCell(int iK, int is) const
{
  _notWritten("localIndexOfSideInCell");
}

int MeshInterface::localIndexOfEdgeInCell(int iK, int ie) const
{
  _notWritten("localIndexOfEdgeInCell");
}

/*--------------------------------------------------------------------------*/

// const alat::SparsityPatternFixArray<2>& MeshInterface::getCellNeighbours() const
// {
//   _notWritten("getCellNeighbours");
// }

FadalightMeshEnums::meshtype MeshInterface::getType() const
{
  _notWritten("getType");
}

void MeshInterface::getLocalIndicesOfSidesInCell(alat::armaivec& sideindex_a, alat::armaivec& sideindex_e) const
{
  _notWritten("getLocalIndicesOfSidesInCell");
}

void MeshInterface::getLocalIndicesOfSidesAndDiagonalsInCell(alat::armaivec& sideindex_a, alat::armaivec& sideindex_e) const
{
  _notWritten("getLocalIndicesOfSidesAndDiagonalsInCell");
}
