#include  "FadalightMesh/hexahedralmesh.hpp"
#include  "FadalightMesh/multilevelmesh.hpp"
#include  "FadalightMesh/quadrilateralmesh.hpp"
#include  "FadalightMesh/quadtotri.hpp"
#include  <cassert>
#include  <fstream>
#include  <sstream>

using namespace FadalightMesh;

/*--------------------------------------------------------------------------*/

MultiLevelMesh::~MultiLevelMesh()
{}

MultiLevelMesh::MultiLevelMesh(std::string type) : FadalightMesh::MeshInterface(), _type(type), _activelevel(0)
{}

MultiLevelMesh::MultiLevelMesh( const MultiLevelMesh& multilevelmesh) : FadalightMesh::MeshInterface(multilevelmesh)
{
  assert(0);
}

MultiLevelMesh& MultiLevelMesh::operator=( const MultiLevelMesh& multilevelmesh)
{
  FadalightMesh::MeshInterface::operator=(multilevelmesh);
  assert(0);
  return *this;
}

std::string MultiLevelMesh::getClassName() const
{
  return "MultiLevelMesh";
}

MultiLevelMesh* MultiLevelMesh::clone() const
{
  return new MultiLevelMesh(*this);
}

FadalightMeshEnums::meshtype MultiLevelMesh::getType() const
{
  return _meshes[_activelevel]->getType();
}

alat::Vector<FadalightMesh::MeshInterface*>& MultiLevelMesh::getMeshes()
{
  return _meshes;
}

/*--------------------------------------------------------------------------*/

// bool MultiLevelMesh::cellIsCurved(int iK) const
// {
//   return _meshes[_activelevel]->cellIsCurved(iK);
// }

int MultiLevelMesh::getNLevels() const
{
  return _meshes.size();
}

const FadalightMesh::MeshInterface* MultiLevelMesh::getMesh(int level) const
{
  return _meshes[level];
}

FadalightMesh::MeshInterface* MultiLevelMesh::getMesh(int level)
{
  return _meshes[level];
}

FadalightMesh::MeshInterface*& MultiLevelMesh::getMeshPointer(int level)
{
  return _meshes[level];
}

/*--------------------------------------------------------------------------*/

std::string MultiLevelMesh::getInfo() const
{
  std::string info;
  std::stringstream ss1, ss2, ss3;
  ss1<<getNNodes();
  ss2<<getNCells();
  ss3<<getNSides();
  return ss1.str()+"_"+ss2.str()+"_"+ss3.str();
}

/*--------------------------------------------------------------------------*/
void MultiLevelMesh::setVisuType(const std::string& visutype) const
{
  for(int level = 0; level < getNLevels(); level++)
  {
    getMesh(level)->setVisuType(visutype);
  }
}

bool MultiLevelMesh::isMultilevel() const
{
  return true;
}

bool MultiLevelMesh::geometryObjectExists(std::string name) const
{
  return getMesh(0)->geometryObjectExists(name);
}

const FadalightMesh::GeometryObject* MultiLevelMesh::getGeometryObject(std::string name) const
{
  return getMesh(0)->getGeometryObject(name);
}

int MultiLevelMesh::getDimension() const
{
  return getMesh(0)->getDimension();
}

const FadalightMesh::BoundaryInfo* MultiLevelMesh::getBoundaryInfo() const
{
  return getMesh(_activelevel)->getBoundaryInfo();
}

void MultiLevelMesh::setResolution(int level) const
{
  _activelevel = level;
}
int MultiLevelMesh::getResolution() const
{
  return _activelevel;
}

alat::Node MultiLevelMesh::getNodeOfCell(int iK) const
{
  return getMesh(_activelevel)->getNodeOfCell(iK);
}

alat::Node MultiLevelMesh::getNodeOfSide(int iS) const
{
  return getMesh(_activelevel)->getNodeOfSide(iS);
}

int MultiLevelMesh::getNNodes() const
{
  return getMesh(_activelevel)->getNNodes();
}

int MultiLevelMesh::getNEdges() const
{
  return getMesh(_activelevel)->getNEdges();
}

int MultiLevelMesh::getNSides() const
{
  return getMesh(_activelevel)->getNSides();
}

int MultiLevelMesh::getNCells() const
{
  return getMesh(_activelevel)->getNCells();
}

int MultiLevelMesh::getNNodesPerCell(int iK) const
{
  return getMesh(_activelevel)->getNNodesPerCell(iK);
}

int MultiLevelMesh::getNNodesPerSide(int iS) const
{
  return getMesh(_activelevel)->getNNodesPerSide(iS);
}

int MultiLevelMesh::getNSidesPerCell(int iK) const
{
  return getMesh(_activelevel)->getNSidesPerCell(iK);
}

int MultiLevelMesh::getNEdgesPerCell(int iK) const
{
  return getMesh(_activelevel)->getNEdgesPerCell(iK);
}

int MultiLevelMesh::getNEdgesPerSide(int iS) const
{
  return getMesh(_activelevel)->getNEdgesPerSide(iS);
}

std::string MultiLevelMesh::getCellType() const
{
  return getMesh(_activelevel)->getCellType();
}

int MultiLevelMesh::getNodeIdOfCell(int i, int ii) const
{
  return getMesh(_activelevel)->getNodeIdOfCell(i, ii);
}

int MultiLevelMesh::getNodeIdOfSide(int i, int ii) const
{
  return getMesh(_activelevel)->getNodeIdOfSide(i, ii);
}

int MultiLevelMesh::getNodeIdOfEdge(int i, int ii) const
{
  return getMesh(_activelevel)->getNodeIdOfEdge(i, ii);
}

int MultiLevelMesh::getSideIdOfCell(int i, int ii) const
{
  return getMesh(_activelevel)->getSideIdOfCell(i, ii);
}

int MultiLevelMesh::getEdgeIdOfCell(int i, int ii) const
{
  return getMesh(_activelevel)->getEdgeIdOfCell(i, ii);
}

int MultiLevelMesh::getEdgeIdOfSide(int i, int ii) const
{
  return getMesh(_activelevel)->getEdgeIdOfSide(i, ii);
}

int MultiLevelMesh::getCellIdOfSide(int i, int ii) const
{
  // std::cerr << "_activelevel " << _activelevel << "\n";
  return getMesh(_activelevel)->getCellIdOfSide(i, ii);
}

int MultiLevelMesh::getLocalIndexOfSideInCell(int iK, int iS) const
{
  return getMesh(_activelevel)->getLocalIndexOfSideInCell(iK, iS);
}

void MultiLevelMesh::getNodesOfCell(int iK, alat::Vector<alat::Node>& F) const
{
  getMesh(_activelevel)->getNodesOfCell(iK, F);
}

int MultiLevelMesh::getCouplingOffset(int iS) const
{
  return getMesh(_activelevel)->getCouplingOffset(iS);
}

void MultiLevelMesh::writeH5(H5::H5File& file) const
{
  getMesh(0)->writeH5(file);
}

void MultiLevelMesh::writeMeshInfo(std::string filename) const
{
  assert(0);
  // std::cerr << "MultiLevelMesh::writeMeshInfo() " << getMesh(0)->getNCells() << "\n";
  getMesh(0)->writeMeshInfo(filename);
}

/*--------------------------------------------------------------------------*/

void MultiLevelMesh::readFadalightMesh(const std::string& basefilename)
{
  std::string multimeshdir = basefilename+".fadalightmesh/MultiMesh";
  // std::cerr << "MultiLevelMesh::readFadalightMesh() multimeshdir " << multimeshdir << "\n";
  std::string filename = multimeshdir+"/n";
  std::ifstream file( filename.c_str() );
  if( not file.is_open() )
  {
    _error_string("readFadalightMesh", "cannot open file", filename);
  }
  int nlevels;
  file>>nlevels;
  file.close();
  if(nlevels == 0)
  {
    _error_string("readFadalightMesh", "MultiMesh but 0 levels");
  }
  // std::cerr << "MultiLevelMesh::readFadalightMesh() nlevels " << nlevels << "\n";
  _meshes.set_size(nlevels);
  for(int level = 0; level < nlevels; level++)
  {
    if(_type == "FadalightMesh::QuadrilateralMesh")
    {
      _meshes[level] = new FadalightMesh::QuadrilateralMesh;
    }
    else if(_type == "FadalightMesh::QuadToTri")
    {
      _meshes[level] = new FadalightMesh::QuadToTri;
    }
    else if(_type == "FadalightMesh::HexahedralMesh")
    {
      _meshes[level] = new FadalightMesh::HexahedralMesh;
    }
    else
    {
      _error_string("readFadalightMesh", "unknown meshtype: "+_type);
    }
  }
  _meshes[0]->readFadalightMesh(basefilename);
  for(int level = 1; level < nlevels; level++)
  {
    std::stringstream ss0;
    ss0<<level;
    std::string levelmeshname = multimeshdir+"/mesh."+ss0.str();
    _meshes[level]->readFadalightMesh(levelmeshname);
  }
}

/*--------------------------------------------------------------------------*/
void MultiLevelMesh::writeFadalightMesh(const std::string& basefilename, arma::file_type datatype) const
{
  std::cerr << "MultiLevelMesh::writeFadalightMesh() basefilename " << basefilename << "\n";
  std::string multimeshdir = basefilename+".fadalightmesh/MultiMesh";
  std::string command;
  command = "mkdir -p "+multimeshdir;
  system( command.c_str() );
  std::cerr << "MultiLevelMesh::writeFadalightMesh() multimeshdir " << multimeshdir << "\n";
  std::string filename = multimeshdir+"/n";
  std::ofstream file( filename.c_str() );
  assert( file.is_open() );
  int nlevels;
  file<<nlevels;
  file.close();
  _meshes[0]->writeFadalightMesh(basefilename, datatype);
  for(int level = 1; level < nlevels; level++)
  {
    std::stringstream ss0;
    ss0<<level;
    std::string levelmeshname = multimeshdir+"/mesh."+ss0.str();
    _meshes[level]->writeFadalightMesh(levelmeshname, datatype);
  }
}

/*--------------------------------------------------------------------------*/

int MultiLevelMesh::localIndexOfNodeInCell(int iK, int in) const
{
  return getMesh(_activelevel)->localIndexOfNodeInCell(iK, in);
}

int MultiLevelMesh::localIndexOfSideInCell(int iK, int is) const
{
  return getMesh(_activelevel)->localIndexOfSideInCell(iK, is);
}

int MultiLevelMesh::localIndexOfEdgeInCell(int iK, int ie) const
{
  return getMesh(_activelevel)->localIndexOfEdgeInCell(iK, ie);
}

const alat::Node& MultiLevelMesh::getNodeOfSide(int is, int ii) const
{
  return getMesh(_activelevel)->getNodeOfSide(is, ii);
}

const alat::Node& MultiLevelMesh::getNode(int i) const
{
  return getMesh(_activelevel)->getNode(i);
}
void MultiLevelMesh::getLocalIndicesOfSidesInCell(alat::armaivec& sideindex_a, alat::armaivec& sideindex_e) const
{
  // std::cerr << "@@@@@@@@@@@@@@ _activelevel="<<_activelevel<<"\n";
  getMesh(_activelevel)->getLocalIndicesOfSidesInCell(sideindex_a, sideindex_e);
}
void MultiLevelMesh::getLocalIndicesOfSidesAndDiagonalsInCell(alat::armaivec& sideindex_a, alat::armaivec& sideindex_e) const
{
  // std::cerr << "@@@@@@@@@@@@@@ _activelevel="<<_activelevel<<"\n";
  getMesh(_activelevel)->getLocalIndicesOfSidesAndDiagonalsInCell(sideindex_a, sideindex_e);
}
