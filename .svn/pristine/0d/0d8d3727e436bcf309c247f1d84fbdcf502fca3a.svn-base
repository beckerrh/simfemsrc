// #include  "Alat/iomanager.hpp"
#include  "FadalightMesh/getmeshtype.hpp"
#include  "FadalightMesh/multilevelmesh.hpp"
#include  "FadalightMesh/singlemultilevelmeshcomposition.hpp"
// #include  "H5Cpp.hpp"
#include  <cassert>
#include  <fstream>
#include  <iomanip>
#include  <sstream>

using namespace FadalightMesh;

/*--------------------------------------------------------------------------*/

SingleMultiLevelMeshComposition::~SingleMultiLevelMeshComposition()
{
  if(_mesh)
  {
    delete _mesh;
    _mesh = NULL;
  }
}

SingleMultiLevelMeshComposition::SingleMultiLevelMeshComposition() : FadalightMesh::MeshCompositionInterface()
{}
SingleMultiLevelMeshComposition::SingleMultiLevelMeshComposition( const SingleMultiLevelMeshComposition& singlemeshcomposition) : FadalightMesh::MeshCompositionInterface(singlemeshcomposition)
{
  assert(0);
}

SingleMultiLevelMeshComposition& SingleMultiLevelMeshComposition::operator=( const SingleMultiLevelMeshComposition& singlemeshcomposition)
{
  FadalightMesh::MeshCompositionInterface::operator=(singlemeshcomposition);
  assert(0);
  return *this;
}

std::string SingleMultiLevelMeshComposition::getClassName() const
{
  return "SingleMultiLevelMeshComposition";
}

SingleMultiLevelMeshComposition* SingleMultiLevelMeshComposition::clone() const
{
  return new SingleMultiLevelMeshComposition(*this);
}

/*--------------------------------------------------------------------------*/

std::string SingleMultiLevelMeshComposition::getInfo() const
{
  std::stringstream ss1, ss2, ss3;
  ss1 << getNDomains();
  ss2 << getNCouplingMeshes();
  ss3 << getNCells();
  return ss1.str()+"_"+ss2.str()+"_"+ss3.str();
}

/*--------------------------------------------------------------------------*/
int SingleMultiLevelMeshComposition::getDimension() const
{
  const FadalightMesh::MultiLevelMesh* mlmesh = dynamic_cast<const FadalightMesh::MultiLevelMesh*>( _mesh );
  assert(mlmesh);
  return mlmesh->getMesh(0)->getDimension();
}

/*--------------------------------------------------------------------------*/
int SingleMultiLevelMeshComposition::getNCells() const
{
  const FadalightMesh::MultiLevelMesh* mlmesh = dynamic_cast<const FadalightMesh::MultiLevelMesh*>( _mesh );
  assert(mlmesh);
  return mlmesh->getMesh(0)->getNCells();
}

/*--------------------------------------------------------------------------*/

int SingleMultiLevelMeshComposition::getNDomains() const
{
  return 1;
}

/*--------------------------------------------------------------------------*/

int SingleMultiLevelMeshComposition::getNCouplingMeshes() const
{
  return 0;
}

/*--------------------------------------------------------------------------*/

const FadalightMesh::MeshInterface* SingleMultiLevelMeshComposition::getMesh(int i) const
{
  return _mesh;
}

/*--------------------------------------------------------------------------*/

const FadalightMesh::CouplingMeshInterface* SingleMultiLevelMeshComposition::getCouplingMesh(int i) const
{
  return NULL;
}

/*--------------------------------------------------------------------------*/

void SingleMultiLevelMeshComposition::constructFadalightMesh(const std::string& meshname)
{
  alat::StringPair p = getMeshType(meshname);
  std::string type = p.first;
  std::string datatype = p.second;
  assert(datatype == "ascii" || datatype == "binary");

  _mesh = new FadalightMesh::MultiLevelMesh(type);
}

/*--------------------------------------------------------------------------*/

void SingleMultiLevelMeshComposition::read(const std::string& basefilename)
{
  _mesh->readFadalightMesh(basefilename);
}

/*--------------------------------------------------------------------------*/

void SingleMultiLevelMeshComposition::write(const std::string& basefilename, arma::file_type datatype) const
{
  _mesh->writeFadalightMesh(basefilename, datatype);
}

/*--------------------------------------------------------------------------*/

void SingleMultiLevelMeshComposition::writeMeshInfo(std::string filename, std::string blockfilename) const
{
  // std::cerr << "@@@ SingleMultiLevelMeshComposition::writeMeshInfo()\n";
  const FadalightMesh::MultiLevelMesh* mlmesh = dynamic_cast<const FadalightMesh::MultiLevelMesh*>( _mesh );
  assert(mlmesh);
  int nlevels = _mesh->getNLevels();
  std::ofstream file( filename.c_str() );
  file <<_mesh->getDimension()<< "\n";
  file <<nlevels<< "\n";
  file << getNDomains()<< "\n";
  file.close();
  int iblock = 0;
  for(int level=0;level<nlevels;level++)
  {
    // std::string blockfilename = alat::IoManager::getFileNameOnBlock( filename, iblock, level);
    mlmesh->getMesh(level)->writeMeshInfo(blockfilename);
  }
}

/*--------------------------------------------------------------------------*/

void SingleMultiLevelMeshComposition::writeH5(std::string filename) const
{
  assert(0);
  // const FadalightMesh::MultiLevelMesh* mlmesh = dynamic_cast<const FadalightMesh::MultiLevelMesh*>( _mesh );
  // assert(mlmesh);
  // int ndomains = getNDomains();
  // for(int idomain = 0; idomain < ndomains; idomain++)
  // {
  //   for(int level=0;level<mlmesh->getNLevels();level++)
  //   {
  //     std::string blockfilename = filename;
  //     std::stringstream ss;
  //     ss<< std::setfill('0') << std::setw(4) << idomain<< "_level_" << std::setfill('0') << std::setw(2) << level;
  //     blockfilename += "_block_"+ss.str()+".hpp5";
  //     const H5std_string h5filename(blockfilename);
  //     H5::H5File file( h5filename, H5F_ACC_TRUNC );
  //     mlmesh->getMesh(level)->writeH5(file);
  //     file.close();
  //   }
  // }
}
