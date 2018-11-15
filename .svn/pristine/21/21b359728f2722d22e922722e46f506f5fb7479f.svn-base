// #include "Alat/iomanager.hpp"
#include  "FadalightMesh/getmeshtype.hpp"
#include  "FadalightMesh/hexahedralmesh.hpp"
#include  "FadalightMesh/singlemeshcomposition.hpp"
// #include  "FadalightMesh/tetrahedralmesh.hpp"
#include  "FadalightMesh/quadrilateralmesh.hpp"
#include  "FadalightMesh/trianglemesh.hpp"
// #include  "FadalightMesh/linemesh.hpp"
// #include  "H5Cpp.hpp"
#include  <cassert>
#include  <fstream>
#include  <sstream>
#include  <iomanip>
#include  <sstream>

using namespace FadalightMesh;

/*--------------------------------------------------------------------------*/

SingleMeshComposition::~SingleMeshComposition()
{
  if(_mesh)
  {
    delete _mesh;
    _mesh = NULL;
  }
}

SingleMeshComposition::SingleMeshComposition() : FadalightMesh::MeshCompositionInterface()
{}
SingleMeshComposition::SingleMeshComposition( const SingleMeshComposition& singlemeshcomposition) : FadalightMesh::MeshCompositionInterface(singlemeshcomposition)
{
  assert(0);
}

SingleMeshComposition& SingleMeshComposition::operator=( const SingleMeshComposition& singlemeshcomposition)
{
  FadalightMesh::MeshCompositionInterface::operator=(singlemeshcomposition);
  assert(0);
  return *this;
}

std::string SingleMeshComposition::getClassName() const
{
  return "SingleMeshComposition";
}

SingleMeshComposition* SingleMeshComposition::clone() const
{
  return new SingleMeshComposition(*this);
}

/*--------------------------------------------------------------------------*/
int SingleMeshComposition::getDimension() const
{
  return _mesh->getDimension();
}

/*--------------------------------------------------------------------------*/

int SingleMeshComposition::getNCells() const
{
  return _mesh->getNCells();
}

/*--------------------------------------------------------------------------*/

int SingleMeshComposition::getNDomains() const
{
  return 1;
}

/*--------------------------------------------------------------------------*/

int SingleMeshComposition::getNCouplingMeshes() const
{
  return 0;
}

/*--------------------------------------------------------------------------*/

const FadalightMesh::MeshInterface* SingleMeshComposition::getMesh(int i) const
{
  return _mesh;
}

/*--------------------------------------------------------------------------*/

const FadalightMesh::CouplingMeshInterface* SingleMeshComposition::getCouplingMesh(int i) const
{
  return NULL;
}

/*--------------------------------------------------------------------------*/

void SingleMeshComposition::constructFadalightMesh(const std::string& meshname)
{
  alat::StringPair p = getMeshType(meshname);
  std::string type = p.first;
  std::string datatype = p.second;
  assert(datatype == "ascii" || datatype == "binary");

  // if(type=="FadalightMesh::LineMesh")
  // {
  //   _mesh = new FadalightMesh::LineMesh;
  // }

  // else
  if(type == "FadalightMesh::QuadrilateralMesh")
  {
    _mesh = new FadalightMesh::QuadrilateralMesh;
  }
  else if(type == "FadalightMesh::HexahedralMesh")
  {
    _mesh = new FadalightMesh::HexahedralMesh;
  }
  else if(type=="FadalightMesh::TriangleMesh")
  {
    _mesh = new FadalightMesh::TriangleMesh;
  }
  else
  {
    _error_string("constructFadalightMesh", "unknown meshtype: "+type);
  }
  // else if(type=="FadalightMesh::TetrahedralMesh")
  // {
  //   _mesh = new FadalightMesh::TetrahedralMesh;
  // }

  // ****
  // setGeometryObjectsConstructor(M);
  // delete _mesh->getGeometryObjectsConstructor();
  // _mesh->getGeometryObjectsConstructor()=new GeometryObjectsConstructor;
}

/*--------------------------------------------------------------------------*/

void SingleMeshComposition::read(const std::string& basefilename)
{
  _mesh->readFadalightMesh(basefilename);
}

/*--------------------------------------------------------------------------*/

void SingleMeshComposition::write(const std::string& basefilename, arma::file_type datatype) const
{
  _mesh->writeFadalightMesh(basefilename, datatype);
}

/*--------------------------------------------------------------------------*/

std::string SingleMeshComposition::getInfo() const
{
  std::stringstream ss1, ss2, ss3;
  ss1 << getNDomains();
  ss2 << getNCouplingMeshes();
  ss3 << getNCells();
  return ss1.str()+"_"+ss2.str()+"_"+ss3.str();
}

/*--------------------------------------------------------------------------*/

void SingleMeshComposition::writeMeshInfo(std::string filename, std::string blockfilename) const
{
  std::ofstream file( filename.c_str() );
  file <<_mesh->getDimension()<< "\n";
  int nlevels=1;
  file <<nlevels<< "\n";
  file << getNDomains() << "\n";
  file.close();
  int iblock = 0;
  // std::string blockfilename = alat::IoManager::getFileNameOnBlock( filename, iblock);
  _mesh->writeMeshInfo(blockfilename);
}

/*--------------------------------------------------------------------------*/

void SingleMeshComposition::writeH5(std::string filename) const
{
  assert(0);
  // int ndomains = getNDomains();
  // for(int idomain = 0; idomain < ndomains; idomain++)
  // {
  //   std::string blockfilename = filename;
  //   // char st[5];
  //   // sprintf(st, "%04d", idomain);
  //
  //   std::stringstream ss;
  //   ss<< std::setfill('0') << std::setw(4) << idomain;
  //   blockfilename += "_block_"+ss.str()+".hpp5";
  //   const H5std_string h5filename(blockfilename);
  //   H5::H5File file( h5filename, H5F_ACC_TRUNC );
  //   _mesh->writeH5(file);
  //   file.close();
  // }
}
