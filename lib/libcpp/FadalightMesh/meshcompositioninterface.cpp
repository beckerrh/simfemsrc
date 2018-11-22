#include  "FadalightMesh/meshcompositioninterface.hpp"
#include  <cassert>

using namespace FadalightMesh;

/*--------------------------------------------------------------------------*/

MeshCompositionInterface::~MeshCompositionInterface()
{}

/*--------------------------------------------------------------------------*/

MeshCompositionInterface::MeshCompositionInterface() : alat::InterfaceBase()
{}

/*--------------------------------------------------------------------------*/

MeshCompositionInterface::MeshCompositionInterface( const MeshCompositionInterface& meshcompositioninterface) : alat::InterfaceBase(meshcompositioninterface)
{
  assert(0);
}

/*--------------------------------------------------------------------------*/

MeshCompositionInterface& MeshCompositionInterface::operator=( const MeshCompositionInterface& meshcompositioninterface)
{
  InterfaceBase::operator=(meshcompositioninterface);
  assert(0);
  return *this;
}

/*--------------------------------------------------------------------------*/

std::string MeshCompositionInterface::getInterfaceName() const
{
  return "MeshCompositionInterface";
}

/*--------------------------------------------------------------------------*/

std::string MeshCompositionInterface::getClassName() const
{
  return "MeshCompositionInterface";
}

/*--------------------------------------------------------------------------*/

const FadalightMesh::MeshInterface* MeshCompositionInterface::getMacroMesh() const
{
  _notWritten("getMacroMesh");
  return NULL;
}

/*--------------------------------------------------------------------------*/

MeshCompositionInterface* MeshCompositionInterface::clone() const
{
  assert(0);
// return new MeshCompositionInterface(*this);
}

/*--------------------------------------------------------------------------*/

void MeshCompositionInterface::writeMeshInfo(std::string filename, std::string blockfilename) const
{
  _notWritten("writeMeshInfo");
}

// /*--------------------------------------------------------------------------*/
//
// void MeshCompositionInterface::initCouplingGrids(const alat::armavecModelManagerInterface* _modelmanager)
// {
// _notWritten("initCouplingGrids");
// }

/*--------------------------------------------------------------------------*/

void MeshCompositionInterface::constructFadalightMesh(const std::string& meshname)
{
  _notWritten("constructFadalightMesh");
}
