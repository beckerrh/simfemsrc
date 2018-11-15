#include  "Mesh/meshinterface.hpp"
#include  <cassert>

using namespace mesh;

/*--------------------------------------------------------------------------*/
MeshInterface::~MeshInterface(){}
MeshInterface::MeshInterface() : alat::InterfaceBase(){}
MeshInterface::MeshInterface( const MeshInterface& meshinterface) : alat::InterfaceBase(meshinterface){}
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
void MeshInterface::writeVtk(std::string filename) const
{
  _notWritten("writeVtk");
}
void MeshInterface::writeBoundaryVtk(std::string filename) const
{
  _notWritten("writeBoundaryVtk");
}
void MeshInterface::readGmsh(std::string filename)
{
  _notWritten("readGmsh");
}
