#include  "Mesh/interfacemeshunitinterface.hpp"
#include  <cassert>

using namespace mesh;

/*--------------------------------------------------------------------------*/
InterfaceMeshUnitInterface::~InterfaceMeshUnitInterface() {}
InterfaceMeshUnitInterface::InterfaceMeshUnitInterface(): MeshUnitInterface(){}
InterfaceMeshUnitInterface::InterfaceMeshUnitInterface( const InterfaceMeshUnitInterface& interfacemeshunitinterface): MeshUnitInterface(interfacemeshunitinterface)
{
  (*this).operator=(interfacemeshunitinterface);
}
InterfaceMeshUnitInterface& InterfaceMeshUnitInterface::operator=( const InterfaceMeshUnitInterface& interfacemeshunitinterface)
{
  assert(0);
  MeshUnitInterface::operator=(interfacemeshunitinterface);
  return *this;
}
std::string InterfaceMeshUnitInterface::getClassName() const
{
  return "InterfaceMeshUnitInterface";
}
/*--------------------------------------------------------------------------*/
