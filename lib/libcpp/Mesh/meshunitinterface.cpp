#include  "Mesh/meshunitinterface.hpp"
#include  <cassert>

using namespace mesh;

/*--------------------------------------------------------------------------*/
MeshUnitInterface::~MeshUnitInterface(){}
MeshUnitInterface::MeshUnitInterface() : alat::InterfaceBase(){}
MeshUnitInterface::MeshUnitInterface( const MeshUnitInterface& meshinterface) : alat::InterfaceBase(meshinterface){}
MeshUnitInterface& MeshUnitInterface::operator=( const MeshUnitInterface& meshinterface)
{
  InterfaceBase::operator=(meshinterface);
  assert(0);
  return *this;
}
std::string MeshUnitInterface::getInterfaceName() const
{
  return "MeshUnitInterface";
}
/*--------------------------------------------------------------------------*/
alat::armaivec MeshUnitInterface::_getSideOfCell(int i, int ii) const
{
  _notWritten("_getSideOfCell");
  return alat::armaivec(-1);
}
alat::armaivec MeshUnitInterface::_getEdgeOfCell(int i, int ii) const
{
  _notWritten("_getEdgeOfCell");
  return alat::armaivec(-1);
}

/*--------------------------------------------------------------------------*/
void MeshUnitInterface::writeVtk(std::string filename) const
{
  _notWritten("writeVtk");
}
void MeshUnitInterface::writeBoundaryVtk(std::string filename) const
{
  _notWritten("writeBoundaryVtk");
}
void MeshUnitInterface::readGmsh(std::string filename)
{
  _notWritten("readGmsh");
}

/*--------------------------------------------------------------------------*/
// double MeshUnitInterface::computeMeasureOfCell(int iK) const
// {
//   _notWritten("computeMeasureOfCell");
//   return 0.0;
// }
// void MeshUnitInterface::computeNormal(Normal normal, const alat::Node& xS, const Side S, int iK, int ii, int iK1, int ii1) const
// {
//   _notWritten("computeNormal");
// }
