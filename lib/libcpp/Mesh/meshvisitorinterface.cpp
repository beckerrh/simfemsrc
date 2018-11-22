#include  "Mesh/meshvisitorinterface.hpp"
#include  <cassert>

using namespace mesh;

/*--------------------------------------------------------------------------*/
MeshVisitorInterface::~MeshVisitorInterface() {}
MeshVisitorInterface::MeshVisitorInterface(): alat::InterfaceBase(){}
MeshVisitorInterface::MeshVisitorInterface( const MeshVisitorInterface& meshvisitorinterface): alat::InterfaceBase(meshvisitorinterface) {}
MeshVisitorInterface& MeshVisitorInterface::operator=( const MeshVisitorInterface& meshvisitorinterface)
{
  assert(0);
  alat::InterfaceBase::operator=(meshvisitorinterface);
  return *this;
}
std::string MeshVisitorInterface::getClassName() const
{
  return "MeshVisitorInterface";
}
/*--------------------------------------------------------------------------*/
double MeshVisitorInterface::computeMeasureOfCell(int iK, const mesh::MeshUnitInterface* mesh) const {_notWritten("computeMeasureOfCell"); return 0.;}
void MeshVisitorInterface::computeNormal(MeshUnitInterface::Normal normal, const mesh::MeshUnitInterface* mesh, const alat::Node& xS, const MeshUnitInterface::Side S, int iK0, int ii0, int iK1, int ii1) const{_notWritten("computeNormal");}
