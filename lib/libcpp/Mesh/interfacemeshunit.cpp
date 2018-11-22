#include  "Mesh/interfacemeshunit.hpp"
#include  <cassert>

using namespace mesh;

/*--------------------------------------------------------------------------*/
InterfaceMeshUnit::~InterfaceMeshUnit() {}
InterfaceMeshUnit::InterfaceMeshUnit(): MeshUnit(), _plainmesh(NULL) {}
InterfaceMeshUnit::InterfaceMeshUnit( const InterfaceMeshUnit& interfacemeshunit): MeshUnit(interfacemeshunit)
{
  _plainmesh = interfacemeshunit._plainmesh;
  _interfacemeshinfo = interfacemeshunit._interfacemeshinfo;
}
void InterfaceMeshUnit::init(std::shared_ptr<const mesh::MeshVisitorInterface> visitor, const MeshUnitInterface* plainmesh)
{
  _plainmesh = plainmesh;
  MeshUnit::init(visitor);
}
InterfaceMeshUnit& InterfaceMeshUnit::operator=( const InterfaceMeshUnit& interfacemeshunit)
{
  assert(0);
  MeshUnit::operator=(interfacemeshunit);
  return *this;
}
std::string InterfaceMeshUnit::getClassName() const
{
  return "InterfaceMeshUnit";
}
/*--------------------------------------------------------------------------*/
const MeshUnitInterface* InterfaceMeshUnit::getPlainMesh() const {return _plainmesh;}
const InterfaceMeshInfo& InterfaceMeshUnit::getInterfaceMeshInfo() const {return _interfacemeshinfo;}
InterfaceMeshInfo& InterfaceMeshUnit::getInterfaceMeshInfo() {return _interfacemeshinfo;}

/*---------------------------------------------------------*/
void InterfaceMeshUnit::loadH5(const arma::hdf5_name& spec)
{
  MeshUnit::loadH5(spec);
  // _parentmesh = ???;
  _interfacemeshinfo.loadH5(arma::hdf5_name(spec.filename, spec.dsname+"/"+_interfacemeshinfo.getClassName(), spec.opts));
}

/*---------------------------------------------------------*/
void InterfaceMeshUnit::saveH5(const arma::hdf5_name& spec) const
{
  MeshUnit::saveH5(spec);
  // _parentmesh = ???;
  _interfacemeshinfo.saveH5(arma::hdf5_name(spec.filename, spec.dsname+"/"+_interfacemeshinfo.getClassName(), spec.opts));
}
