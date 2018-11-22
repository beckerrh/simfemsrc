#include  "Mesh/boundarymeshunit.hpp"
#include  <cassert>

using namespace mesh;

/*--------------------------------------------------------------------------*/
BoundaryMeshUnit::~BoundaryMeshUnit() {}
BoundaryMeshUnit::BoundaryMeshUnit(): MeshUnit(), _parentmesh(NULL){}
BoundaryMeshUnit::BoundaryMeshUnit( const BoundaryMeshUnit& boundarymeshunit): MeshUnit(boundarymeshunit)
{
  _parentmesh = boundarymeshunit._parentmesh;
  _boundarymeshinfo = boundarymeshunit._boundarymeshinfo;
}
void BoundaryMeshUnit::init(std::shared_ptr<const mesh::MeshVisitorInterface> visitor, const MeshUnitInterface* parentmesh)
{
  _parentmesh = parentmesh;
  MeshUnit::init(visitor);
}
BoundaryMeshUnit& BoundaryMeshUnit::operator=( const BoundaryMeshUnit& boundarymeshunit)
{
  assert(0);
  MeshUnit::operator=(boundarymeshunit);
  return *this;
}
std::string BoundaryMeshUnit::getClassName() const
{
  return "BoundaryMeshUnit";
}
/*--------------------------------------------------------------------------*/
const MeshUnitInterface* BoundaryMeshUnit::getParentMesh() const {return _parentmesh;}
const BoundaryMeshInfo& BoundaryMeshUnit::getBoundaryMeshInfo() const {return _boundarymeshinfo;}
BoundaryMeshInfo& BoundaryMeshUnit::getBoundaryMeshInfo() {return _boundarymeshinfo;}

/*---------------------------------------------------------*/
void BoundaryMeshUnit::loadH5(const arma::hdf5_name& spec)
{
  MeshUnit::loadH5(spec);
  // _parentmesh = ???;
  _boundarymeshinfo.loadH5(arma::hdf5_name(spec.filename, spec.dsname+"/"+_boundarymeshinfo.getClassName(), spec.opts));
}

/*---------------------------------------------------------*/
void BoundaryMeshUnit::saveH5(const arma::hdf5_name& spec) const
{
  MeshUnit::saveH5(spec);
  // _parentmesh = ???;
  _boundarymeshinfo.saveH5(arma::hdf5_name(spec.filename, spec.dsname+"/"+_boundarymeshinfo.getClassName(), spec.opts));
}
