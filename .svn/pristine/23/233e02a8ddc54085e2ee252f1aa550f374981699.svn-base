#include  "Mesh/meshvisitorpoint.hpp"
#include  <cassert>

using namespace mesh;

/*--------------------------------------------------------------------------*/
MeshVisitorPoint::~MeshVisitorPoint() {}
MeshVisitorPoint::MeshVisitorPoint(): MeshVisitorInterface(){}
MeshVisitorPoint::MeshVisitorPoint( const MeshVisitorPoint& meshvisitorpoint): MeshVisitorInterface(meshvisitorpoint){}
MeshVisitorPoint& MeshVisitorPoint::operator=( const MeshVisitorPoint& meshvisitorpoint)
{
  assert(0);
  MeshVisitorInterface::operator=(meshvisitorpoint);
  return *this;
}
std::string MeshVisitorPoint::getClassName() const
{
  return "MeshVisitorPoint";
}
meshEnums::meshtype MeshVisitorPoint::getType() const
{
  return meshEnums::PointMesh;
}
std::unique_ptr<MeshVisitorInterface> MeshVisitorPoint::newBoundaryVisitor() const
{
  assert(0);
  return std::unique_ptr<MeshVisitorInterface>(nullptr);
}
std::unique_ptr<MeshVisitorInterface> MeshVisitorPoint::clone() const
{
  return std::unique_ptr<mesh::MeshVisitorInterface>(new MeshVisitorPoint(*this));
}

/*--------------------------------------------------------------------------*/
void MeshVisitorPoint::set_dimensions(int& ndim, int& nnodes_per_cell, int& nedges_per_cell, int& nsides_per_cell, int& nnodes_per_side) const
{
  ndim = 0;
  nnodes_per_cell = 1;
  nsides_per_cell = 1;
  nnodes_per_side = 1;
  nedges_per_cell = 0;
}
