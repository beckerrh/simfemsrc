#include  "Alat/vecop.hpp"
#include  "Mesh/meshvisitortetrahedral.hpp"
#include  "Mesh/meshvisitortriangle.hpp"
#include  "Mesh/nodesandnodesofcells.hpp"
#include  <cassert>

using namespace mesh;

/*--------------------------------------------------------------------------*/
MeshVisitorTetrahedral::~MeshVisitorTetrahedral() {}
MeshVisitorTetrahedral::MeshVisitorTetrahedral(): MeshVisitorInterface(){}
MeshVisitorTetrahedral::MeshVisitorTetrahedral( const MeshVisitorTetrahedral& meshvisitortetrahedral): MeshVisitorInterface(meshvisitortetrahedral){}
MeshVisitorTetrahedral& MeshVisitorTetrahedral::operator=( const MeshVisitorTetrahedral& meshvisitortetrahedral)
{
  assert(0);
  MeshVisitorInterface::operator=(meshvisitortetrahedral);
  return *this;
}
std::string MeshVisitorTetrahedral::getClassName() const
{
  return "MeshVisitorTetrahedral";
}
std::unique_ptr<MeshVisitorInterface> MeshVisitorTetrahedral::clone() const
{
  return std::unique_ptr<mesh::MeshVisitorInterface>(new MeshVisitorTetrahedral(*this));
}

meshEnums::meshtype MeshVisitorTetrahedral::getType() const
{
  return meshEnums::TetrahedralMesh;
}
std::unique_ptr<MeshVisitorInterface> MeshVisitorTetrahedral::newBoundaryVisitor() const
{
  return std::unique_ptr<MeshVisitorInterface>(new MeshVisitorTriangle);
}

/*--------------------------------------------------------------------------*/
void MeshVisitorTetrahedral::set_dimensions(int& ndim, int& nnodes_per_cell, int& nedges_per_cell, int& nsides_per_cell, int& nnodes_per_side) const
{
  ndim = 3;
  nnodes_per_cell = 4;
  nsides_per_cell = 4;
  nnodes_per_side = 3;
  nedges_per_cell = 6;
}

/*-------------------------------------------------------*/
double MeshVisitorTetrahedral::computeMeasureOfCell(int iK, const mesh::MeshUnitInterface* mesh) const
{
  const arma::mat& nodes = mesh->getNodesAndNodesOfCells().getNodes();
  const MeshUnitInterface::Cell K = mesh->getCell(iK);
  double dx1 = nodes(0,K[1])-nodes(0,K[0]);
  double dx2 = nodes(0,K[2])-nodes(0,K[0]);
  double dx3 = nodes(0,K[3])-nodes(0,K[0]);
  double dy1 = nodes(1,K[1])-nodes(1,K[0]);
  double dy2 = nodes(1,K[2])-nodes(1,K[0]);
  double dy3 = nodes(1,K[3])-nodes(1,K[0]);
  double dz1 = nodes(2,K[1])-nodes(2,K[0]);
  double dz2 = nodes(2,K[2])-nodes(2,K[0]);
  double dz3 = nodes(2,K[3])-nodes(2,K[0]);
  return ( dx1*( dy2*dz3-dz2*dy3 )-dy1*( dx2*dz3-dz2*dx3 )+dz1*( dx2*dy3-dx3*dy2 ) )/6.0;
}
/*---------------------------------------------------------*/
void MeshVisitorTetrahedral::computeNormal(MeshUnitInterface::Normal normal, const mesh::MeshUnitInterface* mesh, const alat::Node& xS, const MeshUnitInterface::Side S, int iK0, int ii0, int iK1, int ii1) const
{
  const arma::mat& nodes = mesh->getNodesAndNodesOfCells().getNodes();
  const alat::armaimat& nodes_of_cells = mesh->getNodesAndNodesOfCells().getNodesOfCells();
  alat::Node u, v, w;
  u.x() = nodes(0,S[1])-nodes(0,S[0]);
  u.y() = nodes(1,S[1])-nodes(1,S[0]);
  u.z() = nodes(2,S[1])-nodes(2,S[0]);
  v.x() = nodes(0,S[2])-nodes(0,S[0]);
  v.y() = nodes(1,S[2])-nodes(1,S[0]);
  v.z() = nodes(2,S[2])-nodes(2,S[0]);
  alat::vectorProduct(w,v,u);

  alat::Node x0 = mesh->getNodeOfCell(iK0);
  x0 -= xS;
  double test = arma::dot(x0,w);
  if(test<0.0)
  {
    normal = 0.5*w;
  }
  else
  {
    normal = -0.5*w;
 }
}
