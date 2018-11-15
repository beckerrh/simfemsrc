#include  "Alat/vecop.hpp"
#include  "Mesh/meshvisitortriangle.hpp"
#include  "Mesh/meshvisitorline.hpp"
#include  "Mesh/nodesandnodesofcells.hpp"
#include  <cassert>

using namespace mesh;

/*--------------------------------------------------------------------------*/
MeshVisitorTriangle::~MeshVisitorTriangle() {}
MeshVisitorTriangle::MeshVisitorTriangle(): MeshVisitorInterface(){}
MeshVisitorTriangle::MeshVisitorTriangle( const MeshVisitorTriangle& meshvisitortriangle): MeshVisitorInterface(meshvisitortriangle){}
MeshVisitorTriangle& MeshVisitorTriangle::operator=( const MeshVisitorTriangle& meshvisitortriangle)
{
  assert(0);
  MeshVisitorInterface::operator=(meshvisitortriangle);
  return *this;
}
std::string MeshVisitorTriangle::getClassName() const
{
  return "MeshVisitorTriangle";
}
std::unique_ptr<MeshVisitorInterface> MeshVisitorTriangle::clone() const
{
  // std::cerr << "MeshVisitorTriangle::clone()\n";
  return std::unique_ptr<mesh::MeshVisitorInterface>(new MeshVisitorTriangle(*this));
}

/*--------------------------------------------------------------------------*/
void MeshVisitorTriangle::set_dimensions(int& ndim, int& nnodes_per_cell, int& nedges_per_cell, int& nsides_per_cell, int& nnodes_per_side) const
{
  ndim = 2;
  nnodes_per_cell = 3;
  nsides_per_cell = 3;
  nnodes_per_side = 2;
  nedges_per_cell = 3;
}
meshEnums::meshtype MeshVisitorTriangle::getType() const
{
  return meshEnums::TriangleMesh;
}
std::unique_ptr<MeshVisitorInterface> MeshVisitorTriangle::newBoundaryVisitor() const
{
  return std::unique_ptr<MeshVisitorInterface>(new MeshVisitorLine);
}

/*-------------------------------------------------------*/
double MeshVisitorTriangle::computeMeasureOfCell(int iK, const mesh::MeshUnitInterface* mesh) const
{
  const arma::mat& nodes = mesh->getNodesAndNodesOfCells().getNodes();
  const MeshUnitInterface::Cell K = mesh->getCell(iK);
  alat::Node u, v, w;
  u.x() = nodes(0,K[1])-nodes(0,K[0]);
  u.y() = nodes(1,K[1])-nodes(1,K[0]);
  u.z() = nodes(2,K[1])-nodes(2,K[0]);
  v.x() = nodes(0,K[2])-nodes(0,K[0]);
  v.y() = nodes(1,K[2])-nodes(1,K[0]);
  v.z() = nodes(2,K[2])-nodes(2,K[0]);
  alat::vectorProduct(w,u,v);
  return 0.5*w.norm();
}
/*---------------------------------------------------------*/
void MeshVisitorTriangle::computeNormal(MeshUnitInterface::Normal normal, const mesh::MeshUnitInterface* mesh, const alat::Node& xS, const MeshUnitInterface::Side S, int iK0, int ii0, int iK1, int ii1) const
{
  const arma::mat& nodes = mesh->getNodesAndNodesOfCells().getNodes();
  const alat::armaimat& nodes_of_cells = mesh->getNodesAndNodesOfCells().getNodesOfCells();

  alat::Node u;
  u.x() = nodes(0,S[1])-nodes(0,S[0]);
  u.y() = nodes(1,S[1])-nodes(1,S[0]);
  u.z() = nodes(2,S[1])-nodes(2,S[0]);
  double unorm = u.norm();

  alat::Node t0, t1;
  t0.x() = nodes(0, S[0]) - nodes(0, nodes_of_cells(ii0,iK0));
  t0.y() = nodes(1, S[0]) - nodes(1, nodes_of_cells(ii0,iK0));
  t0.z() = nodes(2, S[0]) - nodes(2, nodes_of_cells(ii0,iK0));
  t1.x() = nodes(0, S[1]) - nodes(0, nodes_of_cells(ii0,iK0));
  t1.y() = nodes(1, S[1]) - nodes(1, nodes_of_cells(ii0,iK0));
  t1.z() = nodes(2, S[1]) - nodes(2, nodes_of_cells(ii0,iK0));

  // n0 = t0 - a t1
  // 0 = t0*u - a t1*u
  double a0 = t0*u;
  double a1 = t1*u;
  t0.add(-a0/a1, t1);
  double t0norm = t0.norm();
  assert(t0norm>0.0);

  alat::Node x0 = mesh->getNodeOfCell(iK0);
  x0 -= xS;
  double test = x0*t0;
  if(test<0.0)
  {
    t0.scale(unorm/t0norm);
  }
  else
  {
    t0.scale(-unorm/t0norm);
  }
  // std::cerr << "t0 " << t0 << "\n";
  if(iK1==-1)
  {
    normal = t0;
    return;
  }

  alat::Node t2, t3;
  t2.x() = nodes(0, S[0]) - nodes(0, nodes_of_cells(ii1,iK1));
  t2.y() = nodes(1, S[0]) - nodes(1, nodes_of_cells(ii1,iK1));
  t2.z() = nodes(2, S[0]) - nodes(2, nodes_of_cells(ii1,iK1));
  t3.x() = nodes(0, S[1]) - nodes(0, nodes_of_cells(ii1,iK1));
  t3.y() = nodes(1, S[1]) - nodes(1, nodes_of_cells(ii1,iK1));
  t3.z() = nodes(2, S[1]) - nodes(2, nodes_of_cells(ii1,iK1));

  a0 = t2*u;
  a1 = t3*u;

  t2.add(-a0/a1, t3);
  double t2norm = t2.norm();
  alat::Node x1 = mesh->getNodeOfCell(iK1);
  x1 -= xS;
  test = x1*t2;
  if(test<0.0)
  {
    t2.scale(-unorm/t2norm);
  }
  else
  {
    t2.scale(unorm/t2norm);
  }
  // std::cerr << "t2 " << t2 << "\n";

  normal = 0.5*t0 + 0.5*t2;
  // std::cerr << "u " << u << " normal " << normal << "\n";
}
