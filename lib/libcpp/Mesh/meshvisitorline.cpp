#include  "Mesh/meshvisitorline.hpp"
#include  "Mesh/meshvisitorpoint.hpp"
#include  "Mesh/nodesandnodesofcells.hpp"
#include  <cassert>

using namespace mesh;

/*--------------------------------------------------------------------------*/
MeshVisitorLine::~MeshVisitorLine() {}
MeshVisitorLine::MeshVisitorLine(): MeshVisitorInterface(){}
MeshVisitorLine::MeshVisitorLine( const MeshVisitorLine& meshvisitorline): MeshVisitorInterface(meshvisitorline){}
MeshVisitorLine& MeshVisitorLine::operator=( const MeshVisitorLine& meshvisitorline)
{
  assert(0);
  MeshVisitorInterface::operator=(meshvisitorline);
  return *this;
}
std::string MeshVisitorLine::getClassName() const
{
  return "MeshVisitorLine";
}
meshEnums::meshtype MeshVisitorLine::getType() const
{
  return meshEnums::LineMesh;
}
std::unique_ptr<MeshVisitorInterface> MeshVisitorLine::clone() const
{
  return std::unique_ptr<mesh::MeshVisitorInterface>(new MeshVisitorLine(*this));
}

/*--------------------------------------------------------------------------*/
void MeshVisitorLine::set_dimensions(int& ndim, int& nnodes_per_cell, int& nedges_per_cell, int& nsides_per_cell, int& nnodes_per_side) const
{
  ndim = 1;
  nnodes_per_cell = 2;
  nsides_per_cell = 2;
  nnodes_per_side = 1;
  nedges_per_cell = 1;
}
std::unique_ptr<MeshVisitorInterface> MeshVisitorLine::newBoundaryVisitor() const
{
  return std::unique_ptr<MeshVisitorInterface>(new MeshVisitorPoint);
}

/*-------------------------------------------------------*/
double MeshVisitorLine::computeMeasureOfCell(int iK, const mesh::MeshUnitInterface* mesh) const
{
  // const arma::mat& nodes = MeshUnit::nodesandnodesofcells.nodes;
  const arma::mat& nodes = mesh->getNodesAndNodesOfCells().getNodes();
  // const Tetrahedron& K = getCell(iK);
  const MeshUnitInterface::Cell K = mesh->getCell(iK);
  alat::armavec dx(3);
  dx[0] = nodes(0,K[1])-nodes(0,K[0]);
  dx[1] = nodes(1,K[1])-nodes(1,K[0]);
  dx[2] = nodes(2,K[1])-nodes(2,K[0]);
  return arma::norm(dx);
}
/*---------------------------------------------------------*/
void MeshVisitorLine::computeNormal(MeshUnitInterface::Normal normal, const mesh::MeshUnitInterface* mesh, const alat::Node& xS, const MeshUnitInterface::Side S, int iK0, int ii0, int iK1, int ii1) const
{
  const arma::mat& nodes = mesh->getNodesAndNodesOfCells().getNodes();
  const alat::armaimat& nodes_of_cells = mesh->getNodesAndNodesOfCells().getNodesOfCells();
  // const arma::mat& nodes = MeshUnit::nodesandnodesofcells.nodes;
  // const alat::armaimat& nodes_of_cells = MeshUnit::nodesandnodesofcells.nodes_of_cells;
  alat::armavec dx0(3);
  const MeshUnitInterface::Cell K0 = mesh->getCell(iK0);
  dx0[0] = nodes(0,K0[1])-nodes(0,K0[0]);
  dx0[1] = nodes(1,K0[1])-nodes(1,K0[0]);
  dx0[2] = nodes(2,K0[1])-nodes(2,K0[0]);
  double d0 = 1.0/arma::norm(dx0);
  dx0 *= d0;

  alat::Node x0 = mesh->getNodeOfCell(iK0);
  x0 -= xS;
  double test = arma::dot(x0,dx0);
  if(test>0.0)
  {
    dx0 *= -1.0;
	}
  if(iK1==-1)
  {
  	normal = dx0;
  	return;
  }
  alat::armavec dx1(3);
  const MeshUnitInterface::Cell K1 = mesh->getCell(iK1);
  dx1[0] = nodes(0,K1[1])-nodes(0,K1[0]);
  dx1[1] = nodes(1,K1[1])-nodes(1,K1[0]);
  dx1[2] = nodes(2,K1[1])-nodes(2,K1[0]);
  double d1 = 1.0/arma::norm(dx1);
  dx1 *= d1;
  alat::Node x1 = mesh->getNodeOfCell(iK1);
  x1 -= xS;
  test = arma::dot(x1,dx1);
  if(test<0.0)
  {
    dx1 *= -1.0;
	}
  normal = 0.5*dx0+0.5*dx1;
  // std::cerr << "xS="<<xS << " S="<<S << " iK0="<<iK0 << " iK1="<<iK1<< " --> dx0="<<dx0<< " --> dx1="<<dx1<< " ----> normal="<<normal<<"\n";
}
