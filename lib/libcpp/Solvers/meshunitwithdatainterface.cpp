#include  "Solvers/meshunitwithdatainterface.hpp"
#include  "Solvers/solverinterface.hpp"
#include  <cassert>

using namespace solvers;

/*--------------------------------------------------------------------------*/
MeshInfo::MeshInfo(const arma::mat& in_nodes, const alat::armaimat& in_nodes_of_cells, const alat::armaimat& in_sides_of_cells, const alat::armaimat& in_cells_of_sides, const alat::armaimat& in_edges_of_cells, const alat::armaimat& in_nodes_of_sides, const alat::armaimat& in_nodes_of_edges, const alat::armaicube& in_localnodes_of_edges_of_cell, const arma::vec& in_measure_of_cells, const arma::mat& in_normals, const arma::fmat& in_sigma, double in_dim, int in_nnodes, int in_ncells, int in_nsides, int in_nedges, int in_nnodespercell, int in_nedgespercell, int in_nsidespercell, int in_nnodesperside, const mesh::MeshUnitInterface::BoundaryInformationMap& in_bdrymesheunitsmap) : nodes(in_nodes),  nodes_of_cells(in_nodes_of_cells), sides_of_cells(in_sides_of_cells), cells_of_sides(in_cells_of_sides), edges_of_cells(in_edges_of_cells), nodes_of_sides(in_nodes_of_sides), nodes_of_edges(in_nodes_of_edges), localnodes_of_edges_of_cell(in_localnodes_of_edges_of_cell), measure_of_cells(in_measure_of_cells), normals(in_normals), sigma(in_sigma), dim(in_dim), nnodes(in_nnodes), ncells(in_ncells), nsides(in_nsides), nedges(in_nedges), nnodespercell(in_nnodespercell), nedgespercell(in_nedgespercell), nsidespercell(in_nsidespercell), nnodesperside(in_nnodesperside), bdrymesheunitsmap(in_bdrymesheunitsmap){}

/*--------------------------------------------------------------------------*/
MeshUnitWithDataInterface::~MeshUnitWithDataInterface() {}
MeshUnitWithDataInterface::MeshUnitWithDataInterface(): alat::InterfaceBase(){}
MeshUnitWithDataInterface::MeshUnitWithDataInterface( const MeshUnitWithDataInterface& meshunitwithdatainterface): alat::InterfaceBase(meshunitwithdatainterface)
{
assert(0);
}
MeshUnitWithDataInterface& MeshUnitWithDataInterface::operator=( const MeshUnitWithDataInterface& meshunitwithdatainterface)
{
  assert(0);
  alat::InterfaceBase::operator=(meshunitwithdatainterface);
  return *this;
}
std::string MeshUnitWithDataInterface::getClassName() const
{
  return "MeshUnitWithDataInterface";
}
/*--------------------------------------------------------------------------*/
