#include  "Mesh/edgesandcells.hpp"
#include  <cassert>
#include  <mpi.h>

using namespace mesh;

/*--------------------------------------------------------------------------*/
EdgesAndCells::~EdgesAndCells() {}
EdgesAndCells::EdgesAndCells(): GeometryObject(){}
EdgesAndCells::EdgesAndCells( const EdgesAndCells& edgesandcells): GeometryObject(edgesandcells)
{
  _nodes_of_edges = edgesandcells._nodes_of_edges;
  _edges_of_cells = edgesandcells._edges_of_cells;
}
EdgesAndCells& EdgesAndCells::operator=( const EdgesAndCells& edgesandcells)
{
  _nodes_of_edges = edgesandcells._nodes_of_edges;
  _edges_of_cells = edgesandcells._edges_of_cells;
  return *this;
}
std::string EdgesAndCells::getClassName() const
{
  return "EdgesAndCells";
}
std::unique_ptr<GeometryObject> EdgesAndCells::clone() const
{
  return std::unique_ptr<mesh::GeometryObject>(new EdgesAndCells(*this));
}

/*--------------------------------------------------------------------------*/
const alat::armaimat&  EdgesAndCells::getNodesOfEdges() const {return _nodes_of_edges;}
alat::armaimat&  EdgesAndCells::getNodesOfEdges() {return _nodes_of_edges;}
const alat::armaimat&  EdgesAndCells::getEdgesOfCells() const {return _edges_of_cells;}
alat::armaimat&  EdgesAndCells::getEdgesOfCells() {return _edges_of_cells;}
const alat::armaicube& EdgesAndCells::getLocalnodesOfEdgesOfCells() const{return _localnodes_of_edges_in_cells;}
alat::armaicube& EdgesAndCells::getLocalnodesOfEdgesOfCells(){return _localnodes_of_edges_in_cells;}

/*--------------------------------------------------------------------------*/
alat::armaivec EdgesAndCells::getSizes() const
{
  alat::armaivec sizes(7);
  sizes[0] = _nodes_of_edges.n_rows;
  sizes[1] = _nodes_of_edges.n_cols;
  sizes[2] = _edges_of_cells.n_rows;
  sizes[3] = _edges_of_cells.n_cols;
  sizes[4] = _localnodes_of_edges_in_cells.n_rows;
  sizes[5] = _localnodes_of_edges_in_cells.n_cols;
  sizes[6] = _localnodes_of_edges_in_cells.n_slices;
  // std::cerr << "EdgesAndCells sizes="<<sizes.t()<<"\n";
  return sizes;
}
void EdgesAndCells::setSizes(alat::armaivec::const_iterator sizes)
{
  _nodes_of_edges.set_size(sizes[0], sizes[1]);
  _edges_of_cells.set_size(sizes[2], sizes[3]);
  _localnodes_of_edges_in_cells.set_size(sizes[6], sizes[7], sizes[8]);
}
void EdgesAndCells::send(int neighbor, int tag) const
{
  MPI_Request request;
  MPI_Isend(_nodes_of_edges.begin(), _nodes_of_edges.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Isend(_edges_of_cells.begin(), _edges_of_cells.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Isend(_localnodes_of_edges_in_cells.begin(), _localnodes_of_edges_in_cells.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
}
void EdgesAndCells::recv(int neighbor, int tag)
{
  MPI_Status status;
  MPI_Request request;
  MPI_Irecv(_nodes_of_edges.begin(), _nodes_of_edges.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, &status);
  MPI_Irecv(_edges_of_cells.begin(), _edges_of_cells.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, &status);
  MPI_Irecv(_localnodes_of_edges_in_cells.begin(), _localnodes_of_edges_in_cells.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, &status);
}

/*--------------------------------------------------------------------------*/
void EdgesAndCells::loadH5(const arma::hdf5_name& spec)
{
  _nodes_of_edges.load(arma::hdf5_name(spec.filename, spec.dsname+"/nodes_of_edges", spec.opts));
  _edges_of_cells.load(arma::hdf5_name(spec.filename, spec.dsname+"/edges_of_cells", spec.opts));
  _localnodes_of_edges_in_cells.load(arma::hdf5_name(spec.filename, spec.dsname+"/localnodes_of_edges_in_cells", spec.opts));
}
void EdgesAndCells::saveH5(const arma::hdf5_name& spec) const
{
  _nodes_of_edges.save(arma::hdf5_name(spec.filename, spec.dsname+"/nodes_of_edges", spec.opts));
  _edges_of_cells.save(arma::hdf5_name(spec.filename, spec.dsname+"/edges_of_cells", spec.opts));
  _localnodes_of_edges_in_cells.save(arma::hdf5_name(spec.filename, spec.dsname+"/localnodes_of_edges_in_cells", spec.opts));
}
