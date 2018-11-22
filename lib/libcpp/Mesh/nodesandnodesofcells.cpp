#include  "Alat/sparsitypattern.hpp"
#include  "Mesh/nodesandnodesofcells.hpp"
#include  <cassert>
#include  <mpi.h>

using namespace mesh;

/*--------------------------------------------------------------------------*/
NodesAndNodesOfCells::~NodesAndNodesOfCells() {}
NodesAndNodesOfCells::NodesAndNodesOfCells(): GeometryObject(){}
NodesAndNodesOfCells::NodesAndNodesOfCells( const NodesAndNodesOfCells& nodesandnodesofcells): GeometryObject(nodesandnodesofcells)
{
  _nodes = nodesandnodesofcells._nodes;
  _nodes_of_cells = nodesandnodesofcells._nodes_of_cells;
}
NodesAndNodesOfCells& NodesAndNodesOfCells::operator=( const NodesAndNodesOfCells& nodesandnodesofcells)
{
  _nodes = nodesandnodesofcells._nodes;
  _nodes_of_cells = nodesandnodesofcells._nodes_of_cells;
  return *this;
}
std::string NodesAndNodesOfCells::getClassName() const
{
  return "NodesAndNodesOfCells";
}
std::unique_ptr<GeometryObject> NodesAndNodesOfCells::clone() const
{
  return std::unique_ptr<mesh::GeometryObject>(new NodesAndNodesOfCells(*this));
}

/*--------------------------------------------------------------------------*/
const arma::mat&  NodesAndNodesOfCells::getNodes() const {return _nodes;}
arma::mat&  NodesAndNodesOfCells::getNodes() {return _nodes;}

const alat::armaimat&  NodesAndNodesOfCells::getNodesOfCells() const {return _nodes_of_cells;}
alat::armaimat&  NodesAndNodesOfCells::getNodesOfCells() {return _nodes_of_cells;}

/*--------------------------------------------------------------------------*/
alat::armaivec NodesAndNodesOfCells::getSizes() const
{
  alat::armaivec sizes(4);
  sizes[0] = _nodes.n_rows;
  sizes[1] = _nodes.n_cols;
  sizes[2] = _nodes_of_cells.n_rows;
  sizes[3] = _nodes_of_cells.n_cols;
  // std::cerr << "NodesAndNodesOfCells sizes="<<sizes.t()<<"\n";
  return sizes;
}
void NodesAndNodesOfCells::setSizes(alat::armaivec::const_iterator sizes)
{
  // int nrows = sizes[0];
  // int ncols = sizes[1];
  // std::cerr << "nrows="<<nrows << " ncols="<<ncols<<"\n";
  // nrows = sizes[2];
  // ncols = sizes[3];
  // std::cerr << "nrows="<<nrows << " ncols="<<ncols<<"\n";
  _nodes.set_size(sizes[0], sizes[1]);
  _nodes_of_cells.set_size(sizes[2], sizes[3]);
}
void NodesAndNodesOfCells::send(int neighbor, int tag) const
{
  MPI_Request request;
  MPI_Isend(_nodes.begin(), _nodes.size(), MPI_DOUBLE, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Isend(_nodes_of_cells.begin(), _nodes_of_cells.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
}
void NodesAndNodesOfCells::recv(int neighbor, int tag)
{
  MPI_Status status;
  MPI_Request request;
  MPI_Irecv(_nodes.begin(), _nodes.size(), MPI_DOUBLE, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, &status);
  MPI_Irecv(_nodes_of_cells.begin(), _nodes_of_cells.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, &status);
}

/*--------------------------------------------------------------------------*/
void NodesAndNodesOfCells::loadH5(const arma::hdf5_name& spec)
{
  arma::hdf5_name h5name(spec.filename, spec.dsname+"/nodes", spec.opts);
  arma::hdf5_name h5name2(spec.filename, spec.dsname+"/nodes_of_cells", spec.opts);
  _nodes.load(h5name);
  _nodes_of_cells.load(h5name2);
}
void NodesAndNodesOfCells::saveH5(const arma::hdf5_name& spec) const
{
  arma::hdf5_name h5name(spec.filename, spec.dsname+"/nodes", spec.opts);
  arma::hdf5_name h5name2(spec.filename, spec.dsname+"/nodes_of_cells", spec.opts);
  _nodes.save(h5name);
  _nodes_of_cells.save(h5name2);
}
