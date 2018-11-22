#include  "Mesh/sidesandcells.hpp"
#include  <cassert>
#include  <mpi.h>

using namespace mesh;

/*--------------------------------------------------------------------------*/
SidesAndCells::~SidesAndCells() {}
SidesAndCells::SidesAndCells(): GeometryObject(){}
SidesAndCells::SidesAndCells( const SidesAndCells& sidesandcells): GeometryObject(sidesandcells)
{
  _nodes_of_sides = sidesandcells._nodes_of_sides;
  _sides_of_cells = sidesandcells._sides_of_cells;
  _cells_of_sides = sidesandcells._cells_of_sides;
}
SidesAndCells& SidesAndCells::operator=( const SidesAndCells& sidesandcells)
{
  _nodes_of_sides = sidesandcells._nodes_of_sides;
  _sides_of_cells = sidesandcells._sides_of_cells;
  _cells_of_sides = sidesandcells._cells_of_sides;
  return *this;
}
std::string SidesAndCells::getClassName() const
{
  return "SidesAndCells";
}
std::unique_ptr<GeometryObject> SidesAndCells::clone() const
{
  return std::unique_ptr<mesh::GeometryObject>(new SidesAndCells(*this));
}

/*--------------------------------------------------------------------------*/
const alat::armaimat&  SidesAndCells::getNodesOfSides() const {return _nodes_of_sides;}
alat::armaimat&  SidesAndCells::getNodesOfSides() {return _nodes_of_sides;}
const alat::armaimat&  SidesAndCells::getSidesOfCells() const {return _sides_of_cells;}
alat::armaimat&  SidesAndCells::getSidesOfCells() {return _sides_of_cells;}
const alat::armaimat&  SidesAndCells::getCellsOfSides() const {return _cells_of_sides;}
alat::armaimat&  SidesAndCells::getCellsOfSides() {return _cells_of_sides;}

/*--------------------------------------------------------------------------*/
alat::armaivec SidesAndCells::getSizes() const
{
  alat::armaivec sizes(6);
  sizes[0] = _nodes_of_sides.n_rows;
  sizes[1] = _nodes_of_sides.n_cols;
  sizes[2] = _sides_of_cells.n_rows;
  sizes[3] = _sides_of_cells.n_cols;
  sizes[4] = _cells_of_sides.n_rows;
  sizes[5] = _cells_of_sides.n_cols;
  // std::cerr << "SidesAndCells sizes="<<sizes.t()<<"\n";
  return sizes;
}
void SidesAndCells::setSizes(alat::armaivec::const_iterator sizes)
{
  _nodes_of_sides.set_size(sizes[0], sizes[1]);
  _sides_of_cells.set_size(sizes[2], sizes[3]);
  _cells_of_sides.set_size(sizes[4], sizes[5]);
}
void SidesAndCells::send(int neighbor, int tag) const
{
  MPI_Request request;
  MPI_Isend(_nodes_of_sides.begin(), _nodes_of_sides.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Isend(_sides_of_cells.begin(), _sides_of_cells.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Isend(_cells_of_sides.begin(), _cells_of_sides.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
}
void SidesAndCells::recv(int neighbor, int tag)
{
  MPI_Status status;
  MPI_Request request;
  MPI_Irecv(_nodes_of_sides.begin(), _nodes_of_sides.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, &status);
  MPI_Irecv(_sides_of_cells.begin(), _sides_of_cells.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, &status);
  MPI_Irecv(_cells_of_sides.begin(), _cells_of_sides.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, &status);
}

/*--------------------------------------------------------------------------*/
void SidesAndCells::loadH5(const arma::hdf5_name& spec)
{
  _nodes_of_sides.load(arma::hdf5_name(spec.filename, spec.dsname+"/nodes_of_sides", spec.opts));
  _sides_of_cells.load(arma::hdf5_name(spec.filename, spec.dsname+"/sides_of_cells", spec.opts));
  _cells_of_sides.load(arma::hdf5_name(spec.filename, spec.dsname+"/cells_of_sides", spec.opts));
}
void SidesAndCells::saveH5(const arma::hdf5_name& spec) const
{
  _nodes_of_sides.save(arma::hdf5_name(spec.filename, spec.dsname+"/nodes_of_sides", spec.opts));
  _sides_of_cells.save(arma::hdf5_name(spec.filename, spec.dsname+"/sides_of_cells", spec.opts));
  _cells_of_sides.save(arma::hdf5_name(spec.filename, spec.dsname+"/cells_of_sides", spec.opts));
}
