#include  "Mesh/boundaryinformation.hpp"
#include  <cassert>
#include  <mpi.h>

using namespace mesh;

/*--------------------------------------------------------------------------*/
BoundaryInformation::~BoundaryInformation() {}
BoundaryInformation::BoundaryInformation(): GeometryObject(){}
BoundaryInformation::BoundaryInformation( const BoundaryInformation& boundarymeshinformation): GeometryObject(boundarymeshinformation)
{
  _cells_on_bdry_of_plain = boundarymeshinformation._cells_on_bdry_of_plain;
}
BoundaryInformation& BoundaryInformation::operator=( const BoundaryInformation& boundarymeshinformation)
{
  _cells_on_bdry_of_plain = boundarymeshinformation._cells_on_bdry_of_plain;
  return *this;
}
std::string BoundaryInformation::getClassName() const
{
  return "BoundaryInformation";
}
std::unique_ptr<GeometryObject> BoundaryInformation::clone() const
{
  return std::unique_ptr<mesh::GeometryObject>(new BoundaryInformation(*this));
}

/*--------------------------------------------------------------------------*/
void BoundaryInformation::set_size(int n)
{
  _cells_on_bdry_of_plain.set_size(3,n);
}

/*--------------------------------------------------------------------------*/
const alat::armaimat& BoundaryInformation::getCellsOnBdryOfPlain() const {return _cells_on_bdry_of_plain;}
alat::armaimat& BoundaryInformation::getCellsOnBdryOfPlain() {return _cells_on_bdry_of_plain;}
int BoundaryInformation::size() const {return _cells_on_bdry_of_plain.n_cols;}
int BoundaryInformation::cellid(int i) const {return _cells_on_bdry_of_plain(0,i);}
int BoundaryInformation::sideid(int i) const {return _cells_on_bdry_of_plain(1,i);}
int BoundaryInformation::sidelocalid(int i) const {return _cells_on_bdry_of_plain(2,i);}

/*--------------------------------------------------------------------------*/
alat::armaivec BoundaryInformation::getSizes() const
{
  alat::armaivec sizes(2);
  sizes[0] = _cells_on_bdry_of_plain.n_rows;
  sizes[1] = _cells_on_bdry_of_plain.n_cols;
  // std::cerr << "BoundaryInformation::getSizes = " << sizes.t();
  return sizes;
}
void BoundaryInformation::setSizes(alat::armaivec::const_iterator sizes)
{
  assert(_cells_on_bdry_of_plain.size()==0);
  // std::cerr << "BoundaryInformation::setSizes = " << sizes[0] << " " << sizes[1] << "\n";
  _cells_on_bdry_of_plain.set_size(sizes[0], sizes[1]);
}
void BoundaryInformation::send(int neighbor, int tag) const
{
  MPI_Request request;
  MPI_Isend(_cells_on_bdry_of_plain.begin(), _cells_on_bdry_of_plain.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
}
void BoundaryInformation::recv(int neighbor, int tag)
{
  MPI_Status status;
  MPI_Request request;
  MPI_Irecv(_cells_on_bdry_of_plain.begin(), _cells_on_bdry_of_plain.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, &status);
}

/*--------------------------------------------------------------------------*/
void BoundaryInformation::loadH5(const arma::hdf5_name& spec)
{
  arma::hdf5_name h5name(spec.filename, spec.dsname+"/cells_on_bdry_of_plain", spec.opts);
  _cells_on_bdry_of_plain.load(h5name);
}
void BoundaryInformation::saveH5(const arma::hdf5_name& spec) const
{
  arma::hdf5_name h5name(spec.filename, spec.dsname+"/cells_on_bdry_of_plain", spec.opts);
  _cells_on_bdry_of_plain.save(h5name);
}

/*--------------------------------------------------------------------------*/
std::ostream& mesh::operator<<(std::ostream& os, const BoundaryInformation& boundaryinformation)
{
  os << boundaryinformation.getCellsOnBdryOfPlain();
  return os;
}
