#include  "Mesh/boundarymeshinfo.hpp"
#include  <cassert>
#include  <mpi.h>

using namespace mesh;

/*--------------------------------------------------------------------------*/
BoundaryMeshInfo::~BoundaryMeshInfo() {}
BoundaryMeshInfo::BoundaryMeshInfo(): GeometryObject(){}
BoundaryMeshInfo::BoundaryMeshInfo( const BoundaryMeshInfo& boundarymeshinfo): GeometryObject(boundarymeshinfo)
{
  _nodes_to_parent = boundarymeshinfo._nodes_to_parent;
  _cells_to_parent = boundarymeshinfo._cells_to_parent;
}
BoundaryMeshInfo& BoundaryMeshInfo::operator=( const BoundaryMeshInfo& boundarymeshinfo)
{
  assert(0);
  GeometryObject::operator=(boundarymeshinfo);
  return *this;
}
std::string BoundaryMeshInfo::getClassName() const
{
  return "BoundaryMeshInfo";
}
std::unique_ptr<GeometryObject> BoundaryMeshInfo::clone() const
{
  return std::unique_ptr<mesh::GeometryObject>(new BoundaryMeshInfo(*this));
}

/*--------------------------------------------------------------------------*/
const alat::armaivec&  BoundaryMeshInfo::getNodesToParent() const {return _nodes_to_parent;}
alat::armaivec&  BoundaryMeshInfo::getNodesToParent() {return _nodes_to_parent;}
const alat::armaimat&  BoundaryMeshInfo::getCellsToParent() const {return _cells_to_parent;}
alat::armaimat&  BoundaryMeshInfo::getCellsToParent() {return _cells_to_parent;}

/*--------------------------------------------------------------------------*/
alat::armaivec BoundaryMeshInfo::getSizes() const
{
  alat::armaivec sizes(3);
  sizes[0] = _nodes_to_parent.size();
  sizes[1] = _cells_to_parent.n_rows;
  sizes[2] = _cells_to_parent.n_cols;
  return sizes;
}
void BoundaryMeshInfo::setSizes(alat::armaivec::const_iterator sizes)
{
  _nodes_to_parent.set_size(sizes[0]);
  _cells_to_parent.set_size(sizes[1], sizes[2]);
}
void BoundaryMeshInfo::send(int neighbor, int tag) const
{
  MPI_Request request;
  MPI_Isend(_nodes_to_parent.begin(), _nodes_to_parent.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Isend(_cells_to_parent.begin(), _cells_to_parent.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
}
void BoundaryMeshInfo::recv(int neighbor, int tag)
{
  MPI_Status status;
  MPI_Request request;
  MPI_Irecv(_nodes_to_parent.begin(), _nodes_to_parent.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, &status);
  MPI_Irecv(_cells_to_parent.begin(), _cells_to_parent.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, &status);
}

/*--------------------------------------------------------------------------*/
void BoundaryMeshInfo::loadH5(const arma::hdf5_name& spec)
{
  _nodes_to_parent.load(arma::hdf5_name(spec.filename, spec.dsname+"/nodes_to_parent", spec.opts));
  _cells_to_parent.load(arma::hdf5_name(spec.filename, spec.dsname+"/cells_to_parent", spec.opts));
}
void BoundaryMeshInfo::saveH5(const arma::hdf5_name& spec) const
{
  _nodes_to_parent.save(arma::hdf5_name(spec.filename, spec.dsname+"/nodes_to_parent", spec.opts));
  _cells_to_parent.save(arma::hdf5_name(spec.filename, spec.dsname+"/cells_to_parent", spec.opts));
}
