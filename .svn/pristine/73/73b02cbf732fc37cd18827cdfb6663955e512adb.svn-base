#include  "Mesh/interfacemeshinfo.hpp"
#include  <cassert>
#include  <mpi.h>

using namespace mesh;

/*--------------------------------------------------------------------------*/
InterfaceMeshInfo::~InterfaceMeshInfo() {}
InterfaceMeshInfo::InterfaceMeshInfo(): GeometryObject(){}
InterfaceMeshInfo::InterfaceMeshInfo( const InterfaceMeshInfo& interfacemeshinfo): GeometryObject(interfacemeshinfo)
{
  _nodes_to_plain = interfacemeshinfo._nodes_to_plain;
  _sides_to_plain = interfacemeshinfo._sides_to_plain;
  _cells_to_plain = interfacemeshinfo._cells_to_plain;
}
InterfaceMeshInfo& InterfaceMeshInfo::operator=( const InterfaceMeshInfo& interfacemeshinfo)
{
  assert(0);
  GeometryObject::operator=(interfacemeshinfo);
  return *this;
}
std::string InterfaceMeshInfo::getClassName() const
{
  return "InterfaceMeshInfo";
}
std::unique_ptr<GeometryObject> InterfaceMeshInfo::clone() const
{
  return std::unique_ptr<mesh::GeometryObject>(new InterfaceMeshInfo(*this));
}

/*--------------------------------------------------------------------------*/
const alat::armaivec&  InterfaceMeshInfo::getNodesToPlain() const {return _nodes_to_plain;}
alat::armaivec&  InterfaceMeshInfo::getNodesToPlain() {return _nodes_to_plain;}
const alat::armaivec&  InterfaceMeshInfo::getSidesToPlain() const {return _sides_to_plain;}
alat::armaivec&  InterfaceMeshInfo::getSidesToPlain() {return _sides_to_plain;}
const alat::armaivec&  InterfaceMeshInfo::getCellsToPlain() const {return _cells_to_plain;}
alat::armaivec&  InterfaceMeshInfo::getCellsToPlain() {return _cells_to_plain;}

/*--------------------------------------------------------------------------*/
alat::armaivec InterfaceMeshInfo::getSizes() const
{
  alat::armaivec sizes(3);
  sizes[0] = _nodes_to_plain.size();
  sizes[1] = _sides_to_plain.size();
  sizes[2] = _cells_to_plain.size();
  return sizes;
}
void InterfaceMeshInfo::setSizes(alat::armaivec::const_iterator sizes)
{
  _nodes_to_plain.set_size(sizes[0]);
  _sides_to_plain.set_size(sizes[1]);
  _cells_to_plain.set_size(sizes[2]);
}
void InterfaceMeshInfo::send(int neighbor, int tag) const
{
  MPI_Request request;
  MPI_Isend(_nodes_to_plain.begin(), _nodes_to_plain.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Isend(_sides_to_plain.begin(), _sides_to_plain.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Isend(_cells_to_plain.begin(), _cells_to_plain.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
}
void InterfaceMeshInfo::recv(int neighbor, int tag)
{
  MPI_Status status;
  MPI_Request request;
  MPI_Irecv(_nodes_to_plain.begin(), _nodes_to_plain.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, &status);
  MPI_Irecv(_sides_to_plain.begin(), _sides_to_plain.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, &status);
  MPI_Irecv(_cells_to_plain.begin(), _cells_to_plain.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, &status);
}

/*--------------------------------------------------------------------------*/
void InterfaceMeshInfo::loadH5(const arma::hdf5_name& spec)
{
  _nodes_to_plain.load(arma::hdf5_name(spec.filename, spec.dsname+"/nodes_to_plain", spec.opts));
  _sides_to_plain.load(arma::hdf5_name(spec.filename, spec.dsname+"/sides_to_plain", spec.opts));
  _cells_to_plain.load(arma::hdf5_name(spec.filename, spec.dsname+"/sides_to_plain", spec.opts));
}
void InterfaceMeshInfo::saveH5(const arma::hdf5_name& spec) const
{
  _nodes_to_plain.save(arma::hdf5_name(spec.filename, spec.dsname+"/nodes_to_plain", spec.opts));
  _sides_to_plain.save(arma::hdf5_name(spec.filename, spec.dsname+"/sides_to_plain", spec.opts));
  _cells_to_plain.save(arma::hdf5_name(spec.filename, spec.dsname+"/cells_to_plain", spec.opts));
}
