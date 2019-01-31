#include  "Mesh/measureofcell.hpp"
#include  "Mesh/meshvisitorinterface.hpp"
#include  "Mesh/enums.hpp"
#include  <cassert>
#include  <mpi.h>

using namespace mesh;

/*--------------------------------------------------------------------------*/
MeasureOfCell::~MeasureOfCell(){}
MeasureOfCell::MeasureOfCell() : GeometryObject(){}
MeasureOfCell::MeasureOfCell( const MeasureOfCell& geometryobject) : GeometryObject(geometryobject)
{
  _measureofcell = geometryobject._measureofcell;
  // _centerofcell = geometryobject._centerofcell;
}
MeasureOfCell& MeasureOfCell::operator=( const MeasureOfCell& geometryobject)
{
	InterfaceBase::operator=(geometryobject);
	assert(0);
	return *this;
}
std::string MeasureOfCell::getClassName() const
{
	return "MeasureOfCell";
}
std::unique_ptr<GeometryObject> MeasureOfCell::clone() const
{
  return std::unique_ptr<mesh::GeometryObject>(new MeasureOfCell(*this));
}

/*--------------------------------------------------------------------------*/
const alat::armavec& MeasureOfCell::getMeasureOfCell() const
{
	return _measureofcell;
}

/*--------------------------------------------------------------------------*/
// const arma::mat& MeasureOfCell::getCenterOfCell() const
// {
// 	return _centerofcell;
// }

/*--------------------------------------------------------------------------*/
alat::armaivec MeasureOfCell::getSizes() const
{
  alat::armaivec sizes(1);
  sizes[0] = _measureofcell.size();
  // sizes[1] = _centerofcell.n_cols;
  return sizes;
}
void MeasureOfCell::setSizes(alat::armaivec::const_iterator sizes)
{
  _measureofcell.set_size(sizes[0]);
  // _centerofcell.set_size(3,sizes[1]);
}
void MeasureOfCell::send(int neighbor, int tag) const
{
  MPI_Request request;
  MPI_Isend(_measureofcell.begin(), _measureofcell.size(), MPI_DOUBLE, neighbor, tag, MPI_COMM_WORLD, &request);
  // MPI_Isend(_centerofcell.begin(), _centerofcell.size(), MPI_DOUBLE, neighbor, tag, MPI_COMM_WORLD, &request);
}
void MeasureOfCell::recv(int neighbor, int tag)
{
  MPI_Status status;
  MPI_Request request;
  MPI_Irecv(_measureofcell.begin(), _measureofcell.size(), MPI_DOUBLE, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, &status);
  // MPI_Irecv(_centerofcell.begin(), _centerofcell.size(), MPI_DOUBLE, neighbor, tag, MPI_COMM_WORLD, &request);
  // MPI_Wait(&request, &status);
}

/*--------------------------------------------------------------------------*/
void MeasureOfCell::loadH5(const arma::hdf5_name& spec)
{
  _measureofcell.load(arma::hdf5_name(spec.filename, spec.dsname+"/measureofcell", spec.opts));
  // _centerofcell.load(arma::hdf5_name(spec.filename, spec.dsname+"/centerofcell", spec.opts));
}
void MeasureOfCell::saveH5(const arma::hdf5_name& spec) const
{
  _measureofcell.save(arma::hdf5_name(spec.filename, spec.dsname+"/measureofcell", spec.opts));
  // _centerofcell.save(arma::hdf5_name(spec.filename, spec.dsname+"/centerofcell", spec.opts));
}

void MeasureOfCell::construct(const mesh::MeshUnitInterface* mesh, const GeometryConstructorInterface* geometryconstructor)
{
  assert(mesh);
  // std::cerr << "MeasureOfCell::construct() mesh=" << mesh->getClassName() << "\n";
  // std::cerr << "MeasureOfCell::construct() mesh=" << mesh->getInfo() << "\n";
  int ncells = mesh->getNCells();
  int nnodespercell = mesh->getNNodesPerCell();
  _measureofcell.set_size(ncells);
  // _centerofcell.set_size(3,ncells);
  // _centerofcell.zeros();
  // double d = 1.0/ (double) nnodespercell;
	for(int iK=0;iK<ncells;iK++)
	{
		_measureofcell[iK] = mesh->getVisitor()->computeMeasureOfCell(iK, mesh);
    // for(int ii=0;ii<nnodespercell;ii++)
    // {
    //   int iN = mesh->getNodeOfCell(ii,iK);
    //   _centerofcell.col(iK) += d*mesh->getNode(iN);
    // }
	}
  // std::cerr << "MeasureOfCell::construct()" << _measureofcell.t();
}
