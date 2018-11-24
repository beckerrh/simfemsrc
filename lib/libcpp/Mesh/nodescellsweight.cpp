#include  "Mesh/nodescellsweight.hpp"
#include  "Mesh/meshunitinterface.hpp"
#include  "Mesh/nodesandnodesofcells.hpp"
#include  <cassert>
#include  <mpi.h>

using namespace mesh;

/*--------------------------------------------------------------------------*/
NodesCellsWeight::~NodesCellsWeight() {}
NodesCellsWeight::NodesCellsWeight(): GeometryObject(){}
NodesCellsWeight::NodesCellsWeight( const NodesCellsWeight& nodescellsweight): GeometryObject(nodescellsweight)
{
  assert(0);
}
NodesCellsWeight& NodesCellsWeight::operator=( const NodesCellsWeight& nodescellsweight)
{
  assert(0);
  GeometryObject::operator=(nodescellsweight);
  return *this;
}
std::string NodesCellsWeight::getClassName() const
{
  return "NodesCellsWeight";
}
std::unique_ptr<GeometryObject> NodesCellsWeight::clone() const
{
  return std::unique_ptr< GeometryObject>(new NodesCellsWeight(*this));
}

/*--------------------------------------------------------------------------*/
const alat::armavec& NodesCellsWeight::get() const {return _weights;}

alat::armaivec NodesCellsWeight::getSizes() const
{
  alat::armaivec sizes(1);
  sizes[0] = _weights.size();
  return sizes;
}
void NodesCellsWeight::setSizes(alat::armaivec::const_iterator sizes)
{
  _weights.set_size(sizes[0]);
}
void NodesCellsWeight::send(int neighbor, int tag) const
{
  MPI_Request request;
  MPI_Isend(_weights.begin(), _weights.size(), MPI_DOUBLE, neighbor, tag, MPI_COMM_WORLD, &request);
}
void NodesCellsWeight::recv(int neighbor, int tag)
{
  MPI_Status status;
  MPI_Request request;
  MPI_Irecv(_weights.begin(), _weights.size(), MPI_DOUBLE, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, &status);
}
void NodesCellsWeight::loadH5(const arma::hdf5_name& spec)
{
  arma::hdf5_name h5name(spec.filename, spec.dsname+"/nodescellsweight", spec.opts);
  _weights.load(h5name);
}
void NodesCellsWeight::saveH5(const arma::hdf5_name& spec) const
{
  arma::hdf5_name h5name(spec.filename, spec.dsname+"/nodescellsweight", spec.opts);
  _weights.save(h5name);
}
void NodesCellsWeight::construct(const mesh::MeshUnitInterface* mesh, const GeometryConstructorInterface* geometryconstructor)
{
  int ncells = mesh->getNCells();
  int nnodes = mesh->getNNodes();
  int nnodespercell = mesh->getNNodesPerCell();
	_weights.set_size(nnodes);
  _weights.fill(arma::fill::zeros);
  const alat::armaimat&  nodes_of_cells = mesh->getNodesAndNodesOfCells().getNodesOfCells();
	for(int iK=0;iK<ncells;iK++)
	{
    for(int ii=0;ii<nnodespercell;ii++)
		{
			int iN = nodes_of_cells(ii, iK);
      _weights[iN] += 1.0;
    }
	}
  for(int iN=0;iN<nnodes;iN++) { _weights[iN] = 1.0/_weights[iN];}
}
