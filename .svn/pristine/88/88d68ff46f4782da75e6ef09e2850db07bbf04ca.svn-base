#include  "Mesh/meshvisitorinterface.hpp"
#include  "Mesh/normals.hpp"
#include  "Mesh/enums.hpp"
#include  "Mesh/sidesandcells.hpp"
#include  <cassert>
#include  <mpi.h>

using namespace mesh;

/*--------------------------------------------------------------------------*/
Normals::~Normals(){}
Normals::Normals() : GeometryObject(){}
Normals::Normals( const Normals& geometryobject) : GeometryObject(geometryobject)
{
	_normals = geometryobject._normals;
	_sigma = geometryobject._sigma;
}
Normals& Normals::operator=( const Normals& geometryobject)
{
	InterfaceBase::operator=(geometryobject);
	assert(0);
	return *this;
}
std::string Normals::getClassName() const
{
	return "Normals";
}
std::unique_ptr<GeometryObject> Normals::clone() const
{
  return std::unique_ptr<mesh::GeometryObject>(new Normals(*this));
}

 /*--------------------------------------------------------------------------*/
const arma::mat& Normals::getNormals() const
{
	return _normals;
}
const arma::fmat& Normals::getSigma() const
{
	return _sigma;
}
/*--------------------------------------------------------------------------*/
alat::armaivec Normals::getSizes() const
{
  alat::armaivec sizes(4);
  sizes[0] = _normals.n_rows;
  sizes[1] = _normals.n_cols;
  sizes[2] = _sigma.n_rows;
  sizes[3] = _sigma.n_cols;
  return sizes;
}
void Normals::setSizes(alat::armaivec::const_iterator sizes)
{
  _normals.set_size(sizes[0], sizes[1]);
  _sigma.set_size(sizes[2], sizes[3]);
}
void Normals::send(int neighbor, int tag) const
{
  MPI_Request request;
  MPI_Isend(_normals.begin(), _normals.size(), MPI_DOUBLE, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Isend(_sigma.begin(), _sigma.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
}
void Normals::recv(int neighbor, int tag)
{
  MPI_Status status;
  MPI_Request request;
  MPI_Irecv(_normals.begin(), _normals.size(), MPI_DOUBLE, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, &status);
  MPI_Irecv(_sigma.begin(), _sigma.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, &status);
	// std::cerr << "Normals::recv() " << _normals.t() << "\n";
}

/*--------------------------------------------------------------------------*/
void Normals::loadH5(const arma::hdf5_name& spec)
{
	_normals.load(arma::hdf5_name(spec.filename, spec.dsname+"/normals", spec.opts));
	_sigma.load(arma::hdf5_name(spec.filename, spec.dsname+"/sigma", spec.opts));
}
void Normals::saveH5(const arma::hdf5_name& spec) const
{
	_normals.save(arma::hdf5_name(spec.filename, spec.dsname+"/normals", spec.opts));
	_sigma.save(arma::hdf5_name(spec.filename, spec.dsname+"/sigma", spec.opts));
}

/*--------------------------------------------------------------------------*/
void Normals::construct(const mesh::MeshUnitInterface* mesh, const GeometryConstructorInterface* geometryconstructor)
{
	const alat::armaimat&  _sides_of_cells = mesh->getSidesAndCells().getSidesOfCells();
	const alat::armaimat&  _cells_of_sides = mesh->getSidesAndCells().getCellsOfSides();

	int nsides = mesh->getNSides();
	int ncells = mesh->getNCells();
	assert(_sides_of_cells.n_cols==ncells);
	int nsidespercell = mesh->getNSidesPerCell();
	_normals.set_size(3, nsides);
	_sigma.set_size(nsidespercell,ncells);
	for(int iK=0;iK<ncells;iK++)
	{
		for(int ii=0;ii<nsidespercell;ii++)
		{
			int iS = _sides_of_cells(ii, iK); //mesh->getSideIdOfCell(iK,ii);
			// if(mesh->getCellIdOfSide(iS,0)==iK)
			if(_cells_of_sides(0,iS)==iK)
			{
				_sigma(ii,iK) = 1.0;
				// int iK1 = mesh->getCellIdOfSide(iS,1);
				int iK1 = _cells_of_sides(1,iS);
				int ii1=-1;
				if (iK1>=0)
				{
					for(int ii2=0; ii2<mesh->getNSidesPerCell(); ii2++)
					{
						// if(mesh->getSideIdOfCell(iK1, ii2)==iS)
						if(_sides_of_cells(ii2, iK1)==iS)
						{
							ii1 = ii2;
							break;
						}
					}
				}
        // mesh->computeNormal(_normals.col(iS), mesh->getNodeOfSide(iS), mesh->getSide(iS), iK, ii, iK1, ii1);
        mesh->getVisitor()->computeNormal(_normals.col(iS), mesh, mesh->getNodeOfSide(iS), mesh->getSide(iS), iK, ii, iK1, ii1);
				if(_debug_level==2) std::cerr << _debug_level << " normal="<<_normals.col(iS).t();
			}
			else
			{
				_sigma(ii,iK) = -1.0;
			}
		}
	}
	// for(int iS=0;iS<nsides;iS++)
	// {
	// 	arma::subview_col<double> normal = _normals.col(iS);
	// 	// arma::vec& normal = _normals.col(iS);
	// 	mesh->computeNormal(normal, iS);
	// 	// std::cerr << "normal = "<< normal << "\n";
	// }
}
