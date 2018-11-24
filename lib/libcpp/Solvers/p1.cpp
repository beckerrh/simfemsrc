#include  "Alat/matrixonevariable.hpp"
#include  "Alat/matrixallvariables.hpp"
#include  "Alat/vectoronevariable.hpp"
#include  "Mesh/nodesandnodesofcells.hpp"
#include  "Solvers/p1.hpp"
#include  <cassert>

using namespace solvers;

/*--------------------------------------------------------------------------*/
P1::~P1() {}
P1::P1(): solvers::Fem(){}
P1::P1( const P1& P1): solvers::Fem(P1)
{
  (*this).operator=(P1);
}

P1& P1::operator=( const P1& P1)
{
  solvers::Fem::operator=(P1);
  return *this;
}
std::string P1::getClassName() const
{
  return "P1";
}
std::unique_ptr<FemInterface> P1::clone() const
{
  return std::unique_ptr<solvers::FemInterface>(new P1(*this));
}
solverEnums::fem::femtype P1::getType() const {return solverEnums::fem::P1;}

/*--------------------------------------------------------------------------*/
const arma::uvec& P1::getDofIsBdry() const {return _dofisbdry;}
void P1::setCellIsBdry(arma::uvec& cellisbdry)
{
  _dofisbdry.set_size(_meshinfo->nnodes);
  _dofisbdry.fill(arma::fill::zeros);
  for(mesh::MeshUnitInterface::BoundaryInformationMap::const_iterator p=_meshinfo->bdrymesheunitsmap.begin();p!=_meshinfo->bdrymesheunitsmap.end();p++)
  {
    int color = p->first;
    const alat::armaimat& cells_on_bdry = p->second.getCellsOnBdryOfPlain();
    for(int i = 0; i < cells_on_bdry.n_cols; i++)
    {
      int iK = cells_on_bdry(0,i);
      int iS = cells_on_bdry(1,i);
      int iil = cells_on_bdry(2,i);
      for(int ii=0;ii<_meshinfo->nnodesperside;ii++)
      {
        _dofisbdry[_meshinfo->nodes_of_sides(ii,iS)] = color;
      }
    }
  }
  for(int i=0;i<_meshinfo->ncells;i++)
  {
    for(int ii=0;ii<_meshinfo->nnodespercell;ii++)
    {
      if(_dofisbdry[_meshinfo->nodes_of_cells(ii,i)]) {cellisbdry[i]=true;}
    }
  }
}

/*--------------------------------------------------------------------------*/
int P1::getN() const
{
  return _mesh->getNNodes();
}
int P1::getNPerCell(int iK) const
{
  return _mesh->getNNodesPerCell();
}
void P1::indicesOfCell(int iK, alat::armaivec& indices) const
{
  int n = _mesh->getNNodesPerCell();
  // assert(indices.size()==n);
  indices.set_size(n);
  for(int ii=0;ii<n;ii++)
  {
    indices[ii] = _meshinfo->nodes_of_cells(ii, iK);
  }
}
/*--------------------------------------------------------------------------*/
void P1::computeMatrices(int iK)
{
  arma::mat& mass = _femdata.mass;
  arma::mat& laplace = _femdata.laplace;
  alat::armavec& mass_lumped = _femdata.mass_lumped;
  double moc=_meshinfo->measure_of_cells[iK];
  double d = _meshinfo->dim;
  double scalediff = 1.0/(d*d*moc);
  double scalemass = moc/(d+1.0)/(d+2.0);
  double scalemass_lumped = moc/(d+1.0);
  for(int ii=0; ii<_meshinfo->nnodespercell;ii++)
  {
    int iS = _meshinfo->sides_of_cells(ii,iK);
    for(int jj=0; jj<_meshinfo->nnodespercell;jj++)
    {
      int jS = _meshinfo->sides_of_cells(jj,iK);
      double dot = arma::dot(_meshinfo->normals.col(iS), _meshinfo->normals.col(jS));
      laplace(ii,jj) = dot*scalediff*_meshinfo->sigma(ii,iK)*_meshinfo->sigma(jj,iK);
      mass(ii,jj) = scalemass;
    }
    mass(ii,ii) *= 2.0;
    mass_lumped[ii] = scalemass_lumped;
  }
}

/*--------------------------------------------------------------------------*/
void P1::setCell(int iK)
{
  _femdata.J = _meshinfo->measure_of_cells[iK];
  double scale = 1.0/_femdata.J/_meshinfo->dim;
  // std::cerr << "P1::setCell() _femdata.J="<< _femdata.J << " scale=" << scale<<"\n";
  for(int ii=0;ii<_meshinfo->nnodespercell;ii++)
  {
    double d = -_meshinfo->sigma(ii,iK)*scale;
    int iS = _meshinfo->sides_of_cells(ii,iK);
    _femdata.dphi.col(ii) = _meshinfo->normals.col(iS)*d;
    // std::cerr << "P1::setCell() "<< iK << " " << ii <<" normal=" << _meshinfo->normals.col(iS)<<"\n";
  }
  _trafob = _meshinfo->nodes.col(_meshinfo->nodes_of_cells(0,iK));
  for(int ii=1;ii<_meshinfo->nnodespercell;ii++)
  {
    int iN = _meshinfo->nodes_of_cells(ii,iK);
    _trafoA.col(ii-1) = _meshinfo->nodes.col(iN)-_trafob;
  }
  computeMatrices(iK);
}
void P1::setCellBdry(int iK, int iS, int iil)
{
  setCell(iK);
  assert(_meshinfo->sigma(iil,iK)==1.0);
  _femdata.G = arma::norm(_meshinfo->normals.col(iS));
  _femdata.normal = _meshinfo->normals.col(iS)/_femdata.G;
  // _femdata.G /= _meshinfo->dim;
  _iil = iil;
  _femdata.iil=iil;
}
void P1::setIsi(int iK)
{
  // setCell(iK);
  _femdata.isI.fill(arma::fill::ones);
  for(int ii=0;ii<_meshinfo->nnodespercell;ii++)
  {
    if(_dofisbdry[_meshinfo->nodes_of_cells(ii,iK)]) {_femdata.isI[ii] = false;}
  }
}

const FemData& P1::referencePoint(const alat::Node& vhat, double weight)
{
  _femdata.x=_trafob[0];
  _femdata.y=_trafob[1];
  _femdata.z=_trafob[2];
  for(int ii=0;ii<_meshinfo->nnodespercell-1;ii++)
  {
    _femdata.x += _trafoA(0,ii)*vhat[ii];
    _femdata.y += _trafoA(1,ii)*vhat[ii];
    _femdata.z += _trafoA(2,ii)*vhat[ii];
  }
  _femdata.phi[0]=1.0;
  for(int ii=0;ii<_meshinfo->nnodespercell-1;ii++)
  {
    _femdata.phi[0] -= vhat[ii];
    _femdata.phi[ii+1] = vhat[ii];
  }
  // std::cerr << "_femdata.x="<<_femdata.x<< " _femdata.y=" << _femdata.y<< " _femdata.z=" << _femdata.z<<"\n";
  // std::cerr << "vhat="<<vhat<< " _femdata.phi=" << _femdata.phi<<"\n";
  // std::cerr << " _femdata.J="<< _femdata.J<< " weight=" << weight<<"\n";
  _femdata.weight = _femdata.J*weight;
  return _femdata;
}
const FemData& P1::referencePointBdry(const alat::Node& vhat, double weight)
{
  // ????
  if(_mesh->getDimension()==1)
  {
    if(_iil==0) {_vhat.x()=1.0; _vhat.y()=0.0; _vhat.z()=0.0;}
    else if(_iil==1) {_vhat.x()=0.0; _vhat.y()=1.0; _vhat.z()=0.0;}
    else {assert(0);}
  }
  if(_mesh->getDimension()==2)
  {
    if(_iil==0) {_vhat.x()=1.0-vhat.x(); _vhat.y()=vhat.x(); _vhat.z()=0.0;}
    else if(_iil==1) {_vhat.x()=0.0; _vhat.y()=1.0-vhat.x(); _vhat.z()=0.0;}
    else if(_iil==2) {_vhat.x()=vhat.x(); _vhat.y()=0.0; _vhat.z()=0.0;}
    else {assert(0);}
  }
  if(_mesh->getDimension()==3)
  {
    if(_iil==0) {_vhat.x()=1.0-vhat.x()-vhat.y(); _vhat.y()=vhat.x(); _vhat.z()=vhat.y();}
    else if(_iil==1) {_vhat.x()=0.0; _vhat.y()=vhat.x(); _vhat.z()=vhat.y();}
    else if(_iil==2) {_vhat.x()=vhat.x(); _vhat.y()=0.0; _vhat.z()=vhat.y();}
    else if(_iil==3) {_vhat.x()=vhat.x(); _vhat.y()=vhat.y(); _vhat.z()=0.0;}
    else {assert(0);}
  }
  // ????
  referencePoint(_vhat, weight);
  _femdata.weight = _femdata.G*weight;
  return _femdata;
}
/*--------------------------------------------------------------------------*/
void P1::strongDirichlet(int ivar, alat::MatrixAllVariables& A, const alat::IntSet& dircolors)const
{
  for(alat::IntSet::const_iterator p= dircolors.begin(); p!=dircolors.end();p++)
  {
    int color = *p;
    const alat::armaimat& cells_on_bdry = _meshinfo->bdrymesheunitsmap[color].getCellsOnBdryOfPlain();
    for(int i = 0; i < cells_on_bdry.n_cols; i++)
    {
      int iK = cells_on_bdry(0,i);
      int iS = cells_on_bdry(1,i);
      int iil = cells_on_bdry(2,i);
      for(int ii=0;ii<_meshinfo->nnodesperside;ii++)
      {
        int iN = _meshinfo->nodes_of_sides(ii,iS);
        for(int icomp=0;icomp<_ncomp;icomp++)
        {
          int index = icomp*_meshinfo->nnodes + iN;
          A.get(ivar,ivar)->rowIdentity(index);
        }
      }
    }
  }
}
void P1::strongDirichletZero(alat::VectorOneVariableInterface* u, const alat::IntSet& dircolors)const
{
  alat::VectorOneVariable* uv = dynamic_cast<alat::VectorOneVariable*>(u); assert(uv);
  for(alat::IntSet::const_iterator p= dircolors.begin(); p!=dircolors.end();p++)
  {
    int color = *p;
    const alat::armaimat& cells_on_bdry = _meshinfo->bdrymesheunitsmap[color].getCellsOnBdryOfPlain();
    for(int i = 0; i < cells_on_bdry.n_cols; i++)
    {
      int iK = cells_on_bdry(0,i);
      int iS = cells_on_bdry(1,i);
      int iil = cells_on_bdry(2,i);
      for(int ii=0;ii<_meshinfo->nnodesperside;ii++)
      {
        int iN = _meshinfo->nodes_of_sides(ii,iS);
        for(int icomp=0;icomp<_ncomp;icomp++)
        {
          int index = icomp*_meshinfo->nnodes + iN;
          (*uv)[index] = 0.0;
        }
      }
    }
  }
}
void P1::strongDirichlet(alat::VectorOneVariableInterface* u, const solvers::DirichletInterface& dirichlet, const alat::IntSet& dircolors)const
{
  alat::VectorOneVariable* uv = dynamic_cast<alat::VectorOneVariable*>(u); assert(uv);
  // std::cerr << "P1::strongDirichlet() u=" << *uv << "\n";
  alat::armavec udir(_ncomp);
  for(alat::IntSet::const_iterator p= dircolors.begin(); p!=dircolors.end();p++)
  {
    int color = *p;
    const alat::armaimat& cells_on_bdry = _meshinfo->bdrymesheunitsmap[color].getCellsOnBdryOfPlain();
    for(int i = 0; i < cells_on_bdry.n_cols; i++)
    {
      int iK = cells_on_bdry(0,i);
      int iS = cells_on_bdry(1,i);
      int iil = cells_on_bdry(2,i);
      for(int ii=0;ii<_meshinfo->nnodesperside;ii++)
      {
        int iN = _meshinfo->nodes_of_sides(ii,iS);
        // std::cerr << "iN="<<iN << " xyz="<<_meshinfo->nodes(0,iN) << " " << _meshinfo->nodes(1,iN) << " "<< _meshinfo->nodes(2,iN)<<"\n";
        dirichlet(udir, _meshinfo->nodes(0,iN), _meshinfo->nodes(1,iN), _meshinfo->nodes(2,iN));
        // std::cerr << "P1::strongDirichlet() udir=" << udir;
        for(int icomp=0;icomp<_ncomp;icomp++)
        {
          int index = icomp*_meshinfo->nnodes + iN;
          (*uv)[index] = udir[icomp];
        }
      }
    }
  }
}
