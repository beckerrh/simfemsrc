#include  "Alat/matrixonevariable.hpp"
#include  "Alat/matrixallvariables.hpp"
#include  "Alat/vectoronevariable.hpp"
#include  "Mesh/nodesandnodesofcells.hpp"
#include  "Solvers/cr1.hpp"
#include  <cassert>

using namespace solvers;

/*--------------------------------------------------------------------------*/
CR1::~CR1() {}
CR1::CR1(): solvers::Fem(){}
CR1::CR1( const CR1& c1): solvers::Fem(c1)
{
  (*this).operator=(c1);
}

CR1& CR1::operator=( const CR1& c1)
{
  solvers::Fem::operator=(c1);
  return *this;
}
std::string CR1::getClassName() const
{
  return "CR1";
}
std::unique_ptr<FemInterface> CR1::clone() const
{
  return std::unique_ptr<solvers::FemInterface>(new CR1(*this));
}
solverEnums::fem::femtype CR1::getType() const {return solverEnums::fem::CR1;}

/*--------------------------------------------------------------------------*/
const arma::uvec& CR1::getDofIsBdry() const {return _dofisbdry;}
void CR1::setCellIsBdry(arma::uvec& cellisbdry)
{
  _dofisbdry.set_size(_meshinfo->nsides);
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
      _dofisbdry[iS] = color;
    }
  }
  for(int i=0;i<_meshinfo->ncells;i++)
  {
    for(int ii=0;ii<_meshinfo->nsidespercell;ii++)
    {
      if(_dofisbdry[_meshinfo->sides_of_cells(ii,i)]) {cellisbdry[i]=true;}
    }
  }
}

/*--------------------------------------------------------------------------*/
int CR1::getN() const
{
  return _mesh->getNSides();
}
int CR1::getNPerCell(int iK) const
{
  return _mesh->getNSidesPerCell();
}
void CR1::indicesOfCell(int iK, alat::armaivec& indices) const
{
  int n = _mesh->getNSidesPerCell();
  // assert(indices.size()==n);
  indices.set_size(n);
  for(int ii=0;ii<n;ii++)
  {
    indices[ii] = _meshinfo->sides_of_cells(ii, iK);
  }
}
/*--------------------------------------------------------------------------*/
void CR1::computeMatrices(int iK)
{
  // assert(0);
}

/*--------------------------------------------------------------------------*/
void CR1::setCell(int iK)
{
  _femdata.J = _meshinfo->measure_of_cells[iK];
  double scale = -1.0/_femdata.J;
  // std::cerr << "CR1::setCell() _femdata.J="<< _femdata.J << " scale=" << scale<<"\n";
  for(int ii=0;ii<_meshinfo->nsidespercell;ii++)
  {
    double d = -_meshinfo->sigma(ii,iK)*scale;
    int iS = _meshinfo->sides_of_cells(ii,iK);
    _femdata.dphi.col(ii) = _meshinfo->normals.col(iS)*d;
    // std::cerr << "CR1::setCell() "<< iK << " " << ii <<" normal=" << _meshinfo->normals.col(iS)<<"\n";
  }
  _trafob = _meshinfo->nodes.col(_meshinfo->nodes_of_cells(0,iK));
  for(int ii=1;ii<_meshinfo->nnodespercell;ii++)
  {
    int iN = _meshinfo->nodes_of_cells(ii,iK);
    _trafoA.col(ii-1) = _meshinfo->nodes.col(iN)-_trafob;
  }
  computeMatrices(iK);
}
void CR1::setCellBdry(int iK, int iS, int iil)
{
  setCell(iK);
  assert(_meshinfo->sigma(iil,iK)==1.0);
  _femdata.G = arma::norm(_meshinfo->normals.col(iS));
  _femdata.normal = _meshinfo->normals.col(iS)/_femdata.G;
  // _femdata.G /= _meshinfo->dim;
  _iil = iil;
  _femdata.iil=iil;
}
void CR1::setIsi(int iK)
{
  // setCell(iK);
  _femdata.isI.fill(arma::fill::ones);
  for(int ii=0;ii<_meshinfo->nsidespercell;ii++)
  {
    if(_dofisbdry[_meshinfo->sides_of_cells(ii,iK)]) {_femdata.isI[ii] = false;}
  }
}

const FemData& CR1::referencePoint(const alat::Node& vhat, double weight)
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
  alat::armavec phip1(_meshinfo->nnodespercell);
  phip1[0]=1.0;
  for(int ii=0;ii<_meshinfo->nnodespercell-1;ii++)
  {
    phip1[0] -= vhat[ii];
    phip1[ii+1] = vhat[ii];
  }
  for(int ii=0;ii<_meshinfo->nnodespercell;ii++)
  {
    _femdata.phi[ii] = 1.0 - _meshinfo->dim*phip1[ii];
  }
  // std::cerr << "_femdata.x="<<_femdata.x<< " _femdata.y=" << _femdata.y<< " _femdata.z=" << _femdata.z<<"\n";
  // std::cerr << "vhat="<<vhat<< " _femdata.phi=" << _femdata.phi<<"\n";
  // std::cerr << " _femdata.J="<< _femdata.J<< " weight=" << weight<<"\n";
  _femdata.weight = _femdata.J*weight;
  return _femdata;
}
const FemData& CR1::referencePointBdry(const alat::Node& vhat, double weight)
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
void CR1::strongDirichlet(int ivar, alat::MatrixAllVariables& A, const alat::IntSet& dircolors)const
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
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        int index = icomp*_meshinfo->nsides + iS;
        A.get(ivar,ivar)->rowIdentity(index);
      }
    }
  }
}
void CR1::strongDirichletZero(alat::VectorOneVariableInterface* u, const alat::IntSet& dircolors)const
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
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        int index = icomp*_meshinfo->nsides + iS;
        (*uv)[index] = 0.0;
      }
    }
  }
}
void CR1::strongDirichlet(alat::VectorOneVariableInterface* u, const solvers::DirichletInterface& dirichlet, const alat::IntSet& dircolors)const
{
  alat::VectorOneVariable* uv = dynamic_cast<alat::VectorOneVariable*>(u); assert(uv);
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
      alat::Node xS = _mesh->getNodeOfSide(iS);
      // std::cerr << "iN="<<iN << " xyz="<<nodes(0,iN)<<" "<<nodes(1,iN)<<" "<<nodes(2,iN)<<" "<<nodes.col(iN)<<"\n";
      dirichlet(udir, xS.x(), xS.y(), xS.z());
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        int index = icomp*_meshinfo->nsides + iS;
        (*uv)[index] = udir[icomp];
      }
    }
  }
}
