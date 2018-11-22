#include  "Alat/matrixonevariable.hpp"
#include  "Alat/matrixallvariables.hpp"
#include  "Alat/vectoronevariable.hpp"
#include  "Mesh/nodesandnodesofcells.hpp"
#include  "Solvers/p2.hpp"
#include  <cassert>
#include  "Solvers/integrationformulapoint.hpp"
#include  "Solvers/integrationformulaline.hpp"
#include  "Solvers/integrationformulatetrahedral.hpp"
#include  "Solvers/integrationformulatriangle.hpp"

using namespace solvers;

/*--------------------------------------------------------------------------*/
P2::~P2() {}
P2::P2(): solvers::Fem(){}
P2::P2( const P2& c1): solvers::Fem(c1)
{
  (*this).operator=(c1);
}

P2& P2::operator=( const P2& c1)
{
  solvers::Fem::operator=(c1);
  return *this;
}
std::string P2::getClassName() const
{
  return "P2";
}
std::unique_ptr<FemInterface> P2::clone() const
{
  return std::unique_ptr<solvers::FemInterface>(new P2(*this));
}
solverEnums::fem::femtype P2::getType() const {return solverEnums::fem::P2;}

/*--------------------------------------------------------------------------*/
std::unique_ptr<solvers::IntegrationFormulaInterface> P2::newFormula()
{
  if(_mesh->getDimension()==1)
  {
    return std::unique_ptr<solvers::IntegrationFormulaInterface>(new LineGauss2());
  }
  else if(_mesh->getDimension()==2)
  {
    return std::unique_ptr<solvers::IntegrationFormulaInterface>(new Triangle7());
  }
  else
  {
    return std::unique_ptr<solvers::IntegrationFormulaInterface>(new TetrahedralFormula5());
  }
}
std::unique_ptr<solvers::IntegrationFormulaInterface> P2::newFormulaBdry()
{
  if(_mesh->getDimension()==1)
  {
    return std::unique_ptr<solvers::IntegrationFormulaInterface>(new IntegrationFormulaPoint());
  }
  else if(_mesh->getDimension()==2)
  {
    return std::unique_ptr<solvers::IntegrationFormulaInterface>(new LineGauss2());
  }
  else{return std::unique_ptr<solvers::IntegrationFormulaInterface>(new Triangle12());}
}
std::unique_ptr<solvers::IntegrationFormulaInterface> P2::newFormulaErrors()
{
  if(_mesh->getDimension()==1)
  {
    return std::unique_ptr<solvers::IntegrationFormulaInterface>(new LineGauss4());
  }
  else if(_mesh->getDimension()==2)
  {
    return std::unique_ptr<solvers::IntegrationFormulaInterface>(new TriangleSimpson());
  }
  else
  {
    return std::unique_ptr<solvers::IntegrationFormulaInterface>(new TetrahedralFormula5());
  }
}

/*--------------------------------------------------------------------------*/
void P2::initData()
{
  _phip1.set_size(_dim+1);
  _phip1grad.set_size(3,_dim+1);
}

/*--------------------------------------------------------------------------*/
const arma::uvec& P2::getDofIsBdry() const {return _dofisbdry;}
void P2::setCellIsBdry(arma::uvec& cellisbdry)
{
  _dofisbdry.set_size(getN());
  _dofisbdry.fill(arma::fill::zeros);
  int nnodes = _meshinfo->nnodes;
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
  for(mesh::MeshUnitInterface::BoundaryInformationMap::const_iterator p=_meshinfo->bdrymesheunitsmap.begin();p!=_meshinfo->bdrymesheunitsmap.end();p++)
  {
    int color = p->first;
    const alat::armaimat& cells_on_bdry = p->second.getCellsOnBdryOfPlain();
    for(int i = 0; i < cells_on_bdry.n_cols; i++)
    {
      int iK = cells_on_bdry(0,i);
      int iS = cells_on_bdry(1,i);
      int iil = cells_on_bdry(2,i);
      for(int ii=0;ii<_meshinfo->nedgespercell;ii++)
      {
        int iE = _meshinfo->edges_of_cells(ii,i);
        if(_dofisbdry[_meshinfo->nodes_of_edges(0,iE)] and _dofisbdry[_meshinfo->nodes_of_edges(1,iE)])
        {
          _dofisbdry[nnodes+iE] = color;
        }
      }
    }
  }
  for(int i=0;i<_meshinfo->ncells;i++)
  {
    for(int ii=0;ii<_meshinfo->nnodespercell;ii++)
    {
      if(_dofisbdry[_meshinfo->nodes_of_cells(ii,i)]) {cellisbdry[i]=true;}
    }
    for(int ii=0;ii<_meshinfo->nedgespercell;ii++)
    {
      if(_dofisbdry[nnodes+_meshinfo->edges_of_cells(ii,i)]) {cellisbdry[i]=true;}
    }
  }
}

/*--------------------------------------------------------------------------*/
int P2::getN() const
{
  return _meshinfo->nnodes + _meshinfo->nedges;
}
int P2::getNPerCell(int iK) const
{
  return _meshinfo->nnodespercell + _meshinfo->nedgespercell;
}
void P2::indicesOfCell(int iK, alat::armaivec& indices) const
{
  assert(indices.size()==getNPerCell());
  int n = _mesh->getNNodesPerCell();
  for(int ii=0;ii<n;ii++)
  {
    indices[ii] = _meshinfo->nodes_of_cells(ii, iK);
  }
  for(int ii=0;ii<_mesh->getNEdgesPerCell();ii++)
  {
    indices[n+ii] = _meshinfo->nnodes + _meshinfo->edges_of_cells(ii,iK);
  }
}
/*--------------------------------------------------------------------------*/
void P2::computeMatrices(int iK)
{
  // assert(0);
}

/*--------------------------------------------------------------------------*/
void P2::setCell(int iK)
{
  _iK=iK;
  _femdata.J = _meshinfo->measure_of_cells[iK];
  double scale = 1.0/_femdata.J/_meshinfo->dim;
  // std::cerr << "P2::setCell() _femdata.J="<< _femdata.J << " scale=" << scale<<"\n";
  for(int ii=0;ii<_meshinfo->nnodespercell;ii++)
  {
    double d = -_meshinfo->sigma(ii,iK)*scale;
    int iS = _meshinfo->sides_of_cells(ii,iK);
    _phip1grad.col(ii) = _meshinfo->normals.col(iS)*d;
    // std::cerr << "P2::setCell() "<< iK << " " << ii <<" normal=" << _meshinfo->normals.col(iS)<<"\n";
  }
  _trafob = _meshinfo->nodes.col(_meshinfo->nodes_of_cells(0,iK));
  for(int ii=1;ii<_meshinfo->nnodespercell;ii++)
  {
    int iN = _meshinfo->nodes_of_cells(ii,iK);
    _trafoA.col(ii-1) = _meshinfo->nodes.col(iN)-_trafob;
  }
  computeMatrices(iK);
}
void P2::setCellBdry(int iK, int iS, int iil)
{
  setCell(iK);
  assert(_meshinfo->sigma(iil,iK)==1.0);
  _femdata.G = arma::norm(_meshinfo->normals.col(iS));
  _femdata.normal = _meshinfo->normals.col(iS)/_femdata.G;
  // _femdata.G /= _meshinfo->dim;
  _iil = iil;
  _femdata.iil=iil;
}
void P2::setIsi(int iK)
{
  _femdata.isI.fill(arma::fill::ones);
  for(int ii=0;ii<_meshinfo->nnodespercell;ii++)
  {
    if(_dofisbdry[_meshinfo->nodes_of_cells(ii,iK)]) {_femdata.isI[ii] = false;}
  }
  for(int ii=0;ii<_meshinfo->nedgespercell;ii++)
  {
    if(_dofisbdry[_meshinfo->nnodes+_meshinfo->edges_of_cells(ii,iK)]) {_femdata.isI[ii] = false;}
  }
}

const FemData& P2::referencePoint(const alat::Node& vhat, double weight)
{
  _femdata.x=_trafob[0];
  _femdata.y=_trafob[1];
  _femdata.z=_trafob[2];
  int nnpc = _meshinfo->nnodespercell;
  for(int ii=0;ii<nnpc-1;ii++)
  {
    _femdata.x += _trafoA(0,ii)*vhat[ii];
    _femdata.y += _trafoA(1,ii)*vhat[ii];
    _femdata.z += _trafoA(2,ii)*vhat[ii];
  }
  _phip1[0]=1.0;
  for(int ii=0;ii<nnpc-1;ii++)
  {
    _phip1[0] -= vhat[ii];
    _phip1[ii+1] = vhat[ii];
  }
  //
  for(int ii=0;ii<_meshinfo->nedgespercell;ii++)
  {
    int i0 = _meshinfo->localnodes_of_edges_of_cell(0,ii,_iK);
    int i1 = _meshinfo->localnodes_of_edges_of_cell(1,ii,_iK);
    _femdata.phi[nnpc+ii] = 4.0*_phip1[i0]*_phip1[i1];
    _femdata.dphi.col(nnpc+ii) = 4.0*_phip1grad.col(i0)*_phip1[i1] + 4.0*_phip1[i0]*_phip1grad.col(i1);
  }
  for(int ii=0;ii<nnpc;ii++)
  {
    _femdata.phi[ii] = _phip1[ii];
    _femdata.dphi.col(ii) = _phip1grad.col(ii);
  }
  for(int ii=0;ii<_meshinfo->nedgespercell;ii++)
  {
    int i0 = _meshinfo->localnodes_of_edges_of_cell(0,ii,_iK);
    int i1 = _meshinfo->localnodes_of_edges_of_cell(1,ii,_iK);
    _femdata.phi[i0] -= 0.5* _femdata.phi[nnpc+ii];
    _femdata.phi[i1] -= 0.5* _femdata.phi[nnpc+ii];
    _femdata.dphi.col(i0) -= 0.5 * _femdata.dphi.col(nnpc+ii);
    _femdata.dphi.col(i1) -= 0.5 * _femdata.dphi.col(nnpc+ii);
  }
  //
  _femdata.weight = _femdata.J*weight;
  return _femdata;
}
const FemData& P2::referencePointBdry(const alat::Node& vhat, double weight)
{
  assert(0);
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
void P2::strongDirichlet(int ivar, alat::MatrixAllVariables& A, const alat::IntSet& dircolors)const
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
void P2::strongDirichletZero(alat::VectorOneVariableInterface* u, const alat::IntSet& dircolors)const
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
void P2::strongDirichlet(alat::VectorOneVariableInterface* u, const solvers::DirichletInterface& dirichlet, const alat::IntSet& dircolors)const
{
  alat::VectorOneVariable* uv = dynamic_cast<alat::VectorOneVariable*>(u); assert(uv);
  arma::vec udir(_ncomp);
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
        // std::cerr << "iN="<<iN << " xyz="<<nodes(0,iN)<<" "<<nodes(1,iN)<<" "<<nodes(2,iN)<<" "<<nodes.col(iN)<<"\n";
        dirichlet(udir, _meshinfo->nodes(0,iN), _meshinfo->nodes(1,iN), _meshinfo->nodes(2,iN));
        for(int icomp=0;icomp<_ncomp;icomp++)
        {
          int index = icomp*_meshinfo->nnodes + iN;
          (*uv)[index] = udir[icomp];
        }
      }
    }
  }
}
