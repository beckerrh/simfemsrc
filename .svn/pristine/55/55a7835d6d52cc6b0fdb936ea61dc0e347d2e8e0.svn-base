#include  "Alat/vectoronevariable.hpp"
#include  "Solvers/fem.hpp"
#include  "Solvers/integrationformulapoint.hpp"
#include  "Solvers/integrationformulaline.hpp"
#include  "Solvers/integrationformulatetrahedral.hpp"
#include  "Solvers/integrationformulatriangle.hpp"
#include  <cassert>
#include  <cmath>

using namespace solvers;

/*--------------------------------------------------------------------------*/
Fem::~Fem() {}
Fem::Fem() : solvers::FemInterface() {}
Fem::Fem( const Fem& fem): solvers::FemInterface(fem)
{
  (*this).operator=(fem);
}
Fem& Fem::operator=( const Fem& fem)
{
  solvers::FemInterface::operator=(fem);
  return *this;
}
/*--------------------------------------------------------------------------*/
std::unique_ptr<solvers::IntegrationFormulaInterface> Fem::newFormulaRhs()
{
  return newFormula();
}
std::unique_ptr<solvers::IntegrationFormulaInterface> Fem::newFormula()
{
  if(_mesh->getDimension()==1)
  {
    return std::unique_ptr<solvers::IntegrationFormulaInterface>(new LineTrapez());
    return std::unique_ptr<solvers::IntegrationFormulaInterface>(new LineGauss2());
  }
  else if(_mesh->getDimension()==2)
  {
    return std::unique_ptr<solvers::IntegrationFormulaInterface>(new TriangleRotatedTrapez());
  }
  else
  {
    return std::unique_ptr<solvers::IntegrationFormulaInterface>(new TetrahedralFormula5());
  }
}
std::unique_ptr<solvers::IntegrationFormulaInterface> Fem::newFormulaBdry()
{
  if(_mesh->getDimension()==1)
  {
    return std::unique_ptr<solvers::IntegrationFormulaInterface>(new IntegrationFormulaPoint());
  }
  else if(_mesh->getDimension()==2)
  {
    return std::unique_ptr<solvers::IntegrationFormulaInterface>(new LineGauss2());
  }
  else{return std::unique_ptr<solvers::IntegrationFormulaInterface>(new TriangleRotatedTrapez());}
}
std::unique_ptr<solvers::IntegrationFormulaInterface> Fem::newFormulaErrors()
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
void Fem::initFem(int ivar, const mesh::MeshUnitInterface* mesh, const MeshInfo* meshinfo, int ncomp)
{
  FemInterface::initFem(ivar, mesh, meshinfo, ncomp);
  _formula = newFormula();
  _formulaerrors = newFormulaErrors();
  _formulabdry = newFormulaBdry();
  _formularhs = newFormulaRhs();
  int nlocal=getNPerCell();
  _trafob.set_size(3);
  _trafoA.set_size(3,nlocal-1);
  _femdata.phi.set_size(nlocal);
  _femdata.dphi.set_size(3,nlocal);
  _femdata.u.set_size(ncomp);
  _femdata.ugrad.set_size(3,ncomp);
  _femdata.normal.set_size(3);
  _femdata.uB.set_size(ncomp);
  _femdata.uI.set_size(ncomp);
  _femdata.uBgrad.set_size(3,ncomp);
  _femdata.uIgrad.set_size(3,ncomp);
  _femdata.isI.set_size(nlocal);
  _femdata.ncomp = ncomp;
  _femdata.nlocal = nlocal;
  _femdata.mass.set_size(nlocal,nlocal);
  _femdata.laplace.set_size(nlocal,nlocal);
  _femdata.mass_lumped.set_size(nlocal);
  initData();
}

std::string Fem::getClassName() const
{
  return "Fem";
}

/*--------------------------------------------------------------------------*/
const solvers::IntegrationFormulaInterface* Fem::getFormula() const {return _formula.get();}
const solvers::IntegrationFormulaInterface* Fem::getFormulaRhs() const {return _formularhs.get();}
const solvers::IntegrationFormulaInterface* Fem::getFormulaErrors() const {return _formulaerrors.get();}
const solvers::IntegrationFormulaInterface* Fem::getFormulaBdry() const {return _formulabdry.get();}
const solvers::FemData& Fem::getFemdata() const {return _femdata;}
void Fem::computeGrad(arma::mat& ugrad, const arma::mat& uloc) const
{
  ugrad.fill(arma::fill::zeros);
  for(int ii=0;ii<_femdata.phi.size();ii++)
  {
    for(int icomp=0;icomp<uloc.n_rows;icomp++)
    {
      ugrad.col(icomp) += _femdata.dphi.col(ii)*uloc(icomp,ii);
    }
  }
}
void Fem::computeFunction(arma::vec& u, const arma::mat& uloc) const
{
  u.fill(arma::fill::zeros);
  for(int ii=0;ii<_femdata.phi.size();ii++)
  {
    for(int icomp=0;icomp<uloc.n_rows;icomp++)
    {
      u[icomp] += _femdata.phi[ii]*uloc(icomp,ii);
    }
  }
}

void Fem::_computeData(const arma::mat& uloc)
{
  _femdata.u.fill(arma::fill::zeros);
  _femdata.ugrad.fill(arma::fill::zeros);
  for(int ii=0;ii<_femdata.phi.size();ii++)
  {
    for(int icomp=0;icomp<uloc.n_rows;icomp++)
    {
      _femdata.u[icomp] += _femdata.phi[ii]*uloc(icomp,ii);
      _femdata.ugrad.col(icomp) += _femdata.dphi.col(ii)*uloc(icomp,ii);
    }
  }
}

const FemData& Fem::referencePointWithData(const alat::Node& vhat, double weight, const arma::mat& uloc)
{
  referencePoint(vhat, weight);
  _computeData(uloc);
  return _femdata;
}
const FemData& Fem::referencePointBdryWithData(const alat::Node& vhat, double weight, const arma::mat& uloc)
{
  referencePointBdry(vhat, weight);
  _computeData(uloc);
  return _femdata;
}
const FemData& Fem::referencePointBdryCellWithData(const alat::Node& vhat, double weight, const arma::mat& uloc)
{
  referencePoint(vhat, weight);
  _computeData(uloc);
  _femdata.uBgrad.fill(arma::fill::zeros);
  _femdata.uIgrad.fill(arma::fill::zeros);
  _femdata.uB.fill(arma::fill::zeros);
  _femdata.uI.fill(arma::fill::zeros);
  for(int ii=0;ii<_femdata.phi.size();ii++)
  {
    if(_femdata.isI[ii])
    {
      for(int icomp=0;icomp<uloc.n_rows;icomp++)
      {
        _femdata.uI[icomp] += _femdata.phi[ii]*uloc(icomp,ii);
        _femdata.uIgrad.col(icomp) += _femdata.dphi.col(ii)*uloc(icomp,ii);
      }
    }
    else
    {
      for(int icomp=0;icomp<uloc.n_rows;icomp++)
      {
        _femdata.uB[icomp] += _femdata.phi[ii]*uloc(icomp,ii);
        _femdata.uBgrad.col(icomp) += _femdata.dphi.col(ii)*uloc(icomp,ii);
      }
    }
  }
  return _femdata;
}


/*--------------------------------------------------------------------------*/
void Fem::setVectorIndices(int iK, alat::armaimat& vec_i)const
{
  int nglob = getN();
  alat::armaivec indices(vec_i.n_cols);
  indicesOfCell(iK, indices);
  vec_i.set_size(_ncomp, indices.size());
  for(int ii=0; ii<indices.size();ii++)
  {
    int i = indices[ii];
    for(int icomp=0;icomp<_ncomp;icomp++)
    {
      vec_i(icomp,ii) = icomp*nglob +i;
    }
  }
  // std::cerr << "Fem::setVectorIndices() " << indices.t();
}
/*--------------------------------------------------------------------------*/
void Fem::toP1(alat::VectorOneVariableInterface* uc1, const alat::VectorOneVariableInterface* u)
{
  alat::VectorOneVariable* uc1v = dynamic_cast<alat::VectorOneVariable*>(uc1); assert(uc1v);
  const alat::VectorOneVariable* uv = dynamic_cast<const alat::VectorOneVariable*>(u); assert(uv);
  int nglob = getN();
  int nloccell = getNPerCell();
  int nnodes = _meshinfo->nnodes;
  int nnodespercell = _meshinfo->nnodespercell;
  int ncells = _meshinfo->ncells;
  // std::cerr << "nglob="<<nglob<<" uv->size()="<<uv->size()<<"\n";
  // std::cerr << "nnodes="<<nnodes<<" _meshinfo->nedges="<<_meshinfo->nedges<<" _meshinfo->nsides="<<_meshinfo->nsides<<"\n";
  assert(uv->size()==_ncomp*nglob);
  assert(uv->n()==nglob);
  assert(uv->ncomp()==_ncomp);

  assert(uc1v->size()==_ncomp*nnodes);
  assert(uc1v->n()==nnodes);
  assert(uc1v->ncomp()==_ncomp);

  std::unique_ptr<const solvers::IntegrationFormulaInterface> IF;
  if(_mesh->getDimension()==1) {IF = std::unique_ptr<const solvers::IntegrationFormulaInterface>(new LineTrapez());}
  else if(_mesh->getDimension()==2) {IF = std::unique_ptr<const solvers::IntegrationFormulaInterface>(new TriangleTrapez());}
  else {IF = std::unique_ptr<const solvers::IntegrationFormulaInterface>(new TetrahedronTrapez());}
  alat::armaivec count(nnodes, arma::fill::zeros);
  arma::mat uloc(_ncomp, nloccell);
  alat::armaivec indices(nnodespercell);
  alat::armaimat vec_i(_ncomp, nloccell);
  uc1v->fill(arma::fill::zeros);
  assert(IF->n()==nnodespercell);

  // std::cerr << "uv="<< *uv << "\n";

  for(int iK=0; iK<_meshinfo->ncells;iK++)
  {
    setCell(iK);
    setVectorIndices(iK, vec_i);
    uv->extract(vec_i, uloc);
    for(int k = 0; k < IF->n(); k++)
    {
      const FemData& fem = referencePointWithData(IF->point(k),IF->weight(k), uloc);
      int iN = _meshinfo->nodes_of_cells(k,iK);
      assert(fabs(fem.x-_meshinfo->nodes(0,iN)) < 1e-10);
      assert(fabs(fem.y-_meshinfo->nodes(1,iN)) < 1e-10);
      assert(fabs(fem.z-_meshinfo->nodes(2,iN)) < 1e-10);
      count[iN]++;
      for(int icomp = 0; icomp < _ncomp; icomp++)
      {
        int i = icomp*nnodes +iN;
        (*uc1v)[i] += fem.u(icomp);
      }
    }
  }
  // std::cerr << "count="<< count.t();
  uc1v->scaleIntVector(count);

  // std::cerr << "uc1v="<< *uc1v << "\n";
}
/*--------------------------------------------------------------------------*/
void Fem::fromP1(alat::VectorOneVariableInterface* u, const alat::VectorOneVariableInterface* uc1)
{
  assert(0);
}

/*--------------------------------------------------------------------------*/
void Fem::interpolate(alat::VectorOneVariableInterface* u, const solvers::FunctionInterface& function)
{
  alat::VectorOneVariable* uv = dynamic_cast<alat::VectorOneVariable*>(u); assert(uv);
  int nglob = getN();
  int nloccell = getNPerCell();
  assert(uv->size()==_ncomp*nglob);
  assert(uv->n()==nglob);
  assert(uv->ncomp()==_ncomp);
  arma::mat Mloc(nloccell,nloccell);
  arma::vec uinterpol(_ncomp);
  arma::mat uloc(_ncomp, nloccell), floc(nloccell,_ncomp);
  // massMatrix(0, Mloc, true);
  // Mloc = arma::inv(Mloc);
  const solvers::IntegrationFormulaInterface* IF = getFormula();
  alat::armaivec count(nglob, arma::fill::zeros);
  alat::armaivec indices(nloccell);
  alat::armaimat vec_i(_ncomp, nloccell);
  uv->fill(arma::fill::zeros);
  for(int iK=0; iK<_meshinfo->ncells;iK++)
  {
    setCell(iK);
    indicesOfCell(iK, indices);
    setVectorIndices(iK, vec_i);
    floc.fill(arma::fill::zeros);
    Mloc.fill(arma::fill::zeros);
    for(int k = 0; k < IF->n(); k++)
    {
      const FemData& fem = referencePoint(IF->point(k),IF->weight(k));
      function(uinterpol, fem.x, fem.y, fem.z, _meshinfo->dim);
      // std::cerr << "uinterpol="<<uinterpol.t()<<" " << fem.weight<<"\n";
      for(int ii = 0; ii < nloccell; ii++)
      {
        for(int jj = 0; jj < nloccell; jj++)
        {
          Mloc(ii, jj) += fem.weight*fem.phi[jj]*fem.phi[ii];
        }
        for(int icomp = 0; icomp < _ncomp; icomp++)
        {
          floc(ii,icomp) += fem.weight*uinterpol[icomp]*fem.phi[ii];
        }
      }
    }
    // std::cerr << "Mloc="<<Mloc<<"\n";
    Mloc = arma::inv(Mloc);
    uloc = arma::trans(Mloc*floc);
    // std::cerr << "MlocInv="<<Mloc<<"\n";
    // std::cerr << "floc="<<floc.t()<<"\n";
    // std::cerr << "uloc="<<uloc.t()<<"\n";
    for(int ii = 0; ii < nloccell; ii++)
    {
      count[indices[ii]]++;
    }
    uv->assemble(vec_i, uloc);
  }
  // std::cerr << "uv="<<uv->t()<<"\n";
  // std::cerr << "count="<<count.t()<<"\n";
  uv->scaleIntVector(count);
  // std::cerr << "uv="<<uv->t()<<"\n";
}
/*--------------------------------------------------------------------------*/
void Fem::computeErrors(int iK, solvers::ErrorsMap& errormaps, const arma::mat& uloc, const solvers::FunctionInterface& exactsolutions)
{
  arma::vec uex(_ncomp), uex_x(_ncomp), uex_y(_ncomp), uex_z(_ncomp);
  arma::mat ugradex(3,_ncomp);
  int nloccell = getNPerCell(iK);
  const solvers::IntegrationFormulaInterface* IF = getFormulaErrors();
  setCell(iK);
  for(int k = 0; k < IF->n(); k++)
  {
    const FemData& fem = referencePointWithData(IF->point(k),IF->weight(k), uloc);
    exactsolutions(uex, fem.x, fem.y, fem.z);
    assert(uex.is_finite());
    // std::cerr << "Fem::computeErrors() fem.x="<<fem.x<<" fem.u="<<fem.u<<" uex="<<uex;
    uex -= fem.u;
    if(errormaps.hasKey("L1"))
    {
      errormaps["L1"] += fem.weight*arma::norm(uex, 1);
    }
    if(errormaps.hasKey("Linf"))
    {
      for(int icomp=0;icomp<errormaps["Linf"].size();icomp++)
      {
        errormaps["Linf"][icomp] = fmax(errormaps["Linf"][icomp], fabs(uex[icomp]));
      }
    }
    if(errormaps.hasKey("L2"))
    {
      // std::cerr << fem.weight <<  " : " << uex*uex << "\n";
      errormaps["L2"] += fem.weight*(uex*uex);
    }
    if(errormaps.hasKey("H1"))
    {
      exactsolutions.x(uex_x, fem.x, fem.y, fem.z);
      exactsolutions.y(uex_y, fem.x, fem.y, fem.z);
      exactsolutions.z(uex_z, fem.x, fem.y, fem.z);
      ugradex.row(0) = uex_x;
      ugradex.row(1) = uex_y;
      ugradex.row(2) = uex_z;
      ugradex -= fem.ugrad;
      errormaps["H1"] += fem.weight*arma::dot(ugradex,ugradex);
    }
  }
}
