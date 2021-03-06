#include  "traditionalintegration.hpp"

/*--------------------------------------------------------------------------*/
TraditionalIntegration::~TraditionalIntegration() {}
TraditionalIntegration::TraditionalIntegration(alat::StringList vars): PdePart(vars){}
TraditionalIntegration::TraditionalIntegration( const TraditionalIntegration& pdepartwithfemtraditional): PdePart(pdepartwithfemtraditional)
{
  assert(0);
}
TraditionalIntegration& TraditionalIntegration::operator=( const TraditionalIntegration& pdepartwithfemtraditional)
{
  assert(0);
  PdePart::operator=(pdepartwithfemtraditional);
  return *this;
}
std::string TraditionalIntegration::getClassName() const
{
  return "TraditionalIntegration";
}
/*--------------------------------------------------------------------------*/
const solver_options::pdepart::opts TraditionalIntegration::setOptions()
{
  return solver_options::pdepart::cell;
}
/*--------------------------------------------------------------------------*/
void TraditionalIntegration::rhsCell(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems)const
{
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  _localmodel->beta(_beta,fem.x, fem.y, fem.z, _mesh->getDimension());
  // std::cerr << "TraditionalIntegration::rhsCell() _beta="<<_beta;
  double delta = _supg(fem.J, arma::norm(_beta));
  for(int ii=0;ii<_ivars.size();ii++)
  {
    int ivar = _ivars[ii];
    int ncomp = fem.ncomp;
    alat::armavec f(ncomp,arma::fill::zeros);
    _application->getRightHandSide(ivar)(f, fem.x, fem.y, fem.z);
    // std::cerr << "TraditionalIntegration::rhsCell() f="<<f;
    int nlocal = fem.nlocal;
    for(int ii=0; ii<nlocal;ii++)
    {
      double betadphi = delta * arma::dot(fem.dphi.col(ii),_beta);
      for(int icomp=0;icomp<ncomp;icomp++)
      {
        floc[ivar][icomp*_nlocal + ii] += f[icomp]*fem.weight*fem.phi[ii];
        floc[ivar][icomp*_nlocal + ii] += f[icomp]*fem.weight*betadphi;
      }
    }
  }
}
/*--------------------------------------------------------------------------*/
void TraditionalIntegration::residualCell(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems) const
{
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  assert(floc.n_rows==1);
  double d = fem.weight*_localmodel->diffusion(fem.x, fem.y, fem.z);
  alat::armavec Fu(1);
  _localmodel->reaction(Fu, fem.u);
  _localmodel->beta(_beta,fem.x, fem.y, fem.z, _mesh->getDimension());
  double delta = _supg(fem.J, arma::norm(_beta));
  double advection=arma::dot(fem.ugrad.col(0),_beta);
  for(int ii=0;ii<_nlocal;ii++)
  {
    double betadphi = delta * arma::dot(fem.dphi.col(ii),_beta);
    floc[_ivar][ii] += fem.weight* fem.phi(ii)*Fu[0];
    floc[_ivar][ii] += fem.weight* fem.phi(ii)*advection;
    floc[_ivar][ii] += fem.weight* betadphi*Fu[0];
    floc[_ivar][ii] += fem.weight* betadphi*advection;
    floc[_ivar][ii] += d* arma::dot(fem.dphi.col(ii),fem.ugrad.col(0));
  }
}
/*--------------------------------------------------------------------------*/
void TraditionalIntegration::matrixCell(solvers::PdePartData::mat& mat, const solvers::FemDatas& fems)const
{
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  double d = fem.weight*_localmodel->diffusion(fem.x, fem.y, fem.z);
  alat::armamat Fu(_ncomp,_ncomp);
  _localmodel->reaction_d(Fu, fem.u);
  _localmodel->beta(_beta,fem.x, fem.y, fem.z, _mesh->getDimension());
  double delta = _supg(fem.J, arma::norm(_beta));
  for(int ii=0; ii<_nlocal;ii++)
  {
    double betadphi = delta * arma::dot(fem.dphi.col(ii),_beta);
    for(int jj=0; jj<_nlocal;jj++)
    {
      double betadphj = arma::dot(fem.dphi.col(jj),_beta);
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        mat(_ivar,_ivar)(icomp*_nlocal+ii, icomp*_nlocal+jj) += fem.weight* fem.phi(ii)*betadphj;
        mat(_ivar,_ivar)(icomp*_nlocal+ii, icomp*_nlocal+jj) += fem.weight* betadphi*betadphj;
        mat(_ivar,_ivar)(icomp*_nlocal+ii, icomp*_nlocal+jj) += d* arma::dot(fem.dphi.col(ii),fem.dphi.col(jj));
        for(int jcomp=0;jcomp<_ncomp;jcomp++)
        {
          mat(_ivar,_ivar)(icomp*_nlocal+ii, jcomp*_nlocal+jj) += fem.weight* fem.phi(ii)*fem.phi(jj)*Fu(icomp, jcomp);
          mat(_ivar,_ivar)(icomp*_nlocal+ii, jcomp*_nlocal+jj) += fem.weight* betadphi*fem.phi(jj)*Fu(icomp, jcomp);
        }
      }
    }
  }
}
