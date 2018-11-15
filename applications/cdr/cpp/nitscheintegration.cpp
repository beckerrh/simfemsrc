#include  "nitscheintegration.hpp"
#include  <cassert>

/*--------------------------------------------------------------------------*/
NitscheIntegration::~NitscheIntegration() {}
NitscheIntegration::NitscheIntegration(alat::StringList vars): TraditionalIntegration(vars){}
NitscheIntegration::NitscheIntegration( const NitscheIntegration& nitscheintegration): TraditionalIntegration(nitscheintegration)
{
  assert(0);
}
NitscheIntegration& NitscheIntegration::operator=( const NitscheIntegration& nitscheintegration)
{
  assert(0);
  TraditionalIntegration::operator=(nitscheintegration);
  return *this;
}
std::string NitscheIntegration::getClassName() const
{
  return "NitscheIntegration";
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
const solver_options::pdepart::opts NitscheIntegration::setOptions()
{
  std::cerr << "NitscheIntegration::setOptions()\n";
  return solver_options::pdepart::cell + solver_options::pdepart::bdry;
}
/*--------------------------------------------------------------------------*/
void NitscheIntegration::rhsBdry(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems)const
{
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  assert(floc.n_rows==1);
  arma::vec dir(1);
  _application->getDirichlet(_ivar)(dir, fem.x, fem.y, fem.z);
  double d = fem.weight*_localmodel->diffusion(fem.x, fem.y, fem.z);
  double gamma = _gamma*fem.G/fem.J;
  _localmodel->beta(_beta,fem.x, fem.y, fem.z, _mesh->getDimension());
  double bn = arma::dot(_beta, fem.normal);
  for(int ii=0;ii<_nlocal;ii++)
  {
    double dphin = arma::dot(fem.dphi.col(ii),fem.normal);
    if(_symmetric)
    {
      floc[_ivar](0,ii) -= d* dphin*dir[0];
    }
    floc[_ivar](0,ii) += gamma*d* fem.phi[ii]*dir[0];
    if(bn<0.0)
    {
      floc[_ivar](0,ii) -= fem.weight*bn*fem.phi[ii]*dir[0];
    }
  }
}
void NitscheIntegration::residualBdry(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems)const
{
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  assert(floc.n_rows==1);
  double d = fem.weight*_localmodel->diffusion(fem.x, fem.y, fem.z);
  double gamma = _gamma*fem.G/fem.J;
  _localmodel->beta(_beta,fem.x, fem.y, fem.z, _mesh->getDimension());
  double bn = arma::dot(_beta, fem.normal);
  double dun = arma::dot(fem.normal,fem.ugrad.col(0));
  for(int ii=0;ii<_nlocal;ii++)
  {
    double dphin = arma::dot(fem.dphi.col(ii),fem.normal);
    floc[_ivar](0,ii) -= d* fem.phi[ii]*dun;
    if(_symmetric)
    {
      floc[_ivar](0,ii) -= d* dphin*fem.u[0];
    }
    floc[_ivar](0,ii) += gamma*d* fem.phi[ii]*fem.u[0];
    if(bn<0.0)
    {
      floc[_ivar](0,ii) -= fem.weight*bn*fem.phi[ii]*fem.u[0];
    }
  }
}
void NitscheIntegration::matrixBdry(solvers::PdePartData::mat& mat, const solvers::FemDatas& fems)const
{
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  double d = fem.weight*_localmodel->diffusion(fem.x, fem.y, fem.z);
  int count=0;
  double gamma = _gamma*fem.G/fem.J;
  _localmodel->beta(_beta,fem.x, fem.y, fem.z, _mesh->getDimension());
  double bn = arma::dot(_beta, fem.normal);
  for(int ii=0; ii<_nlocal;ii++)
  {
    double dphin = arma::dot(fem.dphi.col(ii),fem.normal);
    for(int jj=0; jj<_nlocal;jj++)
    {
      double dphjn = arma::dot(fem.normal,fem.dphi.col(jj));
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        mat(_ivar,_ivar)[count] -= d* fem.phi[ii]*dphjn;
        if(_symmetric)
        {
          mat(_ivar,_ivar)[count] -= d* dphin*fem.phi[jj];
        }
        mat(_ivar,_ivar)[count] += gamma*d* fem.phi[ii]*fem.phi[jj];
        if(bn<0.0)
        {
          mat(_ivar,_ivar)[count] -= fem.weight*bn*fem.phi[ii]*fem.phi[jj];
        }
        count++;
      }
    }
  }
}
