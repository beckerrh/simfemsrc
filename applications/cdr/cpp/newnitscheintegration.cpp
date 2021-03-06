#include  "newnitscheintegration.hpp"
#include  "Solvers/variable.hpp"
#include  "Solvers/meshunitwithdatainterface.hpp"
#include  <cassert>
#include  <cmath>

/*--------------------------------------------------------------------------*/
NewNitscheIntegration::~NewNitscheIntegration() {}
NewNitscheIntegration::NewNitscheIntegration(alat::StringList vars): TraditionalIntegration(vars){}
NewNitscheIntegration::NewNitscheIntegration( const NewNitscheIntegration& nitscheintegration): TraditionalIntegration(nitscheintegration)
{
  assert(0);
}
NewNitscheIntegration& NewNitscheIntegration::operator=( const NewNitscheIntegration& nitscheintegration)
{
  assert(0);
  TraditionalIntegration::operator=(nitscheintegration);
  return *this;
}
std::string NewNitscheIntegration::getClassName() const
{
  return "NewNitscheIntegration";
}
/*--------------------------------------------------------------------------*/
const solver_options::pdepart::opts NewNitscheIntegration::setOptions()
{
  // std::cerr << "NewNitscheIntegration::setOptions()\n";
  // assert(0);
  const solver_options::pdepart::opts opts= solver_options::pdepart::cell + solver_options::pdepart::bdry;
  return opts;
}
void NewNitscheIntegration::rhsBdry(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems)const
{
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  assert(floc.n_rows==1);
  _localmodel->beta(_beta,fem.x, fem.y, fem.z, _mesh->getDimension());
  double bn = arma::dot(_beta, fem.normal);
  if(bn<0.0)
  {
    alat::armavec dir(1);
    _application->getDirichlet(_ivar)(dir, fem.x, fem.y, fem.z);
    for(int ii=0;ii<_nlocal;ii++)
    {
      floc[_ivar][ii] -= fem.weight*bn*fem.phi[ii]*dir[0];
    }
  }
}
void NewNitscheIntegration::residualBdry(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems)const
{
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  _localmodel->beta(_beta,fem.x, fem.y, fem.z, _mesh->getDimension());
  double bn = arma::dot(_beta, fem.normal);
  if(bn<0.0)
  {
    for(int ii=0;ii<_nlocal;ii++)
    {
      floc[_ivar][ii] -= fem.weight*bn*fem.phi[ii]*fem.u[0];
    }
  }
}
void NewNitscheIntegration::matrixBdry(solvers::PdePartData::mat& mat, const solvers::FemDatas& fems)const
{
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  _localmodel->beta(_beta,fem.x, fem.y, fem.z, _mesh->getDimension());
  double bn = arma::dot(_beta, fem.normal);
  int count=0;
  for(int ii=0; ii<_nlocal;ii++)
  {
    for(int jj=0; jj<_nlocal;jj++)
    {
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        if(bn<0.0)
        {
          mat(_ivar,_ivar)(icomp*_nlocal+ii, icomp*_nlocal+jj) -= fem.weight*bn*fem.phi[ii]*fem.phi[jj];
        }
        count++;
      }
    }
  }
}

void NewNitscheIntegration::prepareRhsCellBdry(int iK) const
{
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  alat::armavec dir(1);
  _udirloc.set_size(_ncomp*_nlocal);
  _udirgrad.set_size(3, _ncomp);
  _udir.set_size(_ncomp);
  _udirloc.fill(arma::fill::zeros);
  // std::cerr << "femdata.isI="<<femdata.isI.t() << "\n";
  for(int ii=0;ii<fem.isI.size();ii++)
  {
    if(not fem.isI[ii])
    {
      if((*_fems)[_ivar]->getType()==solverEnums::fem::P1)
      {
        int iN = _meshinfo->nodes_of_cells(ii,iK);
        _application->getDirichlet(_ivar)(dir, _meshinfo->nodes(0,iN), _meshinfo->nodes(1,iN), _meshinfo->nodes(2,iN));
        _udirloc[ii] = dir[0];
      }
      else if((*_fems)[_ivar]->getType()==solverEnums::fem::CR1)
      {
        int iS = _meshinfo->sides_of_cells(ii,iK);
        alat::Node xS = _mesh->getNodeOfSide(iS);
        _application->getDirichlet(_ivar)(dir, xS.x(), xS.y(), xS.z());
        _udirloc[ii] = dir[0];
      }
      else
      {
        assert(0);
      }
    }
  }
  // std::cerr << "_udirloc="<<_udirloc << "\n";
}
/*--------------------------------------------------------------------------*/
void NewNitscheIntegration::rhsBdryCell(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems)const
{
  // std::cerr << "NewNitscheIntegration::rhsBdryCell()\n";
  // assert(0);
  rhsCell(floc, fems);

  const solvers::FemData& fem=*fems[_ivar];
  assert(floc.n_rows==1);
  alat::armavec f(_ncomp,arma::fill::zeros);
  (*_fems)[_ivar]->computeGrad(_udirgrad, _udirloc);
  (*_fems)[_ivar]->computeFunction(_udir, _udirloc);
  _application->getRightHandSide(_ivar)(f, fem.x, fem.y, fem.z);
  double d = fem.weight*_localmodel->diffusion(fem.x, fem.y, fem.z);
  for(int ii=0; ii<_nlocal;ii++)
  {
    double dot = d*arma::dot(fem.dphi.col(ii),_udirgrad.col(0));
    if(fem.isI[ii])
    {
      if(_symmetric)
      {
        for(int icomp=0;icomp<_ncomp;icomp++)
        {
          floc[_ivar][icomp*_nlocal + ii] -= dot;
        }
      }
    }
    else
    {
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        if(_symmetric)
        {
          floc[_ivar][icomp*_nlocal + ii] -= f[icomp]*fem.weight*fem.phi[ii];
        }
        floc[_ivar][icomp*_nlocal + ii] += dot;
        floc[_ivar][icomp*_nlocal + ii] += (_gamma-1.0)*dot;
      }
    }
  }
}
/*--------------------------------------------------------------------------*/
void NewNitscheIntegration::residualBdryCell(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems)const
{
  residualCell(floc, fems);

  const solvers::FemData& fem=*fems[_ivar];
  assert(floc.n_rows==1);
  double d = fem.weight*_localmodel->diffusion(fem.x, fem.y, fem.z);
  for(int ii=0;ii<_nlocal;ii++)
  {
    if(fem.isI[ii])
    {
      if(_symmetric)
      {
        floc[_ivar][ii] -= d* arma::dot(fem.dphi.col(ii),fem.uBgrad.col(0));
      }
    }
    else
    {
      floc[_ivar][ii] -= d* arma::dot(fem.dphi.col(ii),fem.uIgrad.col(0));
      floc[_ivar][ii] += (_gamma-1.0)*d* arma::dot(fem.dphi.col(ii),fem.uBgrad.col(0));
    }
  }
}
/*--------------------------------------------------------------------------*/
void NewNitscheIntegration::matrixBdryCell(solvers::PdePartData::mat& mat, const solvers::FemDatas& fems)const
{
  matrixCell(mat, fems);

  const solvers::FemData& fem=*fems[_ivar];
  double d = fem.weight*_localmodel->diffusion(fem.x, fem.y, fem.z);
  int count=0;
  for(int ii=0; ii<_nlocal;ii++)
  {
    for(int jj=0; jj<_nlocal;jj++)
    {
      bool sub;
      if(_symmetric)
      {
        sub = (not fem.isI[ii] and fem.isI[jj]) or (fem.isI[ii] and not fem.isI[jj]);
      }
      else
      {
        sub = (not fem.isI[ii] and fem.isI[jj]);
      }
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        if(sub)
        {
          mat(_ivar,_ivar)(icomp*_nlocal+ii, icomp*_nlocal+jj) -= d* arma::dot(fem.dphi.col(ii),fem.dphi.col(jj));
        }
        if(not fem.isI[ii] and not fem.isI[jj])
        {
          mat(_ivar,_ivar)(icomp*_nlocal+ii, icomp*_nlocal+jj) += (_gamma-1.0)*d* arma::dot(fem.dphi.col(ii),fem.dphi.col(jj));
        }
        count++;
      }
    }
  }
}
