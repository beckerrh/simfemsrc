#include  "Alat/matrixonevariableinterface.hpp"
#include  "Solvers/pdepartwithintegration.hpp"
#include  "Solvers/meshunitwithdatainterface.hpp"
#include  "Solvers/variable.hpp"
#include  <cassert>

using namespace solvers;

/*--------------------------------------------------------------------------*/
PdePartWithIntegration::~PdePartWithIntegration() {}
PdePartWithIntegration::PdePartWithIntegration(alat::StringList vars): PdePartInterface(vars){}
PdePartWithIntegration::PdePartWithIntegration( const PdePartWithIntegration& pdepartwithfem): PdePartInterface(pdepartwithfem)
{
  assert(0);
}
PdePartWithIntegration& PdePartWithIntegration::operator=( const PdePartWithIntegration& pdepartwithfem)
{
  assert(0);
  PdePartInterface::operator=(pdepartwithfem);
  return *this;
}
std::string PdePartWithIntegration::getClassName() const
{
  return "PdePartWithIntegration";
}
/*--------------------------------------------------------------------------*/
void PdePartWithIntegration::rhsCell(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems)const
{
  for(int ii=0;ii<_ivars.size();ii++)
  {
    int ivar = _ivars[ii];
    const solvers::FemData& fem = *fems[ivar];
    int ncomp = fem.ncomp;
    arma::vec f(ncomp,arma::fill::zeros);
    // std::cerr << "PdePartWithIntegration::rhsCell() " << _application->getClassName() <<"\n";
    if(not _application)
    {
      _error_string("rhsCell","no application set");
    }
    _application->getRightHandSide(ivar)(f, fem.x, fem.y, fem.z);
    int nlocal = fem.nlocal;
    for(int ii=0; ii<nlocal;ii++)
    {
      for(int icomp=0;icomp<ncomp;icomp++)
      {
        floc[ivar](icomp,ii) += f[icomp]*fem.weight*fem.phi[ii];
      }
    }
  }
}
/*--------------------------------------------------------------------------*/
void PdePartWithIntegration::computeRhsCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const
{
  const solvers::IntegrationFormulaInterface* IF = _femforintegration->getFormulaRhs();
  if(_cellisbdry[iK])
  {
    for(int ii=0;ii<_ivars.size();ii++)
    {
      (*_fems)[_ivars[ii]]->setIsi(iK);
    }
    prepareRhsCellBdry(iK);
    for(int k = 0; k < IF->n(); k++)
    {
      for(int ii=0;ii<_ivars.size();ii++)
      {
        int ivar = _ivars[ii];
        _femdatas[ivar] = &(*_fems)[ivar]->referencePointBdryCellWithData(IF->point(k),IF->weight(k), uloc[ivar]);
      }
      rhsBdryCell(floc, _femdatas);
    }
  }
  else
  {
    for(int k = 0; k < IF->n(); k++)
    {
      for(int ii=0;ii<_ivars.size();ii++)
      {
        int ivar = _ivars[ii];
        _femdatas[ivar] = &(*_fems)[ivar]->referencePointWithData(IF->point(k),IF->weight(k), uloc[ivar]);
      }
      // std::cerr << "PdePartWithIntegration::computeRhsCell() AVATN prepareRhsCellBdry\n";
      rhsCell(floc, _femdatas);
    }
  }
}
/*--------------------------------------------------------------------------*/
void PdePartWithIntegration::computeResidualCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const
{
  const solvers::IntegrationFormulaInterface* IF = _femforintegration->getFormula();
  if(_cellisbdry[iK])
  {
    for(int ii=0;ii<_ivars.size();ii++)
    {
      (*_fems)[_ivars[ii]]->setIsi(iK);
    }
    // prepareRhsCellBdry(iK);
    for(int k = 0; k < IF->n(); k++)
    {
      for(int ii=0;ii<_ivars.size();ii++)
      {
        int ivar = _ivars[ii];
        _femdatas[ivar] = &(*_fems)[ivar]->referencePointBdryCellWithData(IF->point(k),IF->weight(k), uloc[ivar]);
      }
      residualBdryCell(floc, _femdatas);
    }
  }
  else
  {
    for(int k = 0; k < IF->n(); k++)
    {
      for(int ii=0;ii<_ivars.size();ii++)
      {
        int ivar = _ivars[ii];
        _femdatas[ivar] = &(*_fems)[ivar]->referencePointWithData(IF->point(k),IF->weight(k), uloc[ivar]);
      }
      residualCell(floc, _femdatas);
    }
  }
}
/*--------------------------------------------------------------------------*/
void PdePartWithIntegration::computeMatrixCell(int iK, solvers::PdePartData::mat& mat, solvers::PdePartData::imat& mat_i, solvers::PdePartData::imat& mat_j, const solvers::PdePartData::vec& uloc)const
{
  const solvers::IntegrationFormulaInterface* IF = _femforintegration->getFormula();
  if(_cellisbdry[iK])
  {
    for(int ii=0;ii<_ivars.size();ii++)
    {
      (*_fems)[_ivars[ii]]->setIsi(iK);
    }
    for(int k = 0; k < IF->n(); k++)
    {
      for(int ii=0;ii<_ivars.size();ii++)
      {
        int ivar = _ivars[ii];
        _femdatas[ivar] = &(*_fems)[ivar]->referencePointBdryCellWithData(IF->point(k),IF->weight(k), uloc[ivar]);
      }
      matrixBdryCell(mat, _femdatas);
    }
  }
  else
  {
    for(int k = 0; k < IF->n(); k++)
    {
      for(int ii=0;ii<_ivars.size();ii++)
      {
        int ivar = _ivars[ii];
        _femdatas[ivar] = &(*_fems)[ivar]->referencePointWithData(IF->point(k),IF->weight(k), uloc[ivar]);
      }
      matrixCell(mat, _femdatas);
    }
  }
}
/*--------------------------------------------------------------------------*/
void PdePartWithIntegration::computeRhsBdry(int color, int iK, int iS, int iil, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const
{
  const solvers::IntegrationFormulaInterface* IF = _femforintegration->getFormulaBdry();
  for(int k = 0; k < IF->n(); k++)
  {
    for(int ii=0;ii<_ivars.size();ii++)
    {
      int ivar = _ivars[ii];
      _femdatas[ivar] = &(*_fems)[ivar]->referencePointBdryWithData(IF->point(k),IF->weight(k), uloc[ivar]);
    }
    rhsBdry(floc, _femdatas);
  }
}
/*--------------------------------------------------------------------------*/
void PdePartWithIntegration::computeResidualBdry(int color, int iK, int iS, int iil, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const
{
  const solvers::IntegrationFormulaInterface* IF = _femforintegration->getFormulaBdry();
  for(int k = 0; k < IF->n(); k++)
  {
    for(int ii=0;ii<_ivars.size();ii++)
    {
      int ivar = _ivars[ii];
      _femdatas[ivar] = &(*_fems)[ivar]->referencePointBdryWithData(IF->point(k),IF->weight(k), uloc[ivar]);
    }
    residualBdry(floc, _femdatas);
  }
}
/*--------------------------------------------------------------------------*/
void PdePartWithIntegration::computeMatrixBdry(int color, int iK, int iS, int iil, solvers::PdePartData::mat& mat, solvers::PdePartData::imat& mat_i, solvers::PdePartData::imat& mat_j, const solvers::PdePartData::vec& uloc)const
{
  const solvers::IntegrationFormulaInterface* IF = _femforintegration->getFormulaBdry();
  for(int k = 0; k < IF->n(); k++)
  {
    for(int ii=0;ii<_ivars.size();ii++)
    {
      int ivar = _ivars[ii];
      _femdatas[ivar] = &(*_fems)[ivar]->referencePointBdryWithData(IF->point(k),IF->weight(k), uloc[ivar]);
    }
    matrixBdry(mat, _femdatas);
  }
}
