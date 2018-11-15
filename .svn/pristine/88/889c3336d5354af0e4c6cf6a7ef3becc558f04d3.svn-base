#include  "Alat/vectoronevariable.hpp"
#include  "Alat/sparsitypatternsoft.hpp"
#include  "Mesh/measureofcell.hpp"
#include  "Mesh/normals.hpp"
#include  "Mesh/nodesandnodesofcells.hpp"
#include  "Solvers/variable.hpp"
#include  "Solvers/meshunitwithdata.hpp"
#include  "Solvers/p1.hpp"
#include  "newnitsche.hpp"
#include  "model.hpp"
#include  <cassert>

/*--------------------------------------------------------------------------*/
NewNitsche::~NewNitsche() {}
NewNitsche::NewNitsche(alat::StringList vars): Traditional(vars){}
NewNitsche::NewNitsche( const NewNitsche& nitschenew): Traditional(nitschenew)
{
  assert(0);
}
NewNitsche& NewNitsche::operator=( const NewNitsche& nitschenew)
{
  assert(0);
  Traditional::operator=(nitschenew);
  return *this;
}
std::string NewNitsche::getClassName() const
{
  return "NewNitsche";
}
/*--------------------------------------------------------------------------*/
const solver_options::pdepart::opts NewNitsche::setOptions()
{
  return solver_options::pdepart::cell;
}

/*--------------------------------------------------------------------------*/
void NewNitsche::setData(const alat::Map<std::string, int>& var2index, const solvers::Parameters& parameters)
{
  Traditional::setData(var2index, parameters);
  const solvers::P1* p1 = dynamic_cast<const solvers::P1*>((*_fems)[_ivar].get());
  assert(p1);
  _dofisbdry = &p1->getDofIsBdry();
  assert(_dofisbdry->size()==_meshinfo->nnodes);
  _symmetric = false;
}

void NewNitsche::_setisI(int iK)const
{
  _isI.set_size(_nlocal);
  assert(_isI.size()==_nlocal);
  _isI.fill(arma::fill::ones);
  for(int ii=0;ii<_nlocal;ii++)
  {
    if((*_dofisbdry)[_meshinfo->nodes_of_cells(ii,iK)])
    {
      _isI[ii] = 0;
    }
  }
}
/*--------------------------------------------------------------------------*/
void NewNitsche::computeRhsCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const
{
  Traditional::computeRhsCell(iK, floc, uloc);
  if(not _cellisbdry[iK])
  {
    return;
  }
  _setisI(iK);
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  alat::Node xK = _mesh->getNodeOfCell(iK);
  double diff = _localmodel->diffusion(xK.x(), xK.y(), xK.z());
  arma::vec f(_ncomp);
  arma::mat udir(_ncomp, _nlocal);
  arma::vec uhelp(_ncomp);
  for(int ii=0; ii<_nlocal;ii++)
  {
    int iN = _meshinfo->nodes_of_cells(ii,iK);
    _application->getDirichlet(_ivar)(uhelp, _meshinfo->nodes(0,iN), _meshinfo->nodes(1,iN), _meshinfo->nodes(2,iN));
    udir.col(ii) = uhelp;
  }
  arma::mat FuD(_ncomp,_nlocal);
  for(int ii=0; ii<_nlocal;ii++)
  {
    _localmodel->reaction(FuD.col(ii), udir.col(ii));
  }
  double moc=_meshinfo->measure_of_cells[iK];
  double scalemass = moc/(_meshinfo->dim+1.0);
  for(int ii=0; ii<_nlocal;ii++)
  {
    int iN = _meshinfo->nodes_of_cells(ii,iK);
    int iS = _meshinfo->sides_of_cells(ii,iK);
    if(_isI[ii])
    {
      if(_symmetric)
      {
        for(int jj=0; jj<_nlocal;jj++)
        {
          if(not _isI[jj])
          {
            int jN = _meshinfo->nodes_of_cells(jj,iK);
            int jS = _meshinfo->sides_of_cells(jj,iK);
            // diffusion
            for(int icomp=0;icomp<_ncomp;icomp++)
            {
              floc[_ivar](icomp,ii) -= diff*fem.laplace(ii,jj)*udir(icomp,jj);
              if(_lumpedmass)
              {
                if(ii==jj) floc[_ivar](icomp,ii) -= fem.mass_lumped[ii]*FuD(icomp, ii);
              }
              else
              {
                floc[_ivar](icomp,ii) -= fem.mass(ii,jj)*FuD(icomp, jj);
              }
            }
          }
        }
      }
    }
    else
    {
      _application->getRightHandSide(_ivar)(f, _meshinfo->nodes(0,iN), _meshinfo->nodes(1,iN), _meshinfo->nodes(2,iN));
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        floc[_ivar](icomp,ii) -= f[icomp]*scalemass;
      }
      for(int jj=0; jj<_nlocal;jj++)
      {
        if(not _isI[jj])
        {
          for(int icomp=0;icomp<_ncomp;icomp++)
          {
            floc[_ivar](icomp,ii) += diff*fem.laplace(ii,jj)*udir(icomp,jj);
            if(_lumpedmass)
            {
              if(ii==jj) floc[_ivar](icomp,ii) += fem.mass_lumped[ii]*FuD(icomp, ii);
            }
            else
            {
              floc[_ivar](icomp,ii) += fem.mass(ii,jj)*FuD(icomp, jj);
            }
          }
        }
      }
    }
  }
}
/*--------------------------------------------------------------------------*/
void NewNitsche::computeResidualCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const
{
  Traditional::computeResidualCell(iK, floc, uloc);
  if(not _cellisbdry[iK])
  {
    return;
  }

  _setisI(iK);
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  alat::Node xK = _mesh->getNodeOfCell(iK);
  double diff = _localmodel->diffusion(xK.x(), xK.y(), xK.z());
  arma::mat Fu(_ncomp,_nlocal);
  for(int ii=0; ii<_nlocal;ii++)
  {
    _localmodel->reaction(Fu.col(ii), uloc[_ivar].col(ii));
  }

  double moc=_meshinfo->measure_of_cells[iK];
  for(int ii=0; ii<_nlocal;ii++)
  {
    int iN = _meshinfo->nodes_of_cells(ii,iK);
    int iS = _meshinfo->sides_of_cells(ii,iK);
    for(int jj=0; jj<_nlocal;jj++)
    {
      int jN = _meshinfo->nodes_of_cells(jj,iK);
      int jS = _meshinfo->sides_of_cells(jj,iK);
      bool sub;
      if(_symmetric)
      {
        sub = (not _isI[ii] and _isI[jj]) or (_isI[ii] and not _isI[jj]);
      }
      else
      {
        sub = (not _isI[ii] and _isI[jj]);
      }
      if(sub)
      {
        for(int icomp=0;icomp<_ncomp;icomp++)
        {
          // diffusion
          floc[_ivar](icomp,ii) -= diff*fem.laplace(ii,jj)*uloc[_ivar](icomp, jj);
          // floc[_ivar](icomp,ii) -= d*uloc[_ivar](icomp, jj);
          if(_lumpedmass)
          {
            // if(ii==jj) floc[_ivar](icomp,ii) -= scalemass*Fu(icomp, ii);
            if(ii==jj) floc[_ivar](icomp,ii) += fem.mass_lumped[ii]*Fu(icomp, ii);
          }
          else
          {
            floc[_ivar](icomp,ii) -= fem.mass(ii,jj)*Fu(icomp, jj);
          }
        }
      }
    }
  }
}
/*--------------------------------------------------------------------------*/
void NewNitsche::computeMatrixCell(int iK, solvers::PdePartData::mat& mat, solvers::PdePartData::imat& mat_i, solvers::PdePartData::imat& mat_j, const solvers::PdePartData::vec& uloc)const
{
  Traditional::computeMatrixCell(iK, mat, mat_i, mat_j, uloc);
  if(not _cellisbdry[iK])
  {
    return;
  }

  _setisI(iK);
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  alat::Node xK = _mesh->getNodeOfCell(iK);
  arma::cube DF(_ncomp,_ncomp,_nlocal);
  for(int ii=0; ii<_nlocal;ii++)
  {
    _localmodel->reaction_d(DF.slice(ii), uloc[_ivar].col(ii));
  }
  double diff = _localmodel->diffusion(xK.x(), xK.y(), xK.z());
  double moc=_meshinfo->measure_of_cells[iK];

  int count=0;
  for(int ii=0; ii<_nlocal;ii++)
  {
    int iN = _meshinfo->nodes_of_cells(ii,iK);
    int iS = _meshinfo->sides_of_cells(ii,iK);
    for(int jj=0; jj<_nlocal;jj++)
    {
      bool sub;
      if(_symmetric)
      {
        sub = (not _isI[ii] and _isI[jj]) or (_isI[ii] and not _isI[jj]);
      }
      else
      {
        sub = (not _isI[ii] and _isI[jj]);
      }
      int jN = _meshinfo->nodes_of_cells(jj,iK);
      int jS = _meshinfo->sides_of_cells(jj,iK);
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        for(int jcomp=0;jcomp<_ncomp;jcomp++)
        {
          if(sub)
          {
            // diffusion
            if(icomp==jcomp)
            {
              mat(_ivar,_ivar)[count] -= diff*fem.laplace(ii,jj);
            }
            if(_lumpedmass)
            {
              if(ii==jj)
              {
                mat(_ivar,_ivar)[count] -= fem.mass_lumped[ii]*DF(icomp,jcomp, ii);
              }
            }
            else
            {
              mat(_ivar,_ivar)[count] -= fem.mass(ii,jj)*DF(icomp,jcomp, ii);
            }
          }
          count++;
        }
      }
    }
  }
}
