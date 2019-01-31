#include  "Alat/vectoronevariable.hpp"
#include  "Alat/sparsitypatternsoft.hpp"
#include  "Mesh/measureofcell.hpp"
#include  "Mesh/normals.hpp"
#include  "Mesh/nodesandnodesofcells.hpp"
#include  "Solvers/variable.hpp"
#include  "Solvers/meshunitwithdata.hpp"
#include  "traditional.hpp"
#include  "model.hpp"
#include  <cassert>

/*--------------------------------------------------------------------------*/
Traditional::~Traditional() {}
Traditional::Traditional(alat::StringList vars): PdePart(vars){}
Traditional::Traditional( const Traditional& pdepart): PdePart(pdepart)
{
  assert(0);
}
Traditional& Traditional::operator=( const Traditional& pdepart)
{
  assert(0);
  PdePart::operator=(pdepart);
  return *this;
}
std::string Traditional::getClassName() const
{
  return "Traditional";
}
/*--------------------------------------------------------------------------*/
const solver_options::pdepart::opts Traditional::setOptions()
{
  return solver_options::pdepart::cell;
}
// /*--------------------------------------------------------------------------*/
// void Traditional::setData(const alat::Map<std::string, int>& var2index, const solvers::Parameters& parameters)
// {
//   PdePart::setData(var2index, parameters);
// }

/*--------------------------------------------------------------------------*/
void Traditional::computeRhsCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const
{
  if(not _lumpedmass)
  {
    PdePart::computeRhsCell(iK, floc, uloc);
    return;
  }
  alat::armavec f(_ncomp);
  double moc=_meshinfo->measure_of_cells[iK];
  double scalemass = moc/(_meshinfo->dim+1.0);

  for(int ii=0; ii<_nlocal;ii++)
  {
    int iN = _meshinfo->nodes_of_cells(ii,iK);
    _application->getRightHandSide(_ivar)(f, _meshinfo->nodes(0,iN), _meshinfo->nodes(1,iN), _meshinfo->nodes(2,iN));
    for(int icomp=0;icomp<_ncomp;icomp++)
    {
      floc[_ivar][icomp*_nlocal + ii] += f[icomp]*scalemass;
    }
  }
  // std::cerr << "computeRhsCell floc[_ivar]="<<floc[_ivar];
}

/*--------------------------------------------------------------------------*/
void Traditional::computeResidualCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc) const
{
  // std::cerr << "computeResidualCell floc[_ivar]="<<floc[_ivar];
  // std::cerr << "computeResidualCell uloc[_ivar]="<<uloc[_ivar];
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  alat::Node xK = _mesh->getNodeOfCell(iK);
  // alat::armavec betavec(3, arma::fill::zeros);
  // (*_betafct)(betavec,xK.x(), xK.y(), xK.z(), _mesh->getDimension());
  // std::cerr << "betavec="<<betavec.t();

  alat::armamat Fu(_ncomp,_nlocal);
  alat::armavec u(_ncomp), f(_ncomp);
  for(int ii=0; ii<_nlocal;ii++)
  {
    for(int icomp=0;icomp<_ncomp;icomp++)
    {
      u[icomp] = uloc[_ivar][icomp*_nlocal + ii];
    }
    _localmodel->reaction(f, u);
    Fu.col(ii) = f;
  }
  // std::cerr << "computeResidualCell Fu="<<Fu;
  double diff = _localmodel->diffusion(xK.x(), xK.y(), xK.z());
  // std::cerr << "computeResidualCell diff="<<diff;
  for(int ii=0; ii<_nlocal;ii++)
  {
    int iS = _meshinfo->sides_of_cells(ii,iK);
    double beta = (*_betavec)[iS]/_meshinfo->dim/(_meshinfo->dim+1)*_meshinfo->sigma(ii,iK);
    for(int icomp=0;icomp<_ncomp;icomp++)
    {
      // reaction
      if(_lumpedmass)
      {
        floc[_ivar][icomp*_nlocal + ii] += fem.mass_lumped[ii]*Fu(icomp, ii);
      }
    }
    for(int jj=0; jj<_nlocal;jj++)
    {
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        // diffusion
        floc[_ivar][icomp*_nlocal + ii] += diff*fem.laplace(ii,jj)*uloc[_ivar][icomp*_nlocal+jj];
        // convection
        floc[_ivar][icomp*_nlocal + ii] += beta*uloc[_ivar][icomp*_nlocal+jj];
        if(not _lumpedmass)
        {
          floc[_ivar][icomp*_nlocal + ii] += fem.mass(ii,jj)*Fu(icomp, jj);
        }
      }
    }
  }
  // std::cerr << "computeResidualCell floc[_ivar]="<<floc[_ivar];
}

/*--------------------------------------------------------------------------*/
void Traditional::computeMatrixCell(int iK, solvers::PdePartData::mat& mat, const solvers::PdePartData::vec& uloc)const
{
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  alat::Node xK = _mesh->getNodeOfCell(iK);
  arma::cube DF(_ncomp,_ncomp,_nlocal);
  alat::armamat DF1(_ncomp,_ncomp);
  alat::armavec u(_ncomp);
  for(int ii=0; ii<_nlocal;ii++)
  {
    for(int icomp=0;icomp<_ncomp;icomp++)
    {
      u[icomp] = uloc[_ivar][icomp*_nlocal + ii];
    }
    _localmodel->reaction_d(DF1, u);
    DF.slice(ii) = DF1;
  }
  double diff = _localmodel->diffusion(xK.x(), xK.y(), xK.z());
  for(int ii=0; ii<_nlocal;ii++)
  {
    int iS = _meshinfo->sides_of_cells(ii,iK);
    double beta = (*_betavec)[iS]/_meshinfo->dim/(_meshinfo->dim+1)*_meshinfo->sigma(ii,iK);
    for(int jj=0; jj<_nlocal;jj++)
    {
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        // diffusion
        mat(_ivar,_ivar)(icomp*_nlocal+ii, icomp*_nlocal+jj) += diff*fem.laplace(ii,jj);
        // convection
        mat(_ivar,_ivar)(icomp*_nlocal+ii, icomp*_nlocal+jj) += beta;
        for(int jcomp=0;jcomp<_ncomp;jcomp++)
        {
          // mass
          if(_lumpedmass)
          {
            if(ii==jj)
            {
              mat(_ivar,_ivar)(icomp*_nlocal+ii, jcomp*_nlocal+jj) += fem.mass_lumped[ii]*DF(icomp,jcomp, ii);
            }
          }
          else
          {
            mat(_ivar,_ivar)(icomp*_nlocal+ii, jcomp*_nlocal+jj) += fem.mass(ii,jj)*DF(icomp,jcomp, ii);
          }
        }
      }
    }
  }
}

/*--------------------------------------------------------------------------*/
void Traditional::computeMatrixCellOld(int iK, solvers::PdePartData::mat& mat, const solvers::PdePartData::vec& uloc)const
{
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  alat::Node xK = _mesh->getNodeOfCell(iK);
  arma::cube DF(_ncomp,_ncomp,_nlocal);
  for(int ii=0; ii<_nlocal;ii++)
  {
    _localmodel->reaction_d(DF.slice(ii), uloc[_ivar].col(ii));
  }
  double diff = _localmodel->diffusion(xK.x(), xK.y(), xK.z());

  double moc=_meshinfo->measure_of_cells[iK];
  double scalediff = 1.0/(_meshinfo->dim*_meshinfo->dim*moc);
  double scalemass = moc/(_meshinfo->dim+1.0);

  int count=0;
  for(int ii=0; ii<_nlocal;ii++)
  {
    int iN = _meshinfo->nodes_of_cells(ii,iK);
    int iS = _meshinfo->sides_of_cells(ii,iK);
    for(int jj=0; jj<_nlocal;jj++)
    {
      int jN = _meshinfo->nodes_of_cells(jj,iK);
      int jS = _meshinfo->sides_of_cells(jj,iK);
      double dot = arma::dot(_meshinfo->normals.col(iS), _meshinfo->normals.col(jS));
      double d = diff*dot*scalediff*_meshinfo->sigma(ii,iK)*_meshinfo->sigma(jj,iK);
      // double d2 = diff*moc*( arma::dot(fem.dphi.col(ii),fem.dphi.col(jj)) );
      // std::cerr << "d="<<d << " d2="<<d2<<"\n";
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        for(int jcomp=0;jcomp<_ncomp;jcomp++)
        {
          // diffusion
          if(icomp==jcomp)
          {
            mat(_ivar,_ivar)(icomp*_nlocal+ii, jcomp*_nlocal+jj) += d;
          }
          // mass
          if(ii==jj)
          {
            mat(_ivar,_ivar)(icomp*_nlocal+ii, jcomp*_nlocal+jj) += scalemass*DF(icomp,jcomp, ii);
          }
          count++;
        }
      }
    }
  }
}
