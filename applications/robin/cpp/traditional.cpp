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
  return solver_options::pdepart::cell + solver_options::pdepart::bdry;
}
/*--------------------------------------------------------------------------*/
void Traditional::setData(const alat::Map<std::string, int>& var2index, const solvers::Parameters& parameters)
{
  PdePart::setData(var2index, parameters);
  // _gamma=0.0;
  _lumpedmass = true;
}

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
}

/*--------------------------------------------------------------------------*/
void Traditional::computeResidualCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc) const
{
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  alat::Node xK = _mesh->getNodeOfCell(iK);
  // alat::armavec betavec(3, arma::fill::zeros);
  // (*_betafct)(betavec,xK.x(), xK.y(), xK.z(), _mesh->getDimension());
  // std::cerr << "betavec="<<betavec.t();

  arma::mat Fu(_ncomp,_nlocal);
  for(int ii=0; ii<_nlocal;ii++)
  {
    _localmodel->reaction(Fu.col(ii), uloc[_ivar].col(ii));
  }
  double diff = _localmodel->diffusion(xK.x(), xK.y(), xK.z());
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
}

/*--------------------------------------------------------------------------*/
void Traditional::computeMatrixCell(int iK, solvers::PdePartData::mat& mat, const solvers::PdePartData::vec& uloc)const
{
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  alat::Node xK = _mesh->getNodeOfCell(iK);
  arma::cube DF(_ncomp,_ncomp,_nlocal);
  for(int ii=0; ii<_nlocal;ii++)
  {
    _localmodel->reaction_d(DF.slice(ii), uloc[_ivar].col(ii));
  }
  double diff = _localmodel->diffusion(xK.x(), xK.y(), xK.z());
  int count=0;
  for(int ii=0; ii<_nlocal;ii++)
  {
    int iS = _meshinfo->sides_of_cells(ii,iK);
    double beta = (*_betavec)[iS]/_meshinfo->dim/(_meshinfo->dim+1)*_meshinfo->sigma(ii,iK);
    for(int jj=0; jj<_nlocal;jj++)
    {
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        for(int jcomp=0;jcomp<_ncomp;jcomp++)
        {
          if(icomp==jcomp)
          {
            // diffusion
            mat(_ivar,_ivar)(icomp*_nlocal+ii, jcomp*_nlocal+jj) += diff*fem.laplace(ii,jj);
            // convection
            mat(_ivar,_ivar)(icomp*_nlocal+ii, jcomp*_nlocal+jj) += beta;
          }
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
          count++;
        }
      }
    }
  }
}

/*--------------------------------------------------------------------------*/
void Traditional::computeRhsBdry(int color, int iK, int iS, int iil, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const
{
  alat::armavec neum(_ncomp), dir(_ncomp);
  const arma::subview_col<double> normal = _meshinfo->normals.col(iS);
  double dn = arma::norm(normal);
  double dns = 1./dn;
  alat::Node xS = _mesh->getNodeOfSide(iS);
  _application->getDirichlet(_ivar)(dir, xS.x(), xS.y(), xS.z());
  _application->getNeumann(_ivar)(neum, dns*normal[0], dns*normal[1], dns*normal[2], xS.x(), xS.y(), xS.z());

  arma::mat udir(_ncomp, _nlocal);
  alat::armavec uhelp(_ncomp);
  for(int ii=0; ii<_nlocal;ii++)
  {
    int iN = _meshinfo->nodes_of_cells(ii,iK);
    _application->getDirichlet(_ivar)(uhelp, _meshinfo->nodes(0,iN), _meshinfo->nodes(1,iN), _meshinfo->nodes(2,iN));
    udir.col(ii) = uhelp;
  }
  alat::Node xK = _mesh->getNodeOfCell(iK);
  double dnormal = dn/_meshinfo->dim;
  double moc=_meshinfo->measure_of_cells[iK];
  double beta = (*_betavec)[iS]/_meshinfo->dim;
  double diff = _localmodel->diffusion(xK.x(), xK.y(), xK.z());
  double h = moc/dn;
  double s =  1.0/(_gamma*h + 1.0/_localmodel->_robin);
  double seps =  s/_localmodel->_robin;
  double reps =  - _gamma*h*s;
  double repseps =  reps/_localmodel->_robin;
  double scalediff = diff/(_meshinfo->dim*moc);

  // seps = repseps = 0.0;

  // std::cerr << "s="<<s << " seps="<<seps << " reps="<<reps << " repseps="<<repseps<<"\n";

  for(int ii=0; ii<_nlocal;ii++)
  {
    int iN = _meshinfo->nodes_of_cells(ii,iK);
    int iSi = _meshinfo->sides_of_cells(ii,iK);
    double dni = scalediff*arma::dot(_meshinfo->normals.col(iSi), normal)*_meshinfo->sigma(ii,iK);
    for(int icomp=0;icomp<_ncomp;icomp++)
    {
      floc[_ivar][icomp*_nlocal + ii] -= reps*dir[icomp]*dni;
      floc[_ivar][icomp*_nlocal + ii] -= repseps*neum[icomp]*dni;
    }
    if(ii==iil) continue;
    for(int icomp=0;icomp<_ncomp;icomp++)
    {
      floc[_ivar][icomp*_nlocal + ii] += seps*neum[icomp]*dnormal;
      floc[_ivar][icomp*_nlocal + ii] += s*udir[icomp*_nlocal + ii]*dnormal;
      if(beta<0.0) floc[_ivar][icomp*_nlocal + ii] -= beta*udir[icomp*_nlocal + ii];
    }
  }
}
/*--------------------------------------------------------------------------*/
void Traditional::computeResidualBdry(int color, int iK, int iS, int iil, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const
{
  alat::Node xK = _mesh->getNodeOfCell(iK);
  double diff = _localmodel->diffusion(xK.x(), xK.y(), xK.z());
  const arma::subview_col<double> normal = _meshinfo->normals.col(iS);
  double dn = arma::norm(normal);
  double dnormal = dn/_meshinfo->dim;
  double moc=_meshinfo->measure_of_cells[iK];
  double beta = (*_betavec)[iS]/_meshinfo->dim;
  double h = moc/dn;
  double s =  1.0/(_gamma*h + 1.0/_localmodel->_robin);
  // double seps =  s/_localmodel->_robin;
  double reps =  - _gamma*h*s;
  double repseps =  reps/_localmodel->_robin;
  double scalediff = diff/(_meshinfo->dim*moc);

  // seps = repseps = 0.0;

  for(int ii=0; ii<_nlocal;ii++)
  {
    int iN = _meshinfo->nodes_of_cells(ii,iK);
    int iSi = _meshinfo->sides_of_cells(ii,iK);
    double dni = scalediff*arma::dot(_meshinfo->normals.col(iSi), normal)*_meshinfo->sigma(ii,iK);
    for(int jj=0; jj<_nlocal;jj++)
    {
      int iSj = _meshinfo->sides_of_cells(jj,iK);
      double dnj = scalediff*arma::dot(_meshinfo->normals.col(iSj), normal)*_meshinfo->sigma(jj,iK);
      if(jj!=iil)
      {
        for(int icomp=0;icomp<_ncomp;icomp++)
        {
          floc[_ivar][icomp*_nlocal + ii] -= reps*dni*uloc[_ivar][icomp*_nlocal+jj]/_meshinfo->dim;
        }
      }
      if(ii!=iil)
      {
        for(int icomp=0;icomp<_ncomp;icomp++)
        {
          floc[_ivar][icomp*_nlocal + ii] -= reps*dnj*uloc[_ivar][icomp*_nlocal+jj]/_meshinfo->dim;
        }
      }
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        floc[_ivar][icomp*_nlocal + ii] += repseps*dni*dnj*uloc[_ivar][icomp*_nlocal+jj]/dn;
      }
      if(not _lumpedmass)
      {
        if(jj!=iil and ii!=iil)
        {
          for(int icomp=0;icomp<_ncomp;icomp++)
          {
          if(ii==jj) floc[_ivar][icomp*_nlocal + ii] += s*uloc[_ivar][icomp*_nlocal+jj]*dn/3.0;
          else floc[_ivar][icomp*_nlocal + ii] += s*uloc[_ivar][icomp*_nlocal+jj]*dn/6.0;
        }
        }
      }
    }
    if(_lumpedmass)
    {
      if(ii==iil) continue;
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        floc[_ivar][icomp*_nlocal + ii] += s*uloc[_ivar][icomp*_nlocal+ii]*dnormal;
      }
    }
  }
}
/*--------------------------------------------------------------------------*/
void Traditional::computeMatrixBdry(int color, int iK, int iS, int iil, solvers::PdePartData::mat& mat, const solvers::PdePartData::vec& uloc)const
{
  alat::Node xK = _mesh->getNodeOfCell(iK);
  double diff = _localmodel->diffusion(xK.x(), xK.y(), xK.z());

  const arma::subview_col<double> normal = _meshinfo->normals.col(iS);
  double dn = arma::norm(normal);
  double dnormal = dn/_meshinfo->dim;
  double moc=_meshinfo->measure_of_cells[iK];
  double beta = (*_betavec)[iS]/_meshinfo->dim;
  double h = moc/dn;
  double s =  1.0/(_gamma*h + 1.0/_localmodel->_robin);
  // double seps =  s/_localmodel->_robin;
  double reps =  - _gamma*h*s;
  double repseps =  reps/_localmodel->_robin;
  double scalediff = diff/(_meshinfo->dim*moc);

  // seps = repseps = 0.0;

  int count=0;
  for(int ii=0; ii<_nlocal;ii++)
  {
    int iN = _meshinfo->nodes_of_cells(ii,iK);
    int iSi = _meshinfo->sides_of_cells(ii,iK);
    double dni = scalediff*arma::dot(_meshinfo->normals.col(iSi), normal)*_meshinfo->sigma(ii,iK);
    for(int jj=0; jj<_nlocal;jj++)
    {
      int jN = _meshinfo->nodes_of_cells(jj,iK);
      int iSj = _meshinfo->sides_of_cells(jj,iK);
      double dnj = scalediff*arma::dot(_meshinfo->normals.col(iSj), normal)*_meshinfo->sigma(jj,iK);
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        if(jj!=iil) mat(_ivar,_ivar)(icomp*_nlocal+ii, jcomp*_nlocal+jj) -= reps*dni/_meshinfo->dim;
        if(ii!=iil) mat(_ivar,_ivar)(icomp*_nlocal+ii, jcomp*_nlocal+jj) -= reps*dnj/_meshinfo->dim;
        mat(_ivar,_ivar)(icomp*_nlocal+ii, jcomp*_nlocal+jj) += repseps*dni*dnj/dn;
        if(_lumpedmass)
        {
          if(ii!=iil && ii==jj)
          {
            mat  (_ivar,_ivar)[count] += s*dnormal;
            if(beta>0.0) mat(_ivar,_ivar)(icomp*_nlocal+ii, jcomp*_nlocal+jj) += beta;
          }
        }
        else
        {
          if(jj!=iil and ii!=iil)
          {
            if(ii==jj)
            {
                mat  (_ivar,_ivar)[count] += s*dn/3.0;
                if(beta>0.0) mat(_ivar,_ivar)(icomp*_nlocal+ii, jcomp*_nlocal+jj) += beta/3.0;
            }
            else
            {
              mat  (_ivar,_ivar)[count] += s*dn/6.0;
              if(beta>0.0) mat(_ivar,_ivar)(icomp*_nlocal+ii, jcomp*_nlocal+jj) += beta/6.0;
            }
          }
        }
        count++;
      }
    }
  }
  // assert(count==sizemat);
}
