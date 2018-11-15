#include  "Alat/vectoronevariable.hpp"
#include  "Alat/sparsitypatternsoft.hpp"
#include  "Mesh/measureofcell.hpp"
#include  "Mesh/normals.hpp"
#include  "Mesh/nodesandnodesofcells.hpp"
#include  "Solvers/variable.hpp"
#include  "Solvers/meshunitwithdata.hpp"
#include  "nitsche.hpp"
#include  "model.hpp"
#include  <cassert>

/*--------------------------------------------------------------------------*/
Nitsche::~Nitsche() {}
Nitsche::Nitsche(alat::StringList vars): Traditional(vars){}
Nitsche::Nitsche( const Nitsche& nitsche): Traditional(nitsche)
{
  assert(0);
}
Nitsche& Nitsche::operator=( const Nitsche& nitsche)
{
  assert(0);
  Traditional::operator=(nitsche);
  return *this;
}
std::string Nitsche::getClassName() const
{
  return "Nitsche";
}
/*--------------------------------------------------------------------------*/
const solver_options::pdepart::opts Nitsche::setOptions()
{
  _gamma = 10.0;
  return solver_options::pdepart::cell + solver_options::pdepart::bdry;
}

/*--------------------------------------------------------------------------*/
void Nitsche::computeRhsBdry(int color, int iK, int iS, int iil, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const
{
  arma::vec dir(_ncomp), dir2(_ncomp);

  const arma::subview_col<double> normal = _meshinfo->normals.col(iS);
  alat::Node xS = _mesh->getNodeOfSide(iS);
  _application->getDirichlet(_ivar)(dir, xS.x(), xS.y(), xS.z());


  arma::mat udir(_ncomp, _nlocal);
  arma::vec uhelp(_ncomp);
  for(int ii=0; ii<_nlocal;ii++)
  {
    int iN = _meshinfo->nodes_of_cells(ii,iK);
    _application->getDirichlet(_ivar)(uhelp, _meshinfo->nodes(0,iN), _meshinfo->nodes(1,iN), _meshinfo->nodes(2,iN));
    udir.col(ii) = uhelp;
  }


  alat::Node xK = _mesh->getNodeOfCell(iK);
  double diff = _localmodel->diffusion(xK.x(), xK.y(), xK.z());

  double moc=_meshinfo->measure_of_cells[iK];
  double scalediff = diff/(_meshinfo->dim*moc);
  double scalegamma = _gamma*diff*arma::dot(normal,normal)/(_meshinfo->dim*_meshinfo->dim*moc);
  double beta = (*_betavec)[iS]/_meshinfo->dim;

  for(int ii=0; ii<_nlocal;ii++)
  {
    int iN = _meshinfo->nodes_of_cells(ii,iK);
    int iSi = _meshinfo->sides_of_cells(ii,iK);
    double dni = scalediff*arma::dot(_meshinfo->normals.col(iSi), normal)*_meshinfo->sigma(ii,iK);
    for(int icomp=0;icomp<_ncomp;icomp++)
    {
      floc[_ivar](icomp,ii) += dir[icomp]*dni;
    }
    if(ii==iil) continue;
    for(int icomp=0;icomp<_ncomp;icomp++)
    {
      floc[_ivar](icomp,ii) += udir(icomp,ii)*scalegamma;
      if(beta<0.0) floc[_ivar](icomp,ii) -= beta*udir(icomp,ii);
    }
  }
}
/*--------------------------------------------------------------------------*/
void Nitsche::computeResidualBdry(int color, int iK, int iS, int iil, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const
{
  alat::Node xK = _mesh->getNodeOfCell(iK);
  double diff = _localmodel->diffusion(xK.x(), xK.y(), xK.z());

  const arma::subview_col<double> normal = _meshinfo->normals.col(iS);

  double moc=_meshinfo->measure_of_cells[iK];
  double scalediff = diff/(_meshinfo->dim*moc);
  double scalegamma = _gamma*diff*arma::dot(normal,normal)/(_meshinfo->dim*_meshinfo->dim*moc);

  double beta = (*_betavec)[iS]/_meshinfo->dim;

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
          floc[_ivar](icomp,ii) += dni*uloc[_ivar](icomp, jj)/_meshinfo->dim;
        }
      }
      if(ii!=iil)
      {
        for(int icomp=0;icomp<_ncomp;icomp++)
        {
          floc[_ivar](icomp,ii) += dnj*uloc[_ivar](icomp, jj)/_meshinfo->dim;
          if(beta>0.0) if(ii==jj) floc[_ivar](icomp,ii) += beta*uloc[_ivar](icomp, ii);
        }
      }
    }
    if(ii==iil) continue;
    for(int icomp=0;icomp<_ncomp;icomp++)
    {
      floc[_ivar](icomp,ii) += uloc[_ivar](icomp, ii)*scalegamma;
    }
  }
}
/*--------------------------------------------------------------------------*/
void Nitsche::computeMatrixBdry(int color, int iK, int iS, int iil, solvers::PdePartData::mat& mat, solvers::PdePartData::imat& mat_i, solvers::PdePartData::imat& mat_j, const solvers::PdePartData::vec& uloc)const
{
  alat::Node xK = _mesh->getNodeOfCell(iK);
  double diff = _localmodel->diffusion(xK.x(), xK.y(), xK.z());

  const arma::subview_col<double> normal = _meshinfo->normals.col(iS);

  double moc=_meshinfo->measure_of_cells[iK];
  double scalediff = diff/(_meshinfo->dim*moc);
  double scalegamma = _gamma*diff*arma::dot(normal,normal)/(_meshinfo->dim*_meshinfo->dim*moc);

  double beta = (*_betavec)[iS]/_meshinfo->dim;

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
        if(jj!=iil) mat(_ivar,_ivar)[count] += dni/_meshinfo->dim;
        if(ii!=iil) mat(_ivar,_ivar)[count] += dnj/_meshinfo->dim;
        if(ii!=iil && ii==jj)
        {
          mat  (_ivar,_ivar)[count] += scalegamma;
          if(beta>0.0) mat(_ivar,_ivar)[count] += beta;
        }
        count++;
      }
    }
  }
  // assert(count==sizemat);
}
