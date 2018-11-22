#include  "Alat/sparsitypatternsoft.hpp"
#include  "Alat/matrixonevariableinterface.hpp"
#include  "Mesh/nodescellsweight.hpp"
#include  "Solvers/pdepartp1.hpp"
#include  "Solvers/meshunitwithdata.hpp"
#include  "Solvers/variable.hpp"
#include  <cassert>

using namespace solvers;

/*--------------------------------------------------------------------------*/
PdePartP1::~PdePartP1() {}
PdePartP1::PdePartP1(alat::StringList vars): PdePartInterface(vars){}
PdePartP1::PdePartP1( const PdePartP1& pdepartc1): PdePartInterface(pdepartc1)
{
  assert(0);
}
PdePartP1& PdePartP1::operator=( const PdePartP1& pdepartc1)
{
  assert(0);
  PdePartInterface::operator=(pdepartc1);
  return *this;
}
std::string PdePartP1::getClassName() const
{
  return "PdePartP1";
}
/*--------------------------------------------------------------------------*/
void PdePartP1::setData(const alat::Map<std::string, int>& var2index, const solvers::Parameters& parameters)
{
  assert(_fems->size()==1);
  _ivar = 0;
  _ncomp = (*_fems)[_ivar]->getNcomp();
}
/*--------------------------------------------------------------------------*/
void PdePartP1::_laplace(int iK, arma::mat& mat)const
{
  double moc=_meshinfo->measure_of_cells[iK];
  double scale = 1.0/(_meshinfo->dim*_meshinfo->dim*moc);
  assert(mat.n_cols==_meshinfo->nnodespercell);
  assert(mat.n_rows==_meshinfo->nnodespercell);
  for(int ii=0; ii<_meshinfo->nnodespercell;ii++)
  {
    int iN = _meshinfo->nodes_of_cells(ii,iK);
    int iS = _meshinfo->sides_of_cells(ii,iK);
    for(int jj=0; jj<_meshinfo->nnodespercell;jj++)
    {
      int jN = _meshinfo->nodes_of_cells(jj,iK);
      int jS = _meshinfo->sides_of_cells(jj,iK);
      double dot = arma::dot(_meshinfo->normals.col(iS), _meshinfo->normals.col(jS));
      mat(ii,jj) = dot*scale*_meshinfo->sigma(ii,iK)*_meshinfo->sigma(jj,iK);
    }
  }
}
/*--------------------------------------------------------------------------*/
void PdePartP1::_massMatrix(int iK, arma::mat& mat, bool unscaled)const
{
  double moc=_meshinfo->measure_of_cells[iK];
  double scale;
  if(unscaled) {scale = 1.0/(_meshinfo->dim+2.0)/(_meshinfo->dim+1.0);}
  else {scale = moc/(_meshinfo->dim+2.0)/(_meshinfo->dim+1.0);}
  for(int ii=0; ii<_meshinfo->nnodespercell;ii++)
  {
    for(int jj=0; jj<_meshinfo->nnodespercell;jj++)
    {
      mat(ii,jj) = scale;
    }
  }
  for(int ii=0; ii<_meshinfo->nnodespercell;ii++) {mat(ii,ii) *= 2.0;}
}
/*--------------------------------------------------------------------------*/
void PdePartP1::computeRhsCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const
{
  if(not _application->hasRightHandSide(_ivar)) return;
  arma::vec f(_ncomp);
  double moc=_meshinfo->measure_of_cells[iK];
  double scale = moc/(_meshinfo->dim+1.0);
  for(int ii=0; ii<_meshinfo->nnodespercell;ii++)
  {
    int iN = _meshinfo->nodes_of_cells(ii,iK);
    _application->getRightHandSide(_ivar)(f, _meshinfo->nodes(0,iN), _meshinfo->nodes(1,iN), _meshinfo->nodes(2,iN));
    for(int icomp=0;icomp<_ncomp;icomp++)
    {
      floc[_ivar](icomp,ii) += f[icomp]*scale;
    }
  }
}

/*--------------------------------------------------------------------------*/
void PdePartP1::computeResidualCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc) const
{
  double moc=_meshinfo->measure_of_cells[iK];
  double scalediff = 1.0/(_meshinfo->dim*_meshinfo->dim*moc);
  for(int ii=0; ii<_meshinfo->nnodespercell;ii++)
  {
    int iN = _meshinfo->nodes_of_cells(ii,iK);
    int iS = _meshinfo->sides_of_cells(ii,iK);
    for(int jj=0; jj<_meshinfo->nnodespercell;jj++)
    {
      int jN = _meshinfo->nodes_of_cells(jj,iK);
      int jS = _meshinfo->sides_of_cells(jj,iK);
      double dot = arma::dot(_meshinfo->normals.col(iS), _meshinfo->normals.col(jS));
      double d = dot*scalediff*_meshinfo->sigma(ii,iK)*_meshinfo->sigma(jj,iK);
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        floc[_ivar](icomp,ii) += d*uloc[_ivar][icomp*_meshinfo->nnodespercell + jj];
      }
    }
  }
  // std::cerr << "PdePartP1::computeResidualCell() floc[_ivar]="<<floc[_ivar]<<"\n";
}


/*--------------------------------------------------------------------------*/
void PdePartP1::computeMatrixCell(int iK, solvers::PdePartData::mat& mat, solvers::PdePartData::imat& mat_i, solvers::PdePartData::imat& mat_j, const solvers::PdePartData::vec& uloc)const
{
  double moc=_meshinfo->measure_of_cells[iK];
  double scalediff = 1.0/(_meshinfo->dim*_meshinfo->dim*moc);
  int count=0;
  for(int ii=0; ii<_meshinfo->nnodespercell;ii++)
  {
    int iN = _meshinfo->nodes_of_cells(ii,iK);
    int iS = _meshinfo->sides_of_cells(ii,iK);
    for(int jj=0; jj<_meshinfo->nnodespercell;jj++)
    {
      int jN = _meshinfo->nodes_of_cells(jj,iK);
      int jS = _meshinfo->sides_of_cells(jj,iK);
      double dot = arma::dot(_meshinfo->normals.col(iS), _meshinfo->normals.col(jS));
      double d = dot*scalediff*_meshinfo->sigma(ii,iK)*_meshinfo->sigma(jj,iK);
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        mat(_ivar,_ivar)[count] = d;
        // mat_i(_ivar,_ivar)[count] = icomp*_meshinfo->nnodes + iN;
        // mat_j(_ivar,_ivar)[count] = icomp*_meshinfo->nnodes + jN;
        count++;
      }
    }
  }
  // std::cerr << "Aloc = " << Aloc.t();
}

/*--------------------------------------------------------------------------*/
void PdePartP1::computeRhsBdry(int color, int iK, int iS, int iil, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const{}
void PdePartP1::computeResidualBdry(int color, int iK, int iS, int iil, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const{}
void PdePartP1::computeMatrixBdry(int color, int iK, int iS, int iil, solvers::PdePartData::mat& mat, solvers::PdePartData::imat& mat_i, solvers::PdePartData::imat& mat_j, const solvers::PdePartData::vec& uloc)const{}
