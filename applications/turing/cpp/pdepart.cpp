#include  "Alat/vectoronevariable.hpp"
#include  "Alat/sparsitypatternsoft.hpp"
#include  "Mesh/measureofcell.hpp"
#include  "Mesh/normals.hpp"
#include  "Mesh/nodescellsweight.hpp"
#include  "Mesh/nodesandnodesofcells.hpp"
#include  "Solvers/variable.hpp"
#include  "Solvers/meshunitwithdata.hpp"
#include  "pdepart.hpp"
#include  "model.hpp"
#include  <cassert>

/*--------------------------------------------------------------------------*/
PdePart::~PdePart() {}
PdePart::PdePart(alat::StringList vars): solvers::PdePartP1(vars){}
PdePart::PdePart( const PdePart& pdepart): solvers::PdePartP1(pdepart)
{
  assert(0);
}
PdePart& PdePart::operator=( const PdePart& pdepart)
{
  assert(0);
  solvers::PdePartP1::operator=(pdepart);
  return *this;
}
std::string PdePart::getClassName() const
{
  return "PdePart";
}

/*--------------------------------------------------------------------------*/
void PdePart::setData(const alat::Map<std::string, int>& var2index, const solvers::Parameters& parameters)
{
  solvers::PdePartP1::setData(var2index, parameters);
  const Model* model = dynamic_cast<const Model*>(_model);
  assert(model);
  _k.set_size(2);
  _a = model->a;
  _b = model->b;
  _k[0] = model->k0;
  _k[1] = model->k1;
}
/*--------------------------------------------------------------------------*/
void PdePart::brussel(arma::subview_col<double> f, const arma::subview_col<double> u)const
{
  f[1] = u[0]*(_b-u[0]*u[1]);
  f[0] = _a - f[1] - u[0];
}
void PdePart::brussel_d(arma::mat& df, const arma::subview_col<double> u)const
{
  df(1,0) = _b-2.0*u[0]*u[1];
  df(1,1) = -u[0]*u[0];
  df(0,0) = -df(1,0)-1.0;
  df(0,1) = -df(1,1);
}

/*--------------------------------------------------------------------------*/
void PdePart::computeRhsCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const
{
  arma::mat Fu(_ncomp,_meshinfo->nnodespercell);
  double moc=_meshinfo->measure_of_cells[iK];
  double scalediff = 1.0/(_meshinfo->dim*_meshinfo->dim*moc);
  double drhs = moc/(_meshinfo->dim+1.0);
  double mass = moc/(_meshinfo->dim+2.0)/(_meshinfo->dim+1.0);
  double lmass = moc/(_meshinfo->dim+1.0);
  double dmass = moc/(_meshinfo->dim+2.0)/(_meshinfo->dim+1.0)/_dt;
  for(int ii=0; ii<_meshinfo->nnodespercell;ii++)
  {
    brussel(Fu.col(ii), uloc[_ivar].col(ii));
  }
  for(int ii=0; ii<_meshinfo->nnodespercell;ii++)
  {
    int iN = _meshinfo->nodes_of_cells(ii,iK);
    int iS = _meshinfo->sides_of_cells(ii,iK);
    for(int icomp=0;icomp<_ncomp;icomp++)
    {
      // mass
      floc[_ivar](icomp,ii) += (lmass/_dt)*uloc[_ivar](icomp, ii);
      // reaction
      floc[_ivar](icomp,ii) += 0.5*lmass*Fu(icomp, ii);

      for(int jj=0; jj<_meshinfo->nnodespercell;jj++)
      {
        int jN = _meshinfo->nodes_of_cells(jj,iK);
        int jS = _meshinfo->sides_of_cells(jj,iK);
        double dot = arma::dot(_meshinfo->normals.col(iS), _meshinfo->normals.col(jS));
        double d = dot*scalediff*_meshinfo->sigma(ii,iK)*_meshinfo->sigma(jj,iK);
        floc[_ivar](icomp,ii) -= 0.5*_k[icomp]*d*uloc[_ivar](icomp, jj);
      }
    }
  }
}

/*--------------------------------------------------------------------------*/
void PdePart::computeResidualCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc) const
{
  // std::cerr << "res _dt="<<_dt << " ";
  // std::cerr << "res uloc="<<uloc[0];
  arma::mat Fu(_ncomp,_meshinfo->nnodespercell);
  double moc=_meshinfo->measure_of_cells[iK];
  double scalediff = 1.0/(_meshinfo->dim*_meshinfo->dim*moc);
  double lmass = moc/(_meshinfo->dim+1.0);
  double mass = moc/(_meshinfo->dim+2.0)/(_meshinfo->dim+1.0);
  double dmass = moc/(_meshinfo->dim+2.0)/(_meshinfo->dim+1.0)/_dt;
  for(int ii=0; ii<_meshinfo->nnodespercell;ii++)
  {
    brussel(Fu.col(ii), uloc[_ivar].col(ii));
  }
  for(int ii=0; ii<_meshinfo->nnodespercell;ii++)
  {
    int iN = _meshinfo->nodes_of_cells(ii,iK);
    int iS = _meshinfo->sides_of_cells(ii,iK);
    // mass
    for(int icomp=0;icomp<_ncomp;icomp++)
    {
      floc[_ivar](icomp,ii) += (lmass/_dt)*uloc[_ivar](icomp, ii);
      // reaction
      floc[_ivar](icomp,ii) -= 0.5*lmass*Fu(icomp, ii);

      for(int jj=0; jj<_meshinfo->nnodespercell;jj++)
      {
        int jN = _meshinfo->nodes_of_cells(jj,iK);
        int jS = _meshinfo->sides_of_cells(jj,iK);
        double dot = arma::dot(_meshinfo->normals.col(iS), _meshinfo->normals.col(jS));
        double d = dot*scalediff*_meshinfo->sigma(ii,iK)*_meshinfo->sigma(jj,iK);
        // mass
        // if(ii==jj) {floc[0](icomp,ii) += 2.0*dmass*uloc[0](icomp, jj);}
        // else {floc[0](icomp,ii) += dmass*uloc[0](icomp, jj);}
        // reaction
        // if(ii==jj) {floc[0](icomp,ii) -= mass*Fu(icomp, jj);}
        // else {floc[0](icomp,ii) -= 0.5*mass*Fu(icomp, jj);}
        // diffusion
        floc[_ivar](icomp,ii) += 0.5*_k[icomp]*d*uloc[_ivar](icomp, jj);
      }
    }
  }
}

/*--------------------------------------------------------------------------*/
void PdePart::computeMatrixCell(int iK, solvers::PdePartData::mat& mat, solvers::PdePartData::imat& mat_i, solvers::PdePartData::imat& mat_j, const solvers::PdePartData::vec& uloc)const
{
  // std::cerr << "mat _dt="<<_dt << " ";
  // std::cerr << "mat uloc="<<uloc[0];
  int nloccell = _meshinfo->nnodespercell;
  int sizemat = _ncomp*nloccell*nloccell + _ncomp*_ncomp*nloccell;
  mat  (_ivar,_ivar).set_size(sizemat);
  mat_i(_ivar,_ivar).set_size(sizemat);
  mat_j(_ivar,_ivar).set_size(sizemat);
  arma::cube DF(_ncomp,_ncomp,_meshinfo->nnodespercell);
  double moc=_meshinfo->measure_of_cells[iK];
  double scalediff = 1.0/(_meshinfo->dim*_meshinfo->dim*moc);
  double lmass = moc/(_meshinfo->dim+1.0);
  double mass = moc/(_meshinfo->dim+2.0)/(_meshinfo->dim+1.0);
  double dmass = moc/(_meshinfo->dim+2.0)/(_meshinfo->dim+1.0)/_dt;
  for(int ii=0; ii<_meshinfo->nnodespercell;ii++)
  {
    brussel_d(DF.slice(ii), uloc[_ivar].col(ii));
  }
  int count=0;
  for(int ii=0; ii<nloccell;ii++)
  {
    int iN = _meshinfo->nodes_of_cells(ii,iK);
    int iS = _meshinfo->sides_of_cells(ii,iK);
    // reaction
    for(int icomp=0;icomp<_ncomp;icomp++)
    {
      for(int jcomp=0;jcomp<_ncomp;jcomp++)
      {
        // mass
        if(icomp==jcomp)
        {
          mat(_ivar,_ivar)[count] += (lmass/_dt);
        }
        mat  (_ivar,_ivar)[count] -= 0.5*lmass*DF(icomp,jcomp, ii);
        mat_i(_ivar,_ivar)[count] = icomp*_meshinfo->nnodes + iN;
        mat_j(_ivar,_ivar)[count] = jcomp*_meshinfo->nnodes + iN;
        count++;
      }
    }

    for(int jj=0; jj<nloccell;jj++)
    {
      int jN = _meshinfo->nodes_of_cells(jj,iK);
      int jS = _meshinfo->sides_of_cells(jj,iK);
      double dot = arma::dot(_meshinfo->normals.col(iS), _meshinfo->normals.col(jS));
      double d = dot*scalediff*_meshinfo->sigma(ii,iK)*_meshinfo->sigma(jj,iK);
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        // diffusion
        mat  (_ivar,_ivar)[count] += 0.5*_k[icomp]*d;
        mat_i(_ivar,_ivar)[count] = icomp*_meshinfo->nnodes + iN;
        mat_j(_ivar,_ivar)[count] = icomp*_meshinfo->nnodes + jN;
        count++;
      }
    }
  }
  assert(count==sizemat);
  // std::cerr << "Aloc = " << Aloc.t();
}
