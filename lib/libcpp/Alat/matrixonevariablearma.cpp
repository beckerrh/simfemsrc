#include  "Alat/matrixonevariablearma.hpp"
#include  "Alat/vectoronevariableinterface.hpp"
#include  <cassert>

using namespace alat;

/*--------------------------------------------------------------------------*/
MatrixOneVariableArma::~MatrixOneVariableArma() {}
MatrixOneVariableArma::MatrixOneVariableArma(): arma::sp_mat(), alat::MatrixOneVariableInterface(){}
MatrixOneVariableArma::MatrixOneVariableArma( const MatrixOneVariableArma& variablematrix): arma::sp_mat(variablematrix), alat::MatrixOneVariableInterface(variablematrix)
{
  (*this).operator=(variablematrix);
}
MatrixOneVariableArma& MatrixOneVariableArma::operator=( const MatrixOneVariableArma& variablematrix)
{
  assert(0);
  arma::sp_mat::operator=(variablematrix);
  return *this;
}
std::string MatrixOneVariableArma::getClassName() const
{
  return "MatrixOneVariableArma";
}
MatrixOneVariableArma* MatrixOneVariableArma::clone() const
{
  return new MatrixOneVariableArma(*this);
}
bool MatrixOneVariableArma::needsConnectivity() const {return false;}

/*--------------------------------------------------------------------------*/
void MatrixOneVariableArma::set_size(int n, int m)
{
  arma::sp_mat::set_size( n, m);
}


/*--------------------------------------------------------------------------*/
void MatrixOneVariableArma::save(std::ostream& os, arma::file_type datatype) const
{
  static_cast<const arma::sp_mat&>(*this).save(os, datatype);
}

/*--------------------------------------------------------------------------*/
void MatrixOneVariableArma::fillzeros()
{
  arma::sp_mat::zeros();
}
/*--------------------------------------------------------------------------*/
void MatrixOneVariableArma::solve(VectorOneVariableInterface* u, const VectorOneVariableInterface* f)
{
  arma::Col<double>* uu = dynamic_cast<arma::Col<double>*>(u);
  assert(uu);
  const arma::Col<double>* ff = dynamic_cast<const arma::Col<double>*>(f);
  assert(ff);
  arma::spsolve(*uu, *this, *ff);
}
/*--------------------------------------------------------------------------*/
void MatrixOneVariableArma::rowIdentity(int index)
{
  for(arma::sp_mat::row_iterator p = arma::sp_mat::begin_row(index);p!=arma::sp_mat::end_row(index);++p)
  {
    if(p.col()==index) *p=1.0;
    else *p=0.0;
  }
}
/*--------------------------------------------------------------------------*/
void MatrixOneVariableArma::assemble(const alat::armamat& Alocal, const alat::armaivec& indicesi, const alat::armaivec& indicesj)
{
  for(int ii=0;ii<Alocal.size();ii++)
  {
    int i = indicesi[ii];
    int j = indicesj[ii];
    ( *this )(i,j) += Alocal[ii];
  }
}
// void MatrixOneVariableArma::assemble(const alat::armaivec& indicesi, const alat::armaivec& indicesj, const arma::mat& Alocal)
// {
//   int nloci = indicesi.size();
//   int nlocj = indicesj.size();
//   for(int ii = 0; ii < nloci; ii++)
//   {
//     int i = indicesi[ii];
//     for(int jj = 0; jj < nlocj; jj++)
//     {
//       int j = indicesj[jj];
//       ( *this )(i,j) += Alocal(ii,jj);
//     }
//   }
// }
/*--------------------------------------------------------------------------*/
void MatrixOneVariableArma::matrixVectorProduct(alat::VectorOneVariableInterface* out, const alat::VectorOneVariableInterface* in, double d) const
{
  alat::armavec* outv = dynamic_cast<alat::armavec*>(out); assert(outv);
  const alat::armavec* inv = dynamic_cast<const alat::armavec*>(in); assert(inv);

  (*outv) += d*static_cast<const arma::sp_mat&>(*this)*(*inv);
}
