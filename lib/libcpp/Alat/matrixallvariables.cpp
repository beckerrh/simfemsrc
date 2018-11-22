#include  "Alat/matrixallvariables.hpp"
#include  "Alat/vectorallvariables.hpp"
#include  "Alat/sparsitypatternsoft.hpp"
#include  "Alat/strings.hpp"
#include  <cassert>

using namespace alat;

/*--------------------------------------------------------------------------*/
MatrixAllVariables::~MatrixAllVariables(){}
MatrixAllVariables::MatrixAllVariables() : alat::Matrix<std::shared_ptr<alat::MatrixOneVariableInterface> >(){}
MatrixAllVariables::MatrixAllVariables(int nvars, int mvars) : alat::Matrix<std::shared_ptr<alat::MatrixOneVariableInterface> >()
{
  set_size(nvars, mvars);
}
MatrixAllVariables::MatrixAllVariables( const MatrixAllVariables& matrixallvariables)
{
  assert(0);
}

MatrixAllVariables& MatrixAllVariables::operator=( const MatrixAllVariables& matrixallvariables)
{
  assert(0);
  return *this;
}

std::string MatrixAllVariables::getClassName() const
{
  return "MatrixAllVariables";
}

MatrixAllVariables* MatrixAllVariables::clone() const
{
  assert(0);
  return NULL;
}
const alat::MatrixOneVariableInterface* MatrixAllVariables::get(int i, int j) const
{
  return (*this)(i,j).get();
}
alat::MatrixOneVariableInterface* MatrixAllVariables::get(int i, int j)
{
  return (*this)(i,j).get();
}

/*--------------------------------------------------------------------------*/
std::ostream& alat::operator<<(std::ostream& os, const MatrixAllVariables& matrix)
{
  matrix.save(os, arma::arma_ascii);
  return os;
}

/*--------------------------------------------------------------------------*/
void MatrixAllVariables::save(std::ostream& os, arma::file_type datatype) const
{
  for(int i=0; i< n(); i++)
  {
    for(int j=0;j< m(); j++)
    {
      os << i << "-" << j << ":\n";
      (*this)(i,j)->save(os, datatype);
    }
  }
}

/*--------------------------------------------------------------------------*/
void MatrixAllVariables::fillzeros()
{
  for(MatrixAllVariables::iterator p = this->begin(); p != this->end(); p++)
  {
    (*p)->fillzeros();
  }
}

/*--------------------------------------------------------------------------*/
void MatrixAllVariables::addMatrix(const MatrixAllVariables& matrix, double d)
{
  assert(n()==matrix.n());
  assert(m()==matrix.m());
  for(int i=0; i< n(); i++)
  {
    for(int j=0;j< m(); j++)
    {
      if((*this)(i,j))
      {
        assert((matrix(i,j)));
        (*this)(i,j)->addMatrix(matrix.get(i,j), d);
      }
    }
  }
}

/*--------------------------------------------------------------------------*/
void MatrixAllVariables::matrixVectorProduct(alat::VectorAllVariables* out, const alat::VectorAllVariables* in, double d) const
{
  assert(n()==out->size());
  assert(m()==in->size());
  for(int i=0; i< n(); i++)
  {
    for(int j=0;j< m(); j++)
    {
      if((*this)(i,j)) (*this)(i,j)->matrixVectorProduct(out->get(i), in->get(j), d);
    }
  }
}

/*--------------------------------------------------------------------------*/
void MatrixAllVariables::reInit(MatrixInOne& matrixinone) const
{
  assert(n()==m());
  int nvars = n();
  matrixinone.offsets.set_size(nvars+1);
  matrixinone.offsets[0] = 0;
  for(int i = 0; i < nvars; i++)
  {
    int ni=-1;
    for(int j = 0; j < nvars; j++)
    {
      if((*this)(i,j))
      {
        int n2 = (*this)(i,j)->getSparsityPattern()->n();
        if(ni!=-1){assert(ni==n2);}
        else{ni = n2;}
      }
    }
    matrixinone.offsets[i+1] = matrixinone.offsets[i] + ni;
  }
  int ntotal = matrixinone.offsets[nvars];
  matrixinone.out.set_size(ntotal);
  matrixinone.in.set_size(ntotal);

  alat::SparsityPatternSoft sparsitypatternsoft;
  sparsitypatternsoft.set_size(ntotal);
  for(int i = 0; i < nvars; i++)
  {
    for(int j = 0; j < nvars; j++)
    {
      if((*this)(i,j)) (*this)(i,j)->addEntriesForDirectSolver(matrixinone.offsets[i], matrixinone.offsets[j], sparsitypatternsoft);
    }
  }
  matrixinone.matrix.initSparsityPattern(sparsitypatternsoft);
}

/*--------------------------------------------------------------------------*/
void MatrixAllVariables::compute(MatrixInOne& matrixinone) const
{
  int nvars = n();
  const alat::SparsityPattern* sparsitypattern = matrixinone.matrix.getSparsityPattern();
  arma::vec&  values = *matrixinone.matrix.getValues();
  values.fill(arma::fill::zeros);
  for(int i = 0; i < nvars; i++)
  {
    for(int j = 0; j < nvars; j++)
    {
      if((*this)(i,j)) (*this)(i,j)->addMatrixForDirectSolver(matrixinone.offsets[i], matrixinone.offsets[j], values, sparsitypattern);
    }
  }
}

/*--------------------------------------------------------------------------*/
// void MatrixAllVariables::solve(alat::VectorAllVariables& out, const alat::VectorAllVariables& in)
// {
//   int n = out.size();
//   _offsets.set_size(n+1);
//   _offsets[0] = 0;
//   for(int i = 0; i < n; i++)
//   {
//     _offsets[i+1] = _offsets[i] + out.get(i)->size();
//   }
//   int ntotal = _offsets[n];
//   _out.set_size(ntotal);
//   _in.set_size(ntotal);
//
//   alat::SparsityPatternSoft sparsitypatternsoft;
//   sparsitypatternsoft.set_size(ntotal);
//   for(int i = 0; i < n; i++)
//   {
//     for(int j = 0; j < n; j++)
//     {
//       if((*this)(i,j)) (*this)(i,j)->addEntriesForDirectSolver(_offsets[i], _offsets[j], sparsitypatternsoft);
//     }
//   }
//   _sparsematrix.initSparsityPattern(sparsitypatternsoft);
//
//   const alat::SparsityPattern* sparsitypattern = _sparsematrix.getSparsityPattern();
//   arma::vec&  values = *_sparsematrix.getValues();
//   values.fill(arma::fill::zeros);
//   for(int i = 0; i < n; i++)
//   {
//     for(int j = 0; j < n; j++)
//     {
//       if((*this)(i,j)) (*this)(i,j)->addMatrixForDirectSolver(_offsets[i], _offsets[j], values, sparsitypattern);
//     }
//   }
//   _umf.reInit(&_sparsematrix);
//   _umf.computeLu();
//   _in.fill(arma::fill::zeros);
//   for(int i = 0; i < n; i++)
//   {
//     in.get(i)->addVectorRhsForDirectSolver(_offsets[i], _in);
//   }
//   // std::cerr << "MatrixAllVariables::solve() in=" << in<<"\n";
//   // std::cerr << "MatrixAllVariables::solve() _in=" << _in.t();
//   _umf.solve(_out, _in);
//   // std::cerr << "MatrixAllVariables::solve() _out=" << _out.t();
//   for(int i = 0; i < n; i++)
//   {
//     out.get(i)->setVectorFromDirectSolver(_offsets[i], _out);
//   }
//   // std::cerr << "MatrixAllVariables::solve() out=" << out<<"\n";
// }
