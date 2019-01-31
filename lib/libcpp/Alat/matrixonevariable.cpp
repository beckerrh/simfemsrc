#include  "Alat/matrixonevariable.hpp"
#include  "Alat/sparsitypatternsoft.hpp"
#include  "Alat/vectoronevariable.hpp"
#include  <cassert>
#include  <iomanip>

using namespace alat;

/*--------------------------------------------------------------------------*/
MatrixOneVariable::~MatrixOneVariable() {}
MatrixOneVariable::MatrixOneVariable(): MatrixOneVariableInterface(){}
MatrixOneVariable::MatrixOneVariable( const MatrixOneVariable& matrixonevariable): MatrixOneVariableInterface(matrixonevariable)
{
  (*this).operator=(matrixonevariable);
}
MatrixOneVariable& MatrixOneVariable::operator=( const MatrixOneVariable& matrixonevariable)
{
  assert(0);
  MatrixOneVariableInterface::operator=(matrixonevariable);
  return *this;
}
std::string MatrixOneVariable::getClassName() const
{
  return "MatrixOneVariable";
}
MatrixOneVariable* MatrixOneVariable::clone() const
{
  return new MatrixOneVariable(*this);
}
/*--------------------------------------------------------------------------*/
const alat::SparsityPattern* MatrixOneVariable::getSparsityPattern() const
{return &_sparsitypattern;}
const alat::armavec*  MatrixOneVariable::getValues() const {return &_values;}
alat::armavec*  MatrixOneVariable::getValues() {return &_values;}
bool MatrixOneVariable::needsConnectivity() const {return true;}
void MatrixOneVariable::initSparsityPattern(const SparsityPatternSoft& sparsitypatternsoft)
{
  // std::cerr << "MatrixOneVariable::initSparsityPattern()" << sparsitypatternsoft << "\n";
  _sparsitypattern.set_size(sparsitypatternsoft);
  _values.set_size(sparsitypatternsoft.ntotal());
}
void MatrixOneVariable::set_size(int n, int m) {}

/*--------------------------------------------------------------------------*/
void MatrixOneVariable::save(std::ostream& os, arma::file_type datatype) const
{
  os<<"sparsitypattern:\n";
  _sparsitypattern.save(os, datatype);
  os<<"values:\n";
  _values.save(os, datatype);
}
void MatrixOneVariable::write(std::ostream& os) const
{
  os << *this;
}

/*--------------------------------------------------------------------------*/
std::ostream& alat::operator<<(std::ostream& os, const MatrixOneVariable& A)
{
  const alat::SparsityPattern* sparsitypattern = A.getSparsityPattern();
  const alat::armavec& values =  *A.getValues();
  int n = sparsitypattern->n();
  assert(values.size()==sparsitypattern->ntotal());
  os << "n = " << n << "\n";
  for(int i = 0; i < n; i++)
  {
    os << "rowlength = " << sparsitypattern->rowsize(i) << "\n";
    for(int pos = sparsitypattern->rowstart(i); pos < sparsitypattern->rowstop(i); pos++)
    {
      int j = sparsitypattern->col(pos);
      os << i << "," << j << " (" << std::setprecision(4) << values[pos] << ") ";
    }
    os << "\n";
  }
  return os;
}

/*--------------------------------------------------------------------------*/
void MatrixOneVariable::fillzeros()
{
  _values.fill(arma::fill::zeros);
}
/*--------------------------------------------------------------------------*/
void MatrixOneVariable::matrixVectorProduct(alat::VectorOneVariableInterface* out, const alat::VectorOneVariableInterface* in, double d) const
{
  alat::VectorOneVariable* outv = dynamic_cast<alat::VectorOneVariable*>(out); assert(outv);
  const alat::VectorOneVariable* inv = dynamic_cast<const alat::VectorOneVariable*>(in); assert(inv);
  for(int index=0;index<_sparsitypattern.n(); index++)
  {
    for(int pos = _sparsitypattern.rowstart(index); pos < _sparsitypattern.rowstop(index); pos++)
    {
      (*outv)[index] += d * _values[pos] * (*inv)[_sparsitypattern.col(pos)];
    }
  }
}
/*--------------------------------------------------------------------------*/
void MatrixOneVariable::addMatrix(const MatrixOneVariableInterface* matrix, double d)
{
  assert(0);
}
/*--------------------------------------------------------------------------*/
void MatrixOneVariable::assemble(const alat::armamat& Alocal, const alat::armaivec& indicesi, const alat::armaivec& indicesj)
{
  // std::cerr << "Alocal = " << Alocal.t();
  // std::cerr << "indicesi = " << indicesi.t();
  // std::cerr << "indicesj = " << indicesj.t();
  for(int ii=0;ii<indicesi.size();ii++)
  {
    int i = indicesi[ii];
    for(int jj=0;jj<indicesj.size();jj++)
    {
      int j = indicesj[jj];
      double val = Alocal(ii,jj);
      int posstart = _sparsitypattern.rowstart(i);
      int posend = _sparsitypattern.rowstop(i);
      bool found = false;
      for(int pos = posstart; pos < posend; pos++)
      {
        if(_sparsitypattern.col(pos) == j)
        {
          _values[pos] += val;
          found=true;
          break;
        }
      }
      assert(found);
    }
  }
}

// void MatrixOneVariable::assemble(const alat::armaivec& indicesi, const alat::armaivec& indicesj, const arma::mat& Alocal)
// {
//   int nloci = indicesi.size();
//   int nlocj = indicesj.size();
//   for(int ii = 0; ii < nloci; ii++)
//   {
//     int i = indicesi[ii];
//     int posstart = _sparsitypattern.rowstart(i);
//     int posend = _sparsitypattern.rowstop(i);
//     for(int jj = 0; jj < nlocj; jj++)
//     {
//       int j = indicesj[jj];
//       for(int pos = posstart; pos < posend; pos++)
//       {
//         if(_sparsitypattern.col(pos) == j)
//         {
//           _values[pos] += Alocal(ii,jj);
//           break;
//         }
//       }
//     }
//   }
// }
/*--------------------------------------------------------------------------*/
void MatrixOneVariable::rowIdentity(int index)
{
  int posstart = _sparsitypattern.rowstart(index);
  int posend = _sparsitypattern.rowstop(index);
  for(int pos = posstart; pos < posend; pos++)
  {
    if(_sparsitypattern.col(pos) == index)
    {
      _values[pos] = 1.0;
    }
    else
    {
      _values[pos] = 0.0;
    }
  }
}
/*--------------------------------------------------------------------------*/
void MatrixOneVariable::rowZero(int index)
{
  int posstart = _sparsitypattern.rowstart(index);
  int posend = _sparsitypattern.rowstop(index);
  for(int pos = posstart; pos < posend; pos++)
  {
    _values[pos] = 0.0;
  }
}
/*--------------------------------------------------------------------------*/
void MatrixOneVariable::solve(VectorOneVariableInterface* u, const VectorOneVariableInterface* f)
{
  _umf.reInit(this);
  _umf.computeLu();
  alat::armavec* ud = dynamic_cast<alat::armavec* >(u); assert(ud);
  const alat::armavec* fd = dynamic_cast<const alat::armavec* >(f); assert(fd);
  _umf.solve(*ud, *fd);
  // std::cerr << "u = " << *ud << "\n";
}
/*--------------------------------------------------------------------------*/
void MatrixOneVariable::addEntriesForDirectSolver(int offsetivar, int offsetjvar, alat::SparsityPatternSoft& sparsitypatternsoft) const
{
  int n = _sparsitypattern.n();
  for(int i = 0; i < n; i++)
  {
    for(int pos = _sparsitypattern.rowstart(i); pos < _sparsitypattern.rowstop(i); pos++)
    {
      int j = _sparsitypattern.col(pos);
      sparsitypatternsoft[offsetivar + i].insert(offsetjvar + j);
    }
  }
}
void MatrixOneVariable::addMatrixForDirectSolver(int offsetivar, int offsetjvar, alat::armavec& matrixvalues, const alat::SparsityPattern* sparsitypattern) const
{
  int n = _sparsitypattern.n();
  for(int i = 0; i < n; i++)
  {
    for(int pos = _sparsitypattern.rowstart(i); pos < _sparsitypattern.rowstop(i); pos++)
    {
      int j = _sparsitypattern.col(pos);
      for(int pos2 = sparsitypattern->rowstart(offsetivar+i); pos2 < sparsitypattern->rowstop(offsetivar+i); pos2++)
      {
        int j2 = sparsitypattern->col(pos2);
        if(j2 == j+offsetjvar)
        {
          matrixvalues[pos2] +=  _values[pos];
        }
      }
    }
  }
}
