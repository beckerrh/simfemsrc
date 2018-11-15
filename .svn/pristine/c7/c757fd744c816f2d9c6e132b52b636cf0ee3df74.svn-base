#include  "Alat/stringvector.hpp"
#include  "Alat/vectoronevariable.hpp"
#include  <cassert>
#include  <cmath>
#include  <iomanip>
#include  <iostream>
#include  <limits>
#include  <sstream>

using namespace alat;

/*----------------------------------------------------------*/
VectorOneVariable::~VectorOneVariable() {}
VectorOneVariable::VectorOneVariable() : arma::vec(), alat::VectorOneVariableInterface(), _ncomp(1), _n(-1){}
VectorOneVariable::VectorOneVariable(int ncomp) : arma::vec(), alat::VectorOneVariableInterface(), _ncomp(ncomp), _n(-1){}
VectorOneVariable::VectorOneVariable(const VectorOneVariable& v) : arma::vec(v), alat::VectorOneVariableInterface(v), _ncomp(v._ncomp), _n(v._n){}

VectorOneVariable& VectorOneVariable::operator=(const VectorOneVariable& v)
{
  assert(0);
  return *this;
}
std::unique_ptr<VectorOneVariableInterface> VectorOneVariable::clone() const
{
  return std::unique_ptr<VectorOneVariableInterface>(new VectorOneVariable(*this));
}
std::string VectorOneVariable::getClassName() const
{
  return "VectorOneVariable";
}
/*----------------------------------------------------------*/
std::ostream& alat::operator<<(std::ostream& os, const VectorOneVariable& g)
{
  const arma::vec& tarma =static_cast<const arma::vec&>(g);
  os << tarma.t();
  return os;
}

/*----------------------------------------------------------*/
void VectorOneVariable::copyFrom(const VectorOneVariable& u)
{
  assert(_ncomp!=-1);
  for(int icomp=0;icomp<_ncomp;icomp++)
  {
    VectorOneVariable::const_iterator q = u.begin(icomp);
    for(iterator p=begin(icomp);p!=end(icomp);p++)
    {
      *p= *q++;
    }
  }
}

/*----------------------------------------------------------*/
void VectorOneVariable::scaleIntVector(const alat::armaivec& count)
{
  assert(_ncomp!=-1);
  for(int icomp=0;icomp<_ncomp;icomp++)
  {
    alat::armaivec::const_iterator q = count.begin();
    for(iterator p=begin(icomp);p!=end(icomp);p++)
    {
      *p /= (double)*q++;
    }
  }
}

/*----------------------------------------------------------*/
void VectorOneVariable::loadhdf5(const arma::hdf5_name& spec)
{
  arma::vec& tarma =static_cast<arma::vec&>(*this);
  tarma.load(spec);
}
void VectorOneVariable::savehdf5(const arma::hdf5_name& spec) const
{
  const arma::vec& tarma =static_cast<const arma::vec&>(*this);
  tarma.save(spec);
}
/*----------------------------------------------------------*/
void VectorOneVariable::save(std::ostream& os, arma::file_type datatype) const
{
  const arma::vec& tarma =static_cast<const arma::vec&>(*this);
  tarma.save(os, datatype);
}
/*----------------------------------------------------------*/
int VectorOneVariable::size() const {return arma::vec::size();}
int VectorOneVariable::ncomp() const {return _ncomp;}
void VectorOneVariable::setNcompAndN(int ncomp, int n) {_ncomp=ncomp; _n=n;}
void VectorOneVariable::set_size(int n)
{
  _n = n;
  arma::vec::set_size( n );
}
/*----------------------------------------------------------*/
VectorOneVariable::const_iterator VectorOneVariable::begin(int icomp) const{return arma::vec::begin() + icomp* _n;}
VectorOneVariable::iterator VectorOneVariable::begin(int icomp){return arma::vec::begin() + icomp* _n;}
VectorOneVariable::iterator VectorOneVariable::end(int icomp){return arma::vec::begin() + ( icomp+1 )* _n;}
VectorOneVariable::const_iterator VectorOneVariable::end(int icomp) const{return arma::vec::begin() + ( icomp+1 )* _n;}
int VectorOneVariable::n() const{return _n;}
/*----------------------------------------------------------*/
void VectorOneVariable::assemble(const alat::armaimat& indices, const arma::mat& local, double d)
{
  for(int ii = 0; ii < indices.n_rows; ii++)
  {
    for(int jj = 0; jj < indices.n_cols; jj++)
    {
      int i = indices(ii,jj);
      ( *this )[i] += d*local(ii,jj);
    }
  }
}
void VectorOneVariable::extract(const alat::armaimat& indices, arma::mat& local) const
{
  local.set_size(indices.n_rows, indices.n_cols);
  for(int ii = 0; ii < indices.n_rows; ii++)
  {
    for(int jj = 0; jj < indices.n_cols; jj++)
    {
      int i = indices(ii,jj);
      local(ii,jj) = ( *this )[i];
    }
  }
}

/*----------------------------------------------------------*/
void VectorOneVariable::scale(double s)
{
  arma::vec::operator*=(s);
}

/*----------------------------------------------------------*/
void VectorOneVariable::fillzeros()
{
  arma::vec::fill(arma::fill::zeros);
}

/*----------------------------------------------------------*/
double VectorOneVariable::norm() const
{
  const arma::vec& tarma =static_cast<const arma::vec&>(*this);
  return arma::norm(tarma);
}

/*----------------------------------------------------------*/
double VectorOneVariable::dot(const alat::VectorOneVariableInterface* v) const
{
  const arma::vec& tarma =static_cast<const arma::vec&>(*this);
  const arma::vec& varma =dynamic_cast<const arma::vec&>(*v);
  return arma::dot( tarma, varma );
}

/*----------------------------------------------------------*/
void VectorOneVariable::equal(double d)
{
  // assert(0);
  arma::vec& tarma =static_cast<arma::vec&>(*this);
  tarma.fill(d);
}
void VectorOneVariable::equal(const alat::VectorOneVariableInterface* v)
{
  arma::vec& tarma =static_cast<arma::vec&>(*this);
  const arma::vec& varma =dynamic_cast<const arma::vec&>(*v);
  tarma = varma;
}

/*----------------------------------------------------------*/
void VectorOneVariable::add(const double& d, const alat::VectorOneVariableInterface* v)
{
  arma::vec& tarma =static_cast<arma::vec&>(*this);
  const arma::vec& varma =dynamic_cast<const arma::vec&>(*v);
  tarma += d*varma;
}

/*----------------------------------------------------------*/
void VectorOneVariable::setVectorFromDirectSolver(int offset, const arma::vec& u)
{
  arma::vec& tarma =static_cast<arma::vec&>(*this);
  tarma = u.subvec(offset, offset+tarma.size()-1);
}
void VectorOneVariable::addVectorRhsForDirectSolver(int offset, arma::vec& f) const
{
  const arma::vec& tarma =static_cast<const arma::vec&>(*this);
  // std::cerr << "tarma " << tarma.t();
  // std::cerr << "offset " << offset << "\n";
  // std::cerr << "offset+tarma.size()-1 " << offset+tarma.size()-1 << "\n";
  f.subvec(offset, offset+tarma.size()-1) += tarma;
}
