#include  "Perulangan/iterativesolvervisitorinterface.hpp"
#include  "Perulangan/preconditionerinterface.hpp"
#include  "Perulangan/richardsonoptimal.hpp"
#include  <cassert>
#include  <sstream>

using namespace perulangan;

/*--------------------------------------------------------------------------*/

RichardsonOptimal::~RichardsonOptimal()
{}

RichardsonOptimal::RichardsonOptimal(const std::string& type, int nvectors, const std::string& solutiontype) : RichardsonRB(type, nvectors, solutiontype)
{}

RichardsonOptimal::RichardsonOptimal( const RichardsonOptimal& richardsonoptimal) : RichardsonRB(richardsonoptimal)
{
  assert(0);
}

RichardsonOptimal& RichardsonOptimal::operator=( const RichardsonOptimal& richardsonoptimal)
{
  RichardsonRB::operator=(richardsonoptimal);
  assert(0);
  return *this;
}

std::string RichardsonOptimal::getClassName() const
{
  std::stringstream ss;
  ss<<"RichardsonOptimal_"<< _type << "_" << _nvectors<<"_"<<_solutiontype;
  return ss.str();
}

/*--------------------------------------------------------------------------*/
int RichardsonOptimal::getNVectors() const
{
  return 2*_nvectors + 2;
}

/*--------------------------------------------------------------------------*/
void RichardsonOptimal::_matrixVectorProduct(int index) const
{
  getVisitor()->vectorZero( getAV(index) );
  getVisitor()->matrixVectorProduct(*_ghostmatrix,  getAV(index), getV(index), 1.0);
}

/*--------------------------------------------------------------------------*/
void RichardsonOptimal::_computeSmallSystem(int index, int nmemory) const
{
  _H.resize(nmemory, nmemory);
  _b.resize(nmemory);
  _x.resize(nmemory);
  alat::GhostVector& r = getMemory(2*_nvectors);


  // std::cerr << "rnorm = " <<getVisitor()->vectorNorm(r) << "\n";

  //-------------------------------------------------
  if(_solutiontype == "gal")
  //-------------------------------------------------
  {
    for(int i = 0; i < nmemory; i++)
    {
      _H(index, i) = getVisitor()->vectorDot( getAV(i), getV(index) );
      if(index != i)
      {
        _H(i, index) = getVisitor()->vectorDot( getAV(index), getV(i) );
      }
      _b[i] = getVisitor()->vectorDot( r, getV(i) );
    }
  }
  //-------------------------------------------------
  else if(_solutiontype == "ls")
  //-------------------------------------------------
  {
    for(int i = 0; i < nmemory; i++)
    {
      double d = getVisitor()->vectorDot( getAV(i), getAV(index) );
      _H(index, i) = d;
      if(index != i)
      {
        _H(i, index) = d;
      }
      _b[i] = getVisitor()->vectorDot( r, getAV(i) );
    }
  }
  //-------------------------------------------------
  else if(_solutiontype == "gals")
  //-------------------------------------------------
  {
    double alpha = 0.01;
    for(int i = 0; i < nmemory; i++)
    {
      //ls
      double d = getVisitor()->vectorDot( getAV(i), getAV(index) );
      _H(index, i) = alpha*d;
      if(index != i)
      {
        _H(i, index) = alpha*d;
      }
      //gal
      _H(index, i) += getVisitor()->vectorDot( getAV(i), getV(index) );
      if(index != i)
      {
        _H(i, index) += getVisitor()->vectorDot( getAV(index), getV(i) );
      }
      //ls
      _b[i] = alpha*getVisitor()->vectorDot( r, getAV(i) );
      //gal
      _b[i] += getVisitor()->vectorDot( r, getV(i) );
    }
  }
  // std::cerr << "nmemory="<<nmemory<<" index="<<index<<"\n";
  // std::cerr << "_H="<<_H<<" _b="<<_b<<"\n";
//  int nrang = arma::rank(_H);
//  if(nrang != _H.n_rows)
//  {
//    std::cerr << getVisitor()->getClassName() << " rang(_H)="<<nrang<<" n="<<_H.n_rows<<" diff="<<_H.n_rows-nrang<<"\n";
//  }
}

/*--------------------------------------------------------------------------*/
void RichardsonOptimal::_addvector(alat::GhostVector& u, int nmemory) const
{
  // std::cerr << "_x="<<_x<<"\n";
  // assert(0);
  // getVisitor()->vectorZero( u);
  for(int i = 0; i < nmemory; i++)
  {
    getVisitor()->vectorAdd( u, _rnorm*_x[i], getV(i) );
  }
}
