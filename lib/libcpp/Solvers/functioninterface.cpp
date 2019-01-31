#include  "Solvers/functioninterface.hpp"
#include  <cassert>
#include  <iostream>
#include  <stdlib.h>

using namespace solvers;

/*--------------------------------------------------------------------------*/
FunctionInterface::~FunctionInterface(){}
FunctionInterface::FunctionInterface() : _eps(1.e-6), _dim(-1), _ncomp(-1){}
FunctionInterface::FunctionInterface( const FunctionInterface& function)
{
  _ncomp = function._ncomp;
  _dim = function._dim;
  _eps = function._eps;
  assert(0);
}
FunctionInterface& FunctionInterface::operator=( const FunctionInterface& function)
{
  _ncomp = function._ncomp;
  _dim = function._dim;
  _eps = function._eps;
  assert(0);
  return *this;
}
std::string FunctionInterface::getClassName() const
{
  return "FunctionInterface";
}
void FunctionInterface::setNcompAndDim(int ncomp, int dim)
{
  _ncomp=ncomp;
  _dim=dim;
  _ul.set_size(ncomp);
  _ur.set_size(ncomp);
  _epsinv = 0.5/_eps;
}


/*----------------------------------------------------------*/
void FunctionInterface::d (alat::armavec& u, int i, double xv, double yv, double zv, double t) const
{
  if(i == 0){x(u, xv, yv, zv, t);}
  else if(i == 1){y(u, xv, yv, zv, t);}
  else{z(u, xv, yv, zv, t);}
}
void FunctionInterface::dd(alat::armavec& u, int i, int j, double xv, double yv, double zv, double t) const
{
  if(i == 0)
  {
    if(j == 0){xx(u, xv, yv, zv, t);}
    else if(j == 1){xy(u, xv, yv, zv, t);}
    else {xz(u, xv, yv, zv, t);}
  }
  else if(i == 1)
  {
    if(j == 0){xy(u, xv, yv, zv, t);}
    else if(j == 1){yy(u, xv, yv, zv, t);}
    else {yz(u, xv, yv, zv, t);}
  }
  else
  {
    if(j == 0){xz(u, xv, yv, zv, t);}
    else if(j == 1){yz(u, xv, yv, zv, t);}
    else {zz(u, xv, yv, zv, t);}
  }
}
void FunctionInterface::x (alat::armavec& u, double xv, double yv, double zv, double t) const
{
  ( *this )( _ur, xv+_eps, yv, zv, t );
  ( *this )( _ul, xv-_eps, yv, zv, t );
  u = _epsinv*(_ur-_ul);
}
void FunctionInterface::y (alat::armavec& u, double xv, double yv, double zv, double t) const
{
  ( *this )( _ur, xv, yv+_eps, zv, t );
  ( *this )( _ul, xv, yv-_eps, zv, t );
  u = _epsinv*(_ur-_ul);
}
void FunctionInterface::z (alat::armavec& u, double xv, double yv, double zv, double t) const
{
  ( *this )( _ur, xv, yv, zv+_eps, t );
  ( *this )( _ul, xv, yv, zv-_eps, t );
  u = _epsinv*(_ur-_ul);
}
void FunctionInterface::xx(alat::armavec& u, double xv, double yv, double zv, double t) const
{
  x( _ur, xv+_eps, yv, zv, t );
  x( _ul, xv-_eps, yv, zv, t );
  u = _epsinv*(_ur-_ul);
}
void FunctionInterface::yy(alat::armavec& u, double xv, double yv, double zv, double t) const
{
  y( _ur, xv, yv+_eps, zv, t );
  y( _ul, xv, yv-_eps, zv, t );
  u = _epsinv*(_ur-_ul);
}
void FunctionInterface::zz(alat::armavec& u, double xv, double yv, double zv, double t) const
{
  z( _ur, xv, yv, zv+_eps, t );
  z( _ul, xv, yv, zv-_eps, t );
  u = _epsinv*(_ur-_ul);
}
void FunctionInterface::xy(alat::armavec& u, double xv, double yv, double zv, double t) const
{
  x( _ur, xv, yv+_eps, zv, t );
  x( _ul, xv, yv-_eps, zv, t );
  u = _epsinv*(_ur-_ul);
}
void FunctionInterface::yz(alat::armavec& u, double xv, double yv, double zv, double t) const
{
  y( _ur, xv, yv, zv+_eps, t );
  y( _ul, xv, yv, zv-_eps, t );
  u = _epsinv*(_ur-_ul);
}
void FunctionInterface::xz(alat::armavec& u, double xv, double yv, double zv, double t) const
{
  x( _ur, xv, yv, zv+_eps, t );
  x( _ul, xv, yv, zv-_eps, t );
  u = _epsinv*(_ur-_ul);
}
void FunctionInterface::t (alat::armavec& u, double xv, double yv, double zv, double t) const
{
  ( *this )( _ur, xv, yv, zv, t+_eps );
  ( *this )( _ul, xv, yv, zv, t-_eps );
  u = _epsinv*(_ur-_ul);
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
RightHandSideConstant::RightHandSideConstant(double d) : solvers::RightHandSideInterface(), _d(d){}
std::string RightHandSideConstant::getClassName()const {return "RightHandSideConstant";}
void RightHandSideConstant::operator()(alat::armavec& f, double x, double y, double z, double t)const{f.fill(_d);}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
DirichletConstant::DirichletConstant(double d) : solvers::DirichletInterface(), _d(d) {}
std::string DirichletConstant::getClassName()const {return "DirichletConstant";}
void DirichletConstant::operator()(alat::armavec& u, double x, double y, double z, double t)const{u.fill(_d);}
DirichletExactSolution::DirichletExactSolution(const solvers::FunctionInterface& solution_in) : solvers::DirichletInterface(), solution(solution_in){}
std::string DirichletExactSolution::getClassName()const {return "DirichletExactSolution";}
void DirichletExactSolution::operator()(alat::armavec& u, double x, double y, double z, double t)const
{
  // std::cerr << "DirichletExactSolution::operator() u="<<u;
  solution(u, x, y, z, t);
  // std::cerr << "DirichletExactSolution::operator() u="<<u;
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
std::string NeumannZero::getClassName()const {return "NeumannZero";}
void NeumannZero::operator()(alat::armavec& u, double nx, double ny, double nz, double x, double y, double z, double t)const{u.zeros();}
