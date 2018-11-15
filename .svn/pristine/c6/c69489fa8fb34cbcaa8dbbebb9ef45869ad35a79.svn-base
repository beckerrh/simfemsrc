#include  "Perulangan/iterativesolvervisitorinterface.hpp"
#include  "Perulangan/preconditionerinterface.hpp"
#include  "Perulangan/gmres.hpp"
#include  <cassert>
#include  <sstream>

using namespace perulangan;

/*--------------------------------------------------------------------------*/
Gmres::~Gmres(){}
Gmres::Gmres(const std::string& type, int nvectors, const std::string& solutiontype) : RichardsonOptimal(type, nvectors, solutiontype){}
Gmres::Gmres( const Gmres& richardsonoptimal) : RichardsonOptimal(richardsonoptimal)
{
  assert(0);
}
Gmres& Gmres::operator=( const Gmres& richardsonoptimal)
{
  RichardsonOptimal::operator=(richardsonoptimal);
  assert(0);
  return *this;
}

std::string Gmres::getClassName() const
{
  std::stringstream ss;
  ss<<"Gmres_"<< _type << "_" << _nvectors<<"_"<<_solutiontype;
  return ss.str();
}

/*--------------------------------------------------------------------------*/
int Gmres::getNVectors() const
{
  return 2*_nvectors + 2;
}

/*--------------------------------------------------------------------------*/
void Gmres::solve(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f) const
{
  _nmemory = 0;
  _nextmemory = 0;
  RichardsonOptimal::solve(status, A, u, f);
}
